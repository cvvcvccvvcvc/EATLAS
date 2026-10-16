from __future__ import annotations

import copy
import gzip
import json
from pathlib import Path

import duckdb
import pytest

from analytics.io.reannotation import prepare_repaired_source
from analytics.io.run_inputs import resolve_source_runs
from analytics.io.variant_source import resolve_variant_table_source, variant_source_sql, variant_sources_by_partition
from genomics.gnomad_cache import GnomadRegionCache
from provenance.evidence_inventory import EVIDENCE_SCOPES, inventory_file_descriptor, write_evidence_inventory
from reannotation.archive import INPUT_FILES, download_snapshot, source_from_archive, upload_snapshot
from reannotation.runner import run_reannotation
from reannotation.store import digest, empty_state, load_snapshot, publish, write_fragment
from run_archiving.archive import archive_run
from test_analysis_inputs import _make_source_run, _table
from test_run_archiving import LocalRemote

WINDOW = {"started_at_utc": "2026-10-16T00:00:00+00:00", "finished_at_utc": "2026-10-16T00:00:01+00:00"}
PROVENANCE = {"git_commit": "b" * 40, "git_dirty": False}


def record(pos=2, af=0.01):
    return {"chrom": "1", "pos": pos, "ref": "A", "alt": "G", "variant_id": f"1-{pos}-A-G",
            "exome": {"af": af}, "genome": None, "joint": None, "consequence": "missense_variant"}


class TrackingRemote(LocalRemote):
    def __init__(self, root):
        super().__init__(root)
        self.downloads = []
        self.uploads = []

    def download_file(self, source, destination):
        self.downloads.append(source)
        return super().download_file(source, destination)

    def copy_file(self, source, destination):
        self.uploads.append(destination)
        return super().copy_file(source, destination)


@pytest.fixture
def archived(tmp_path):
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"reference")
    Path(str(clinvar) + ".tbi").write_bytes(b"index")
    run = _make_source_run(tmp_path / "batch_001", gene_id="1", clinvar_vcf=clinvar, variant_key="1:2:A>G")
    with gzip.open(run / "fetch/sequences/targets/1.fa.gz", "wt") as handle:
        handle.write(">1\nAAAAAAAAAA\n")
    _table(run / "annotation/failures.tsv.gz", ["source", "scope", "chrom", "start", "end"],
           [{"source": "gnomad", "scope": "region", "chrom": "1", "start": 1, "end": 50000}])
    inventory = run / "evidence_inventory.json"
    write_evidence_inventory(inventory, {scope: run / scope for scope in EVIDENCE_SCOPES})
    root = json.loads((run / "run_manifest.json").read_text())
    root.update(exit_status=0, completed_at="2026-10-16T00:00:00Z")
    root["evidence_inventory"] = inventory_file_descriptor(inventory)
    (run / "run_manifest.json").write_text(json.dumps(root))
    remote = TrackingRemote(tmp_path / "remote")
    archive_run(remote, run_dir=run, remote_root="drive:GAPH")
    remote.uploads.clear()
    workspace = tmp_path / "metadata"
    workspace.mkdir()
    binding = source_from_archive(remote, "drive:GAPH", "batch_001", workspace / "batch_001", 1024**3)
    source = resolve_source_runs([run], clinvar_vcf=clinvar)[0]
    return remote, binding, source


def make_snapshot(root, binding, observations):
    state = empty_state()
    state["sources"]["batch_001"] = binding
    if observations:
        write_fragment(root, observations, state)
    return state, publish(root, state, PROVENANCE)


def repaired_row(source, snapshot, tmp_path):
    repaired = prepare_repaired_source(source, snapshot, tmp_path / "analytics")
    variants = resolve_variant_table_source(repaired.variant_annotations_source, required_columns={"variant_key"})
    with duckdb.connect() as con:
        row = con.sql(f"SELECT gnomad_af,gnomad_lookup_status FROM {variant_source_sql(variants)}").fetchone()
    return row, variants


def test_archive_planning_downloads_only_small_inputs(archived):
    remote, binding, source = archived
    assert {x.split('/data/')[1] for x in remote.downloads if '/data/' in x} == set(INPUT_FILES)
    assert binding["tiles"] == ["1:1:25000", "1:25001:50000", "1:50001:75000"]
    assert source.run_dir.exists()


def test_three_versions_overlay_without_copying_or_changing_sources(tmp_path, archived):
    remote, binding, source = archived
    before = {str(p): digest(p) for p in source.run_dir.rglob('*') if p.is_file()}
    root = tmp_path / "supplement"
    state, v1 = make_snapshot(root, binding, {})
    assert repaired_row(source, v1, tmp_path)[0] == ("", "")
    write_fragment(root, {binding["tiles"][0]: ([record()], WINDOW)}, state)
    v2 = publish(root, state, PROVENANCE)
    row, table2 = repaired_row(source, v2, tmp_path)
    assert row == ("0.01", "found")
    first_objects = set(load_snapshot(v2)["objects"])
    write_fragment(root, {key: ([], WINDOW) for key in binding["tiles"][1:]}, state)
    v3 = publish(root, state, PROVENANCE)
    assert load_snapshot(v3)["complete"]
    assert first_objects <= load_snapshot(v3)["objects"].keys()
    assert len(list((root / "objects").glob('*.parquet'))) == 2
    assert repaired_row(source, v1, tmp_path)[0] == ("", "")
    assert repaired_row(source, v2, tmp_path)[0] == ("0.01", "found")
    assert repaired_row(source, v3, tmp_path)[0] == ("0.01", "found")
    assert variant_sources_by_partition(table2, ["partition_000001"])["partition_000001"].identity["supplement"] == v2.stem
    assert before == {str(p): digest(p) for p in source.run_dir.rglob('*') if p.is_file()}
    for snap in (v1, v2, v3):
        upload_snapshot(remote, snap, "drive:GAPH/reannotation")
    assert len([x for x in remote.uploads if '/objects/' in x]) == 2
    assert '/snapshots/' in remote.uploads[-1]
    downloaded = download_snapshot(remote, "drive:GAPH/reannotation", v3.stem, tmp_path / 'download')
    assert load_snapshot(downloaded) == load_snapshot(v3)


def test_empty_response_only_resolves_fully_covered_region(tmp_path, archived):
    _, binding, source = archived
    state, partial = make_snapshot(tmp_path / 'store', binding, {binding['tiles'][0]: ([], WINDOW)})
    assert repaired_row(source, partial, tmp_path)[0] == ('', '')
    write_fragment(tmp_path / 'store', {key: ([], WINDOW) for key in binding['tiles'][1:]}, state)
    complete = publish(tmp_path / 'store', state, PROVENANCE)
    assert repaired_row(source, complete, tmp_path)[0] == ('', 'not_found')


@pytest.mark.parametrize('observations,expected', [({}, 'lookup_failed'),
    ({'1:1:25000': ([record()], WINDOW)}, 'found'),
    ({key: ([], WINDOW) for key in ['1:1:25000', '1:25001:50000', '1:50001:75000']}, 'not_found')])
def test_summary_and_filtering_match_raw_and_materialized_supplements(tmp_path, archived, observations, expected):
    from analytics.analyses.variant_summary_aggregation import aggregate_variant_groups, resolve_variant_aggregation_source
    from analytics.analyses.basic_filtering import build_or_load_filter_score_store
    from analytics.io.archive_sources import prepare_variant_parquet
    _, binding, source = archived
    _, snap = make_snapshot(tmp_path / 'store', binding, observations)
    repaired = prepare_repaired_source(source, snap, tmp_path / 'views')
    prepared = tmp_path / 'prepared'
    prepare_variant_parquet(repaired, prepared)
    support = tmp_path / 'support.tsv.gz'
    _table(support, ['variant_key', 'gene_id', 'strategy', 'alt_support_ortholog_count', 'alt_support_family_count', 'site_aligned_ortholog_count'],
           [{'variant_key': '1:2:A>G', 'gene_id':'1', 'strategy': 's1', 'alt_support_ortholog_count': 2, 'alt_support_family_count': 1, 'site_aligned_ortholog_count': 3}])
    for index, path in enumerate([source.variant_annotations_source, repaired.variant_annotations_source, prepared / 'manifest.json']):
        status = 'lookup_failed' if index == 0 else expected
        result = aggregate_variant_groups(resolve_variant_aggregation_source(path),
                    annotation_failures_path=source.annotation_failures_tsv, threads=1)
        assert result.global_groups['gnomad_status'].tolist() == [status]
        scores, _, _ = build_or_load_filter_score_store(variant_annotations_source=path,
            variant_strategy_support_tsv=support, annotation_failures_tsv=source.annotation_failures_tsv,
            analytics_dir=tmp_path / f'filter{index}', strategies=['s1'])
        with duckdb.connect() as con:
            assert set(con.read_parquet(str(scores)).df()['gnomad_status']) == {status}


def test_failed_upload_does_not_publish_snapshot_and_resume_reuses_objects(tmp_path, archived, monkeypatch):
    remote, binding, _ = archived
    _, snap = make_snapshot(tmp_path / 'store', binding, {binding['tiles'][0]: ([record()], WINDOW)})
    original = remote.verify_checksum
    def fail(*args):
        raise RuntimeError('Drive temporarily unavailable')
    monkeypatch.setattr(remote, 'verify_checksum', fail)
    with pytest.raises(RuntimeError, match='Drive'):
        upload_snapshot(remote, snap, 'drive:GAPH/supplement')
    assert not any('/snapshots/' in x for x in remote.uploads)
    monkeypatch.setattr(remote, 'verify_checksum', original)
    upload_snapshot(remote, snap, 'drive:GAPH/supplement')
    assert remote.uploads[-1].endswith(snap.name)


def test_snapshot_corruption_and_wrong_binding_are_rejected(tmp_path, archived):
    _, binding, source = archived
    wrong = copy.deepcopy(binding)
    wrong['root_manifest_sha256'] = '0' * 64
    _, bad = make_snapshot(tmp_path / 'wrong', wrong, {})
    with pytest.raises(ValueError, match='bind source'):
        prepare_repaired_source(source, bad, tmp_path / 'analytics')
    _, good = make_snapshot(tmp_path / 'store', binding, {binding['tiles'][0]: ([record()], WINDOW)})
    value = load_snapshot(good)
    obj = good.parent.parent / next(iter(value['objects'].values()))['path']
    obj.write_bytes(b'broken')
    with pytest.raises(ValueError, match='changed'):
        load_snapshot(good)


def test_ortholog_evidence_uses_selected_version_and_invalidates_cache(tmp_path):
    from types import SimpleNamespace
    from analytics.io.annotation_support import build_or_load_annotation_support
    from test_annotation_support_cache import _write_source_contract, _read_rows
    run = tmp_path / 'source'
    contract = _write_source_contract(run)
    genes = run / 'fetch/genes.tsv.gz'
    _table(genes, ['gene_id','genomic_accession','chromosome','begin','end'],
           [{'gene_id':'1','genomic_accession':'NC_000001.11','chromosome':'1','begin':100,'end':129}])
    targets = run / 'fetch/targets'
    targets.mkdir()
    with gzip.open(targets / '1.fa.gz', 'wt') as handle:
        handle.write('>1\nAAACGAAAAAAAAAAAAAAAAAAAAAAAAA\n')
    root_manifest = run / 'run_manifest.json'
    root_manifest.write_text('{}')
    source = SimpleNamespace(run_dir=run, source_id='test', root_manifest_json=root_manifest,
        evidence_inventory={'tree_sha256':'e'*64}, variant_annotations_source=contract['source_annotations'],
        genes_tsv=genes, target_sequences_dir=targets)
    # prepare_repaired_source returns a dataclass replacement; use the real SourceRun's fields via a minimal dataclass.
    from dataclasses import make_dataclass
    source = make_dataclass('FixtureSource', [(key, object) for key in vars(source)])(**vars(source))
    binding = {'run_id':'source','root_manifest_sha256':digest(root_manifest), 'evidence_tree_sha256':'e'*64,
               'regions':[['1',104,104]], 'tiles':['1:1:25000']}
    state, partial = make_snapshot(tmp_path / 'store', binding, {})
    write_fragment(tmp_path / 'store', {'1:1:25000':([],WINDOW)}, state)
    absent = publish(tmp_path / 'store', state, PROVENANCE)
    write_fragment(tmp_path / 'store', {'1:1:25000':([{**record(104), 'ref':'G','alt':'A','variant_id':'1-104-G-A'}],WINDOW)}, state)
    found = publish(tmp_path / 'store', state, PROVENANCE)
    failures = []
    found_counts = []
    for snap in [partial, absent, found, partial]:
        repaired = prepare_repaired_source(source, snap, tmp_path / 'views')
        outputs = build_or_load_annotation_support(partition_dirs=[contract['partition']], map_root=contract['map_root'],
            taxonomy=contract['taxonomy'], target_features=contract['target_features'],
            variant_annotations_source=repaired.variant_annotations_source, failures=contract['failures'],
            alignment_manifest=contract['alignment_manifest'], annotation_manifest=contract['annotation_manifest'],
            analytics_dir=tmp_path / 'analytics')
        rows = _read_rows(outputs.ortholog_evidence_summary_tsv)
        failures.append(sum(int(row['gnomad_lookup_failed_count']) for row in rows))
        found_counts.append(sum(int(row['gnomad_found_count']) for row in rows))
    assert failures[0] > 0 and failures[0] == failures[3]
    assert failures[1:3] == [0,0]
    assert found_counts[2] > found_counts[1]


def test_resume_fetches_only_remaining_tiles_and_refresh_replaces_one(tmp_path, archived, monkeypatch):
    remote, binding, _ = archived
    calls = []
    def fetch(chrom, start, end, **kwargs):
        calls.append(start)
        return [record()] if start == 1 else []
    cache = GnomadRegionCache(tmp_path / 'cache', fetcher=fetch)
    options = dict(root=tmp_path / 'store', archive_root='drive:GAPH', output_remote='drive:GAPH/supplement',
                   run_ids=['batch_001'], client=remote, cache=cache, provenance=PROVENANCE)
    first = run_reannotation(**options, max_fetch_tiles=1)
    assert not load_snapshot(first)['complete']
    second = run_reannotation(**options)
    assert load_snapshot(second)['complete']
    assert calls == [1, 25001, 50001]
    run_reannotation(**options)
    assert calls == [1, 25001, 50001]
    from reannotation import runner
    monkeypatch.setattr(runner, 'GnomadRegionCache', lambda path: GnomadRegionCache(path, fetcher=fetch))
    third = run_reannotation(**options, parent=second, refresh=['1:1-10'])
    assert calls == [1, 25001, 50001, 1]
    assert load_snapshot(third)['complete']
    assert load_snapshot(second)['complete']


def test_pipeline_retries_pending_tile_without_refetching_success(tmp_path, monkeypatch):
    from genomics import gnomad_cache
    calls = []
    def fetch(chrom, start, end, **kwargs):
        calls.append(start)
        if start == 25001 and calls.count(start) == 1:
            raise TimeoutError('temporary')
        return []
    monkeypatch.setattr(gnomad_cache.time, 'sleep', lambda _: None)
    cache = GnomadRegionCache(tmp_path / 'cache', fetcher=fetch, group_attempts=1, max_attempts=1)
    assert cache.fetch_region_complete('1', 1, 50000) == []
    assert calls == [1, 25001, 25001]


def test_shared_pacing_and_retry_after(tmp_path, monkeypatch):
    from genomics import gnomad_transport as transport
    monkeypatch.setenv('GAPH_GNOMAD_RATE_STATE', str(tmp_path / 'pacing.json'))
    clock = [1000.0]
    monkeypatch.setattr(transport.time, 'time', lambda: clock[0])
    monkeypatch.setattr(transport.time, 'sleep', lambda delay: clock.__setitem__(0, clock[0] + delay))
    transport.wait_for_request_slot()
    transport.wait_for_request_slot()
    assert clock[0] == 1007.5
    transport.defer_requests('120')
    transport.wait_for_request_slot()
    assert clock[0] == 1127.5
