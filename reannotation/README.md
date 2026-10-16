# Incremental gnomAD reannotation

`python -m reannotation` produces additional scientific evidence for failed
gnomAD regions in completed, archived runs. It is an auxiliary producer, separate
from the end-to-end pipeline and from read-only analytics. It never modifies or
reuploads a source batch. It uses the gnomAD browser API and the same regional
cache, request pacing, normalization semantics, and AF selection as the pipeline.
The shared pacing stays below the provider's published
[ten requests per IP per minute](https://gnomad.broadinstitute.org/news/2025-08-gnomad-2024-user-survey-results/).

## Storage and reproducibility

The planner downloads only the archive-bound run manifest, evidence inventory,
annotation manifest, and failure table. It verifies their hashes and removes its
temporary copies. It does not restore candidate tables, alignments, or FASTA.
Only recorded failed gnomAD regions are queried, with the annotation boundary's
100 bp padding and shared 25 kb cache tiles. Already cached successes are reused.
An empty successful response is evidence of absence; a failed request is pending.

The output directory contains `objects/<SHA256>.parquet`,
`snapshots/<SHA256>.json`, and disposable operational state under `work/`.
Each immutable snapshot binds the original archive and evidence identities,
code commit, dataset selector, observation times, active tiles, and unresolved
tiles. The selector `gnomad_r4` is not a frozen provider release. Browser data
observed today need not equal historical data. The recorded response fragments
make the actual observations reproducible.

Snapshots contain a flattened map to shared fragments. Version 2 or 3 references
unchanged fragments; analytics never walks a chain of versions. A deliberate
correction replaces the active reference for specified tiles and retains old
fragments for old snapshots. This required history is the only retained overlap;
successful source candidate tables are never duplicated in the supplement.
Objects are uploaded and checked on Drive before the snapshot is published.
Source archives remain untouched. One writer is allowed per local output store.

## Run on ITMO

Create the task environment from the repository definition:

```bash
micromamba create --yes -p "$GAPH_ROOT/envs/reannotation" -f envs/reannotation.yml
```

Submit from the Slurm controller, using a clean checkout at freshly fetched
`origin/main`, following the [revision gate](../docs/pipeline_launch.md):

```bash
bash reannotation/slurm/submit_reannotation.sh --expected-commit "$INTENDED_COMMIT" -- \
  --root "$GAPH_ROOT/reannotation/gnomad" \
  --archive-root gdrive:GAPH --output-remote gdrive:GAPH/reannotation/gnomad \
  --gnomad-cache-dir "$GAPH_ROOT/cache/gnomad" \
  --run-id batch_001 --run-id batch_002
```

Repeat `--run-id` for the selected batches. The CLI accepts the existing rclone
configuration; a non-default configuration is a global option before `run`.
Only metadata is staged, sequentially. Successful observations are flushed in
bounded fragments and published after each batch. The local supplement has a
20 GiB default budget and a 1 GiB free-space reserve. The shared regional cache
is additional storage and must fit the user's quota too; filesystem free space
does not measure an account quota. No unrelated files are deleted automatically.

Repeat the same command after interruption: committed tiles and verified uploads
are reused. Network failures remain pending and are retried with delayed waves.
`--max-fetch-tiles N --passes 1` makes a bounded pilot; this may publish a partial
snapshot, which explicitly identifies unresolved coverage. The log prints each
published snapshot ID and remaining tile count. A batch with many API failures
does not prevent attempting other batches in the same wave.

To extend or correct a version, use the same output store and add
`--parent /absolute/path/to/snapshots/<SHA256>.json`. Without a refresh request,
only missing tiles are acquired. Include every source run from the parent in
the run-ID list. Add `--refresh-region 1:START-END` to reacquire
particular tiles, bypassing old cached responses for those tiles only. An
interrupted refresh resumes its own cache. Keep snapshot IDs used in papers;
there is no implicit mutable `latest` selection.

## Consume in analytics

Download an explicitly selected snapshot and its referenced fragments once:

```bash
python -m reannotation download \
  --remote-root gdrive:GAPH/reannotation/gnomad --snapshot "$SNAPSHOT_SHA256" \
  --destination "$GAPH_ROOT/reannotation/gnomad"
```

Pass `--gnomad-reannotation /absolute/path/to/snapshots/<SHA256>.json` to
`analytics.strategy_report` (after `--` in its Slurm launcher). Omitting it
retains the original analysis. It works with local runs and the existing
sequential Drive source loader. All selected source runs must be bound by the
snapshot. No network lookups occur in analytics.

Analytics verifies the snapshot, normalizes provider variants with each source's
target FASTA, and caches a narrow lookup. It joins only original failed regions
with missing gnomAD AF and valid variant keys. Existing successful annotations
remain unchanged. Positive matches can resolve a partially repaired region;
absence is accepted only after all padded tiles for that region succeeded.
Unresolved rows retain failed status. Summary, filtering, and ortholog-support
calculations use the same explicit status. Cache identities and report provenance
include the snapshot ID. Prepared Drive packs materialize the result once, so
subsequent report reads do not repeatedly join the original large tables.

## Checks

`tests/test_reannotation.py` exercises archive-only planning, source immutability,
empty versus unresolved responses, three versions, refresh and resume, transfer
verification, binding rejection, corruption detection, and request pacing.
The existing annotation, aggregation, filtering, support-cache, report, and
launcher tests cover the consumers. Changes to public contracts belong here;
shared transport semantics belong in the Stage 3 contract.
