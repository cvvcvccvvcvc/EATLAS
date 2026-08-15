from __future__ import annotations

import re
import sys
from pathlib import Path
from types import ModuleType

import pytest


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "bin"))

import ensembl_compara_maf  # noqa: E402
import run_bwa_pseudoreads  # noqa: E402
import run_minimap2_alignment  # noqa: E402
import run_nucmer_alignment  # noqa: E402
from alignment_table_schema import (  # noqa: E402
    ALIGNER_OUTPUT_SCHEMAS,
    EVENT_FIELDS,
    FAILURE_FIELDS,
    SEGMENT_FIELDS,
    SUMMARY_FIELDS,
)


STRATEGY_PRODUCERS = {
    "minimap2_asm10": run_minimap2_alignment,
    "minimap2_asm20": run_minimap2_alignment,
    "nucmer": run_nucmer_alignment,
    "bwa_pseudoreads_150_75": run_bwa_pseudoreads,
    "precomputed_ensembl_92_mammals_epo_extended": ensembl_compara_maf,
}


def registered_strategy_names() -> set[str]:
    workflow = (PROJECT_DIR / "main.nf").read_text()
    registry = re.search(
        r"ALIGNMENT_STRATEGY_REGISTRY\s*=\s*\[(.*?)\]\nAVAILABLE_ALIGNMENT_STRATEGIES",
        workflow,
        re.DOTALL,
    )
    assert registry is not None
    registry_body = registry.group(1)
    literal_names = set(re.findall(r"\bname:\s*'([^']+)'", registry_body))
    ensembl_name = re.search(
        r"ENSEMBL_COMPARA_STRATEGY\s*=\s*'([^']+)'",
        workflow,
    )
    assert ensembl_name is not None
    if "name: ENSEMBL_COMPARA_STRATEGY" in registry_body:
        literal_names.add(ensembl_name.group(1))
    return literal_names


def test_contract_covers_every_registered_alignment_strategy() -> None:
    assert set(STRATEGY_PRODUCERS) == registered_strategy_names()


@pytest.mark.parametrize(
    ("strategy", "producer"),
    STRATEGY_PRODUCERS.items(),
)
def test_registered_strategy_uses_canonical_table_schemas(
    strategy: str,
    producer: ModuleType,
) -> None:
    assert producer.SUMMARY_FIELDS is SUMMARY_FIELDS, strategy
    assert producer.SEGMENT_FIELDS is SEGMENT_FIELDS, strategy
    assert producer.EVENT_FIELDS is EVENT_FIELDS, strategy
    assert producer.FAILURE_FIELDS is FAILURE_FIELDS, strategy


def test_aligner_output_schema_registry_is_exact_and_unambiguous() -> None:
    assert ALIGNER_OUTPUT_SCHEMAS == {
        "ortholog_alignment_summary.tsv.gz": SUMMARY_FIELDS,
        "alignment_segments.tsv.gz": SEGMENT_FIELDS,
        "alignment_events.tsv.gz": EVENT_FIELDS,
        "failures.tsv.gz": FAILURE_FIELDS,
    }
    for fields in ALIGNER_OUTPUT_SCHEMAS.values():
        assert len(fields) == len(set(fields))
