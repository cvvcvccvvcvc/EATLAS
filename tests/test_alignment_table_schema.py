from __future__ import annotations

import re
from pathlib import Path
from types import ModuleType

import pytest


PROJECT_DIR = Path(__file__).resolve().parents[1]
from bin import (
    run_bwa_pseudoreads,
    run_minimap2_alignment,
    run_nucmer_alignment,
)
from bin.alignment_table_schema import (
    ALIGNER_OUTPUT_SCHEMAS,
    EVENT_FIELDS,
    FAILURE_FIELDS,
    SEGMENT_FIELDS,
    SUMMARY_FIELDS,
)


STRATEGY_PRODUCERS = {
    "minimap2_asm10": run_minimap2_alignment,
    "minimap2_asm20": run_minimap2_alignment,
    "minimap2_map_ont_pseudoreads_30000_15000": run_minimap2_alignment,
    "nucmer": run_nucmer_alignment,
    "bwa_pseudoreads_150_75": run_bwa_pseudoreads,
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
    return literal_names


def test_contract_covers_every_registered_alignment_strategy() -> None:
    assert set(STRATEGY_PRODUCERS) == registered_strategy_names()


def test_map_ont_long_pseudoread_strategy_is_fixed_and_opt_in() -> None:
    workflow = (PROJECT_DIR / "main.nf").read_text()
    strategy = re.search(
        r"\[\s*name:\s*'minimap2_map_ont_pseudoreads_30000_15000',(.*?)\n\s*\]",
        workflow,
        re.DOTALL,
    )
    assert strategy is not None
    strategy_body = strategy.group(1)
    assert "default_enabled: false" in strategy_body
    assert "minimap2_preset: 'map-ont'" in strategy_body
    assert "minimap2_pseudoread_len: 30000" in strategy_body
    assert "minimap2_pseudoread_step: 15000" in strategy_body


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
