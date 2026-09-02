"""Shared command-line contract for the local Ensembl VEP cache."""

from __future__ import annotations

from pathlib import Path


VEP_ASSEMBLY = "GRCh38"
VEP_SPECIES = "homo_sapiens"


def local_vep_cache_args(*, release: str, cache_dir: str | Path) -> list[str]:
    """Return the cache-selection flags shared by probes and annotations."""

    normalized_release = str(release).strip()
    if not normalized_release:
        raise ValueError("Local VEP requires an explicit release")
    return [
        "--offline",
        "--cache",
        "--refseq",
        "--use_given_ref",
        "--species",
        VEP_SPECIES,
        "--assembly",
        VEP_ASSEMBLY,
        "--cache_version",
        normalized_release,
        "--dir_cache",
        str(cache_dir),
    ]
