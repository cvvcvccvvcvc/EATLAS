"""Reconcile allele-level provider evidence across target-gene contexts."""

from __future__ import annotations

from genomics.variants import ALLELE_ANNOTATION_FIELDS


def allele_evidence_comparison_sql(field: str, *, alias: str = "") -> str:
    """Return the canonical SQL value used to compare allele evidence."""

    if field not in ALLELE_ANNOTATION_FIELDS:
        raise ValueError(f"Unknown allele annotation field: {field}")
    column = f"{alias}.{field}" if alias else field
    value = f"nullif({column}, '')"
    if field == "gnomad_af":
        return f"try_cast({value} AS DOUBLE)"
    return value


def materialize_allele_evidence(
    connection,
    *,
    relation: str,
    key: str,
    fields: tuple[str, ...] = ALLELE_ANNOTATION_FIELDS,
) -> None:
    if "gnomad_af" not in fields:
        raise ValueError("Allele evidence materialization requires gnomad_af")
    resolved = ", ".join(
        f"max(nullif({field}, '')) AS {field}" for field in fields
    )
    counts = ", ".join(
        f"count(DISTINCT {allele_evidence_comparison_sql(field)}) FILTER "
        f"(WHERE {allele_evidence_comparison_sql(field)} IS NOT NULL) "
        f"AS {field}_count"
        for field in fields
    )
    invalid_af = (
        "count(*) FILTER (WHERE nullif(gnomad_af, '') IS NOT NULL AND "
        "try_cast(nullif(gnomad_af, '') AS DOUBLE) IS NULL) "
        "AS gnomad_af_invalid_count"
    )
    connection.execute(
        f"CREATE TEMP TABLE allele_evidence AS SELECT {key}, "
        f"{resolved}, max(try_cast(nullif(gnomad_af, '') AS DOUBLE)) "
        f"AS gnomad_af_value, {counts}, {invalid_af} "
        f"FROM {relation} GROUP BY {key}"
    )

    count_columns = [f"{field}_count" for field in fields]
    conflict_filter = " OR ".join(
        [*(f"{column} > 1" for column in count_columns), "gnomad_af_invalid_count > 0"]
    )
    rows = connection.execute(
        f"SELECT {key}, "
        + ", ".join([*count_columns, "gnomad_af_invalid_count"])
        + f" FROM allele_evidence WHERE {conflict_filter} LIMIT 10"
    ).fetchall()
    conflicts = [
        (field, str(row[0]))
        for row in rows
        for field, count in zip(fields, row[1:-1])
        if int(count) > 1
    ]
    conflicts.extend(
        ("gnomad_af invalid", str(row[0]))
        for row in rows
        if int(row[-1]) > 0
    )
    if conflicts:
        examples = ", ".join(
            f"{variant} ({field})" for field, variant in conflicts[:10]
        )
        raise ValueError(
            "Variant annotations contain conflicting allele-level external "
            f"evidence: {examples}"
        )
