import urllib.request
import json

GNOMAD_GENE_QUERY = """
query VariantsInGene($gene_symbol: String!) {
  gene(gene_symbol: $gene_symbol, reference_genome: GRCh38) {
    gene_id
    symbol
    variants(dataset: gnomad_r4) {
      variant_id
      chrom
      pos
      ref
      alt
      joint { ac an }
    }
  }
}
"""

req = urllib.request.Request(
    "https://gnomad.broadinstitute.org/api",
    method="POST",
    headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"},
    data=json.dumps({"query": GNOMAD_GENE_QUERY, "variables": {"gene_symbol": "BRCA1"}}).encode()
)

with urllib.request.urlopen(req) as response:
    data = json.loads(response.read())
    variants = data.get("data", {}).get("gene", {}).get("variants", [])
    print(f"Found {len(variants)} variants for BRCA1")
