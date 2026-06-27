import pysam

clinvar = pysam.VariantFile("/Users/cvv/code/Bioinformatics/course-work/course-work-code/data/clinvar.vcf.gz")
for rec in clinvar.fetch("1", 1000000, 2000000):
    print(rec.chrom, rec.pos, rec.ref, rec.alts, rec.info.get("CLNSIG"))
    break
