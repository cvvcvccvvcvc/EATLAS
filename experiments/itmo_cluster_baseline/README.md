# ITMO Cluster Baseline Measurements

This file preserves dated measurements that informed resource and infrastructure
decisions. It is historical evidence, not a launch runbook or a source of
current defaults. Reproduce a measurement on the current commit, workload,
cache state, and cluster allocation before using it for tuning.

Current setup and diagnosis are documented in
`docs/itmo_cluster.md`. Current process resources and concurrency live in
`nextflow.config` and `nextflow_schema.json`.

## Infrastructure Snapshot

Verified through 2026-07-16:

- external login host: `ctlab.itmo.ru`;
- controller: `sphinx`;
- home: `/nfs/home/$USER`, with a reported 20 GB quota;
- work allocation: `/mnt/tank/scratch/$USER`, with a reported 250 GB quota;
- `/mnt/tank` was shared between controller and compute nodes;
- test allocations ran on `meduza-1`;
- compute-node outbound DNS/TLS connectivity was present;
- system Java on `sphinx` was OpenJDK 11;
- `normal` QOS reported up to 64 CPU, 512 GB RAM, and 2 GPU per user;
- `high_mem` was not authorized for the account;
- the user association did not report explicit `MaxSubmitJobs` or `MaxJobs`;
- five independent 1 CPU / 1 GB jobs ran concurrently.

The last point disproved an earlier assumption that four jobs was a scheduler
limit.

The `main` partition reported `DefMemPerNode=UNLIMITED`. A task without an
explicit memory request inherited a very large request and was held by
`QOSMaxMemoryPerUser`. The cause was a composed-profile selector replacing
base process resources. The production configuration moved scratch policy to a
`task_scratch` label so the Slurm profile could disable scratch without
replacing CPU, memory, time, retry, or environment directives.

An NCBI `curl -I` probe returned HTTP 403 while still proving routing, DNS, and
TLS. It was not sufficient to validate application access; the one-gene smoke
remained the integration check.

## Alignment CPU Benchmark

A two-gene benchmark used AFDN (GeneID 4301) and BRCA1 (GeneID 672) to compare
two, three, and four CPUs per Minimap2/Nucmer task.

- Three CPUs reduced individual alignment time by about 31–36% relative to two.
- Three CPUs was about 4–6% slower than four.
- Two and three CPUs consumed approximately the same reserved CPU-seconds.
- Four CPUs consumed about 27% more reserved CPU-seconds than three.

The two- and three-CPU outputs had identical summaries, derived feature
coverage, and canonical segment/event evidence. Minimap2 initially differed
only in line-order-derived provenance identifiers; those identifiers were later
made content-derived and stable across thread counts.

This measurement motivated three-CPU Minimap2 and Nucmer requests at the time.
BWA and other strategies were not covered.

## Alignment Memory Tail

Three production runs measured a size-dependent tail in the long-pseudoread
Minimap2 strategy:

- 1,405 of 1,470 logical tasks below 150 million selected-ortholog bases stayed
  below 6.5 GB RSS without retry;
- among 63 genes from 150 to 600 million bases, maximum successful RSS was
  21.5 GB;
- two genes at or above 600 million bases reached 31.8 GB;
- starting all tasks at 8 GB caused 38 failed attempts across 29 genes.

The resulting implementation records total selected-ortholog sequence volume
and selects initial pseudoread-task memory classes from it. Exact current
thresholds remain code/configuration, not this historical note.

## Annotation Memory And Concurrency

In the same measured run, two-gene annotation partitions ranged from 66,921 to
815,531 unique variant contexts and used about 0.38–5.6 GB RSS. Later resource
classes used compact alignment-event counts as the input-size measure.

The account's reported 512 GB per-user memory limit meant that the largest
concurrent first attempts could remain pending even when per-process
`maxForks` allowed them. This is why fork limits must not be read as guaranteed
simultaneous execution.

A 20-gene default-strategy smoke completed without task retries or scheduler
resource failures at alignment fork limit 4. Median measured CPU use was
2.32–2.53 for Minimap2 and 2.04 for Nucmer. It validated the task shape but did
not establish a throughput optimum.

## Production Run Snapshot

The 490-gene
`all_genes_batch_001_20260830_163720_f4_a21_n8` run on commit `8a6f666`
completed 2,347 tasks in 4 h 51 min 49 s without failed tasks or recorded
fetch, alignment, or annotation failures.

Its settings were 4 fetch tasks, 21 tasks per alignment process, and 8 tasks per
annotation process. They became defaults after that run, but this observation
does not make them permanently optimal. Future comparisons must use the same
workload, cache state, scientific inputs, and Slurm constraints.

## Shared gnomAD Cache Snapshot

On 2026-08-09, the cache at
`/mnt/tank/scratch/$USER/gaph_v2/cache/gnomad` contained 3,306 tile files and
occupied about 710 MB for the measured account. The cache is incremental; this
is not an expected current size.

## phyloP Benchmark

Verified on 2026-08-09 against saved target-null annotations from the 590-gene
run:

- all 172,183 allele rows matched the saved remote result, including 458
  positions without a score;
- 172,043 unique positions were read in 638 blocks;
- BigWig read time was 707.361 seconds over UCSC HTTPS and 6.554 seconds from
  the shared local file, about 108× faster;
- the complete local annotation pass took 10.184 seconds with about 269 MiB
  peak RSS.

This comparison covered the phyloP annotation pass only, not the whole report.

## VEP Verification And Benchmark

Verified on 2026-07-29:

- Ensembl image `ensemblorg/ensembl-vep:release_116.0`;
- local image size about 281 MB;
- RefSeq release-116 GRCh38 indexed cache about 25 GB;
- retained cache archive size 26,409,563,680 bytes;
- wrapper installed from `bin/gaph-vep116`.

The official archive was fully validated, extracted outside its final path, and
checked with VEP `--show_cache_info`. Keeping both archive and extracted cache
used about 50 GB.

A six-variant APC/BRCA1 smoke covering SNV, insertion, and deletion matched
Ensembl REST release 116 after normalizing VEP's `-` placeholder to empty.
Consequence terms, selected RefSeq transcript, canonical/MANE fields, impact,
and variant class agreed. Local VEP produced six results in about 3.3 seconds.

A historical 118,549,428-row candidate estimate implied 475 shards at the then
current shard bound. It was a planning estimate, not a completed pipeline
measurement.

A 500,000-row cache benchmark sampled 40 partitions. One-megabase and
five-megabase result-cache tiles both occupied about 1.88 MB and had about
8-second lookup time. The finer tile was retained because it improved sparse
regional pruning without measured cost in that sample.
