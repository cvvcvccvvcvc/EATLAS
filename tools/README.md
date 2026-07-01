Local tool links
================

`tools/bin/` is intentionally ignored. Use it for local symlinks or copies of
external command-line tools that should not be committed to the repository.

The workflow resolves the NCBI Datasets CLI in this order:

1. `--datasets_bin` / `DATASETS_BIN`
2. `tools/bin/datasets`
3. `datasets` on `PATH`
