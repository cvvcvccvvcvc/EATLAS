# Google Drive Report Parity

This experiment compares the new archive-backed report for batches 001–010
with the preserved completed local-source report. It is not a production
launch path. Submit reports with the documented analytics launcher.

Run `python -m experiments.gdrive_report.compare --help` for its interface.
After a successful report, compare scientific scalar counts, all cached
variant-summary tables, nine derived tables, and every Plotly trace (including
support-filter curves). Provenance, timing, generated element IDs, and layout
are excluded. Integer counts compare exactly; floating results allow 1e-12
rounding differences from equivalent aggregation orders.

The result is written atomically, including a concrete failure on any mismatch.
An exit code of zero confirms the full comparison passed. Runtime and memory
come from the report's performance profile and Slurm accounting, not this
comparison's runtime.
