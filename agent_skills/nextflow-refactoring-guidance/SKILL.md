---
name: nextflow-refactoring-guidance
description: Provides guidance for refactoring Nextflow pipelines towards nf-core and modern Software Engineering standards (DSL2, nf-test, dynamic resources). Use when asked to refactor, modularize, or improve a bioinformatics pipeline.
---
# Nextflow Pipeline Refactoring Guidance

This skill outlines the "golden practices" for modernizing Nextflow pipelines (like GAPH v2). Apply these steps when refactoring bioinformatics workflows to ensure modularity, resilience, and reproducibility.

## Contents
- [1. Architecture & Modularity (DSL2)](#1-architecture--modularity-dsl2)
- [2. Configuration & Parameter Management](#2-configuration--parameter-management)
- [3. Error Handling & Resilience](#3-error-handling--resilience)
- [4. Testing & QA (nf-test)](#4-testing--qa-nf-test)
- [5. Workflow Checklist](#5-workflow-checklist)

## 1. Architecture & Modularity (DSL2)
Enforce strict Separation of Concerns. Move away from monolithic `main.nf` scripts.
*   **Modules (`modules/`):** Extract individual processes (e.g., `minimap2`, `bwa`) into single-purpose `.nf` files. Each module must define strictly typed `input` and `output` channels using explicit `emit` names.
*   **Subworkflows (`subworkflows/`):** Group related modules (e.g., `align_and_sort.nf`) to allow testing of larger pipeline components.
*   **Meta Maps:** Always flow data with a metadata dictionary (e.g., `tuple val(meta), path(reads)`) to keep sample properties attached to files throughout execution.

## 2. Configuration & Parameter Management
Remove compute/environment definitions from `.nf` files.
*   **Configuration Hierarchy:** Use `nextflow.config` and nested profiles (`conf/profiles/`) for environments (e.g., `local`, `slurm`).
*   **Schema Validation:** Define pipeline parameters in a `nextflow_schema.json` to enable automatic validation and UI generation, instead of manual `if (!params.foo)` checks.
*   **Tool Configs (`conf/modules.config`):** Isolate tool-specific CLI arguments (using `ext.args`) into configuration blocks.

## 3. Error Handling & Resilience
Design processes to be fault-tolerant and cacheable.
*   **Dynamic Resources:** Implement adaptive retry logic for system limits (e.g., OOM):
    ```nextflow
    process {
        errorStrategy = { task.exitStatus in [104, 134, 137, 140, 143] ? 'retry' : 'finish' }
        maxRetries    = 3
        memory        = { 4.GB * task.attempt }
    }
    ```
*   **Domain Errors:** For biological errors (e.g., "no orthologs found"), catch the error in the script, exit gracefully (return `0`), and emit an empty log to a `failed_samples` channel instead of crashing the pipeline.

## 4. Testing & QA (nf-test)
Adopt `nf-test` as the standard testing framework.
*   **Module & Snapshot Testing:** Write `.nf-test` files for every module to verify outputs using MD5 checksum snapshots.
*   **Test Profile:** Create a `test` profile in `nextflow.config` that automatically fetches a minimal dataset to run the pipeline end-to-end in CI/CD.

## 5. Workflow Checklist
When refactoring a pipeline, follow these steps:
- [ ] **Step 1:** Isolate configuration. Move all `process` directives (cpus, conda, maxForks) from `main.nf` to `nextflow.config` blocks.
- [ ] **Step 2:** Modularize processes. Extract processes from `main.nf` into `modules/local/[process_name].nf`.
- [ ] **Step 3:** Implement `meta` maps. Refactor channels to pass `tuple val(meta), path(file)`.
- [ ] **Step 4:** Add dynamic error handling. Apply the retry-with-multiplier pattern for memory and time.
- [ ] **Step 5:** Write tests. Initialize `nf-test` and write snapshots for core modules.
