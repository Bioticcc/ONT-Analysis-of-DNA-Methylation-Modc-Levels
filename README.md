# ONT 5mC/5hmC analysis

This project contains the reproducible code and pinned command-line tools for
the `Pa2_gDNA-2` Oxford Nanopore methylation analysis. Large inputs,
references, intermediate files, and results are kept on the 20 TB data drive.

## Storage layout

Primary drive (this directory):

```text
ONTAnalysis/
├── bin/                 Stable links to the staged executables
├── config/              Project paths and analysis settings
├── docs/                Workflow and data-layout documentation
├── scripts/             Reproducible setup and analysis scripts
└── tools/               Pinned, versioned tool installations and archives
```

20 TB drive:

```text
/mnt/20TB_A/Adam/ONTAnalysis/
├── Wei_P2a_modSUP/      Original instrument output (preserved unchanged)
├── inputs/
│   ├── reads/           Links to the original run data
│   └── references/      Analysis reference FASTA and indexes
└── outputs/
    └── Pa2_gDNA-2/
        ├── 00_manifests/
        ├── 01_alignment/        Restartable, unsorted BAM chunks
        ├── 02_alignment_qc/
        ├── 03_ont_qc_coverage/
        ├── 04_methylation/
        ├── 05_methylation_exploration/
        ├── logs/
        └── tmp/
```

The primary analysis input is the 74 unaligned pass BAM files. These retain
the Dorado `MM`, `ML`, and `MN` modified-base tags. FASTQ files are not used
for methylation calling because they do not retain those tags.

The canonical scope, directory layout, stage contracts, output requirements,
and development rules are defined in
[`docs/PROJECT_STRUCTURE.md`](docs/PROJECT_STRUCTURE.md). Contributors should
read that guide before changing the pipeline.

## Tool versions

- Dorado 2.0.0
- Minimap2 2.31
- Modkit 0.6.4
- Mosdepth 0.3.14 (supporting coverage utility, pinned by SHA-256)
- Samtools 1.22.1 (supporting utility already installed locally)

Run `scripts/stage_tools.sh` to install or verify the pinned tools, then source
`config/project.env` from analysis scripts.

## Alignment execution and logs

Run the stages in order:

```bash
./scripts/01_align_modbam.sh
./scripts/02_finalize_alignment.sh
./scripts/03_ont_qc_and_coverage.sh
./scripts/04_methylation_analysis.sh
./scripts/05-1_methylation_exploration.sh
```

Use `./scripts/01_align_modbam.sh --check-only` to validate paths, tools, and
the expected 74 source BAMs without starting an alignment.

Stage 01 aligns one source BAM at a time with Dorado `--no-sort`. Each finished
chunk is validated and saved atomically, so an interruption only requires the
active chunk to be repeated. Stage 02 refuses to proceed unless all 74 expected
chunks exist, validates their BAM structure and modification tags, then performs
one controlled coordinate sort.

Timestamped logs, including five-minute CPU/memory/disk heartbeats, are written to the 20 TB drive under
`outputs/Pa2_gDNA-2/logs`. A `*.latest.log` link points to the newest attempt.

Stage 03 produces alignment statistics, genome/chromosome coverage, 100 kb
window breadth at 1x/5x/10x/20x, exact zero/low-depth intervals, and an indexed
run-length depth track. Use `./scripts/03_ont_qc_and_coverage.sh --check-only`
to validate its prerequisites without launching the coverage pass.

Stage 04 validates the joint `C+m`/`C+h` tags and produces an indexed,
strand-combined CpG bedMethyl with separate 5mC and 5hmC records. Its global,
chromosome, and 100 kb tables use weighted call counts rather than averaging
site percentages. Use `./scripts/04_methylation_analysis.sh --check-only`
before launching the full Modkit analysis.

Stage 05 turns the validated Stage 03/04 tables into publication-ready PDF and
PNG figures plus the exact TSV data behind each plot. It includes the
single-sample-compatible global, coverage, genomic-feature, scatter, regional,
heatmap, and Manhattan-style figure families used in SixBase. The SixBase mm10
GENCODE, CpG-island, cCRE, and intergenic annotations are read in place and
recorded in the manifest. Use
`./scripts/05-1_methylation_exploration.sh --check-only` to validate all inputs,
indexes, annotations, R packages, and settings without generating figures.

Correlation, PCA, DMR/DhMR significance, differential heatmaps, and enrichment
are intentionally not generated for this one-sample dataset; Stage 05 records
those scientific boundaries in its report and figure-mapping table.
