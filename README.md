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
        ├── 01_alignment/
        ├── 02_alignment_qc/
        ├── 03_modkit/
        ├── 04_methylation_summary/
        ├── 05_gap_analysis/
        ├── logs/
        └── tmp/
```

The primary analysis input is the 74 unaligned pass BAM files. These retain
the Dorado `MM`, `ML`, and `MN` modified-base tags. FASTQ files are not used
for methylation calling because they do not retain those tags.

## Tool versions

- Dorado 2.0.0
- Minimap2 2.31
- Modkit 0.6.4
- Samtools 1.22.1 (supporting utility already installed locally)

Run `scripts/stage_tools.sh` to install or verify the pinned tools, then source
`config/project.env` from analysis scripts.
