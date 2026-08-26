# Workflow

The current dataset was already basecalled with a SUP model and the joint
all-context 5mC/5hmC modified-base model. The remaining stages are:

1. Validate the 74 pass modBAM chunks and generate a manifest.
2. Align those BAMs to standard mouse GRCm38p6 while preserving `MM`, `ML`,
   and `MN` tags.
3. Coordinate-sort and index the aligned modBAM.
4. Run mapping, coverage, and modification-tag QC.
5. Run Modkit summary and CpG/all-context pileups with separate 5mC and 5hmC
   records.
6. Calculate coverage-weighted global levels and classify missing CpG sites.

POD5 and FASTQ files are retained as source material but are not primary input
to this run of the workflow.
