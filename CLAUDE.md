# cfDNA analysis project

/workspace for you = /mnt/scratch/DM/cfdna on host.

Seminal Snyder 2016 study data re-analysis.

All from the study BAMs are downloaded into /workspace/data/bam

Custom cfdna-finale-snakemake pipeline - wraps FInaleToolkit. Outputs: framgent lengths and 4-mer end motifs across fragment length classes specified in cfdna-finale-snakemake/config.yaml.

Since the config paths are relative to host, you don't run snakemake yourself - it is launched on host.

End-motif-analysis.md defines the research questions.

/workspace/analysis contains the script that answers them.

Refer to README.md for a detailed repository breakdown.