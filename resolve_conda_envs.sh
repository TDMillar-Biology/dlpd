# Compute nodes on HPRC do not have internet connection to resolve conda envs
# So we resolve them on the log in node before running the pipeline
snakemake --use-conda --conda-create-envs-only
