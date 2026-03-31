# minimal placeholder so pipeline runs

with open(snakemake.output[0], "w") as f:
    f.write("QV\n")
    f.write("42\n")