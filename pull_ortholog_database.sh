# TAMU hprc grace env
#module load GCC/13.2.0 OpenMPI/4.1.6 Compleasm/0.2.7
#compleasm.py download -L /scratch/group/chakraborty_lab/trevor/ortholog_databases  --odb odb10 diptera
#compleasm.py download -L /scratch/group/chakraborty_lab/trevor/ortholog_databases  --odb odb12 diptera

module load GCC/12.2.0  OpenMPI/4.1.4 BUSCO/5.7.1

mkdir -p /scratch/group/chakraborty_lab/trevor/ortholog_databases


#busco --download diptera_odb10 \
#  --download_path /scratch/group/chakraborty_lab/trevor/ortholog_databases \
#  -v

#wget "https://busco-data.ezlab.org/v5/data/lineages/diptera_odb10.2024-01-08.tar.gz"
#wget "https://busco-data.ezlab.org/v5/data/lineages/diptera_odb12.2025-07-01.tar.gz"
#tar -xzf diptera_odb10*.tar.gz
#tar -xzf diptera_odb12*.tar.gz


ln -s /scratch/group/chakraborty_lab/trevor/ortholog_databases/diptera_odb10 \
  data/ortholog_databases/diptera_odb10

ln -s /scratch/group/chakraborty_lab/trevor/ortholog_databases/diptera_odb12 \
  data/ortholog_databases/diptera_odb12

