ml purge
ml snakemake/6.1.0-foss-2020b

snakefile="/data/gent/vo/000/gvo00027/vsc42927/scripts/DMRfinder.snakefile"

snakemake -s ${snakefile} --cluster "qsub -V" --default-resources --jobs 100 --rerun-incomplete
