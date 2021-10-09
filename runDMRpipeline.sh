ml purge
ml snakemake/6.1.0-foss-2020b

export LD_LIBRARY_PATH="/apps/gent/CO7/skylake-ib/software/Python/3.8.6-GCCcore-10.2.0/lib"

snakefile="/data/gent/vo/000/gvo00027/vsc42927/scripts/DMRpipeline/DMRfinder.snakefile"

snakemake -s ${snakefile} --cluster "qsub" --jobs 100 --rerun-incomplete
