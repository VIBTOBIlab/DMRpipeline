ml purge
ml snakemake/6.1.0-foss-2020b

snakefile="/data/gent/vo/000/gvo00027/vsc42927/scripts/DMRpipeline/DMRfinder.snakefile"
profile="/data/gent/vo/000/gvo00027/vsc42927/scripts/DMRpipeline/slurm_profile"

snakemake -s ${snakefile} --cluster "sbatch" --jobs 100
