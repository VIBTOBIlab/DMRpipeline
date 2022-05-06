ml purge
ml snakemake/6.1.0-foss-2020b

snakefile="/data/gent/vo/000/gvo00027/vsc42927/scripts/DMRpipeline/DMRfinder.snakefile"

snakemake -s ${snakefile} --cluster "sbatch --time=4:00:00 --cpus-per-task=11" --jobs 100
