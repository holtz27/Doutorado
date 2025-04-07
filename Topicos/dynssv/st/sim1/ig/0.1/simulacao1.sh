#PBS -N ig_0.1
#PBS -l ncpus=1
#PBS -l walltime=200:00:00
#PBS -m ae
#PBS -M brunoeholtz@gmail.com
cd /home/beholtz/topicos/sim1/ig/0.1
module load gcc/4.9.2
module load R/4.3.3
R CMD BATCH simulacao1.R