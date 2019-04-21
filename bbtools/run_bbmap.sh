module load bbtools
bbmap.sh ref=../reference_sequences/GCF_000001635.26_GRCm38.p6_genomic.fna

for f in *.fastq.gz; do
    STEM=$(basename "${f}" .fastq.gz)
    bbmap.sh maxindel=160000 slow=t in="${f}" out="${STEM}".sam
done

cd $BSCRATCH/mus_2019_01_29/bbmap_run

cp ../data/sample_prefixes.txt ./


while IFS= read -r var
do
    echo "$var"
    mkdir -p results/"$var"
    echo bbmap.sh maxindel=160000 slow=t in=../data/"$var"/"$var"_L3_L4_R1_001.fastq.gz in2=../data/"$var"/"$var"_L3_L4_R2_001.fastq.gz out=results/"$var"/"$var"_L3_L4_R1_001.sam > run_"$var"_bbmap.sh
done < "sample_prefixes.txt"

#SBATCH --tasks-per-node=16


nano run_srun_run_JNMV1A_S18_bbmap.sh
#!/bin/bash -l
#SBATCH --time=11:00:00
#SBATCH --nodes=1
#SBATCH --constraint=haswell
srun ./run_JNMV1A_S18_bbmap.sh

sbatch ./run_srun_run_JNMV1A_S18_bbmap.sh

        submit_txt = os.popen('sbatch -D .%s -J %s ./slurm_submission_script.sh'%((' -N %s --mem=%d --time=%s:00:00'%(nodes,int(memory)*1000,time) if no_options == 0 else ''), job_name)).read()+'\n' # --ntasks=1 -n 32
        #subprocess.call('srun --pty /bin/bash && nohup srun -J polyCRACKER --time=%s:00:00 ./runCluster.sh &',shell=True)#'module load uge && sbatch --nodes=%s --time=%s:00:00 -J polyCRACKER -C haswell ./runCluster.sh'%(nodes,time), shell=True)# --qos=jgi -p regular  -L SCRATCH

https://ubccr.freshdesk.com/support/solutions/articles/5000688140-submitting-a-slurm-job-script


regan@cori06:~> sbatch --partition=jgi -A fungalp $(which sbatch_cori.sh) $(which test_hipmer.sh)
Submitted batch job 5956979

regan@cori06:~> sbatch --partition=debug -A fungalp $(which sbatch_cori.sh) $(which test_hipmer.sh)
Submitted batch job 5956980

The former job will schedule rapidly (if the jgi partition is not used) but fail to startup upcrun, the latter will finish with success.

#SBATCH --qos=debug

sbatch first-job.sh

#!/bin/bash
#SBATCH --time=11:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --constraint=haswell

srun check-mpi.intel.cori


sbatch -D .%s -J %s ./slurm_submission_script.sh'%((' -N %s --mem=%d --time=%s:00:00'%(nodes,int(memory)*1000,time) if no_options == 0 else ''), job_name)).read()+'\n' # --ntasks=1 -n 32

qsub -P plant-analysis.p -N run_JNMV1A_S18_bbmap -cwd -b yes -now no -j yes -m abes -M sgordon@lbl.gov -w e -l exclusive.c -l ram.c=100G -l h_rt=11:00:00 ./run_JNMV1A_S18_bbmap.sh


#!/bin/bash#
#SBATCH -p general # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 100 # memory pool for all cores
#SBATCH -t 0-2:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
for i in {1..100000}; do
echo $RANDOM >> SomeRandomNumbers.txt
donesort SomeRandomNumbers.txt


List all current jobs for a user:

squeue -u sgordon






