#!/bin/sh

#SBATCH --job-name="ParPrimeSieve"
#SBATCH --partition=compute
#SBATCH --account=education-eemcs-courses-mastermath
#SBATCH --time=00:10:00
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2024r1
module load openmpi
module load opencoarrays

make clean
make parallel_phantom

#np=64
#npx=512
lambda=0.1
sigma=0.5

# Create the file output.log with as first line "it, error, time, np\n"
echo "it, error, time, np, n" > output.out

# Run the program with, 4, 9, 16, 25, 36, 49, 64 processors
# with npx 512, 1024, 2

for n in 4 9 16 25 36 49 64
do
	for npx in 512 1024 2048 5096
	do
		srun -n $n ./parallel_phantom  -npx $npx -lambda $lambda -sigma $sigma
	done
done
srun -n 64 ./parallel_phantom  -npx 512 -lambda $lambda -sigma $sigma