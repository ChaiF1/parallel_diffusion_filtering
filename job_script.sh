#!/bin/sh

#SBATCH --job-name="ParPrimeSieve"
#SBATCH --partition=compute
#SBATCH --account=education-eemcs-courses-mastermath
#SBATCH --time=00:20:00
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load 2024r1
module load openmpi
module load opencoarrays

make clean
make parallel_phantom

#np=64
npx=512
lambda=0.1
sigma=0.5

# Create the file output.log with as first line "npx, np, total_time, comm_time"
echo "npx, np_row, np_col, total_time, comm_time, comp_time" > output.out


while :
do
	srun -n 48 ./parallel_phantom -npx $npx -lambda $lambda -sigma $sigma -np_row 6 -np_col 8
	npx=$((npx*2))
done

# for np_row in 1 2 3 4 5 6
# do
# 	for np_col in 1 2 3 4 5 6
# 	do
# 		n=$((np_row * np_col))
# 		srun -n $n ./parallel_phantom  -npx 5096 -lambda $lambda -sigma $sigma -np_row $np_row -np_col $np_col
# 	done
# done

# for np_row in 6
# do
# 	for np_col in 3 4 5 6
# 	do
# 		n=$((np_row * np_col))
# 		srun -n $n ./parallel_phantom  -npx 10192 -lambda $lambda -sigma $sigma -np_row $np_row -np_col $np_col
# 	done
# done