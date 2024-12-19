#!/bin/sh

#SBATCH --job-name="ParPrimeSieve"
#SBATCH --partition=compute
#SBATCH --account=education-eemcs-courses-mastermath
#SBATCH --time=00:30:00
#SBATCH --ntasks=48
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

# Create the file output.log with as first line "npx, np, total_time, comm_time"
echo "npx, np_row, np_col, total_time, comm_time, comp_time" > output.out

#for npx in 
# for np_row in 1 2 3 4 5 6
# do
# 	for np_col in 1 2 3 4 5 6
# 	do
# 		n=$((np_row * np_col))
# 		srun -n $n ./parallel_phantom  -npx 8192 -lambda $lambda -sigma $sigma -np_row $np_row -np_col $np_col
# 	done
# done



# for np_row in 1 2 3 4 5 6
# do
# 	for np_col in 1 2 3 4 5 6
# 	do
# 		n=$((np_row * np_col))
# 		srun -n $n ./parallel_phantom  -npx 5096 -lambda $lambda -sigma $sigma -np_row $np_row -np_col $np_col
# 	done
# done


############# ABANDON ALL HOPE YE WHO ENTER HERE ################

# 512
# 725
# 887
# 1024
# 1145
# 1255
# 1355
# 1449
# 1536
# 1620
# 1699
# 1774
# 1847
# 1916
# 1983
# 2048
# 2112
# 2173
# 2232
# 2290
# 2347
# 2402
# 2456
# 2509
# 2560
# 2611
# 2661
# 2710
# 2758
# 2805
# 2851
# 2897
# 2942
# 2986
# 3030
# 3072
# 3115
# 3157
# 3198
# 3239
# 3279
# 3319
# 3358
# 3397
# 3435
# 3473
# 3511
# 3548

# srun -n 1 ./parallel_phantom -npx 512 -lambda $lambda -sigma $sigma -np_row 1 -np_col 1
# srun -n 2 ./parallel_phantom -npx 725 -lambda $lambda -sigma $sigma -np_row 1 -np_col 2
# srun -n 3 ./parallel_phantom -npx 887 -lambda $lambda -sigma $sigma -np_row 1 -np_col 3
# srun -n 4 ./parallel_phantom -npx 1024 -lambda $lambda -sigma $sigma -np_row 2 -np_col 2
# srun -n 5 ./parallel_phantom -npx 1145 -lambda $lambda -sigma $sigma -np_row 1 -np_col 5
# srun -n 6 ./parallel_phantom -npx 1255 -lambda $lambda -sigma $sigma -np_row 2 -np_col 3
# srun -n 7 ./parallel_phantom -npx 1355 -lambda $lambda -sigma $sigma -np_row 1 -np_col 7
# srun -n 8 ./parallel_phantom -npx 1449 -lambda $lambda -sigma $sigma -np_row 2 -np_col 4
# srun -n 9 ./parallel_phantom -npx 1536 -lambda $lambda -sigma $sigma -np_row 3 -np_col 3
# srun -n 10 ./parallel_phantom -npx 1620 -lambda $lambda -sigma $sigma -np_row 2 -np_col 5
# srun -n 11 ./parallel_phantom -npx 1699 -lambda $lambda -sigma $sigma -np_row 1 -np_col 11
# srun -n 12 ./parallel_phantom -npx 1774 -lambda $lambda -sigma $sigma -np_row 3 -np_col 4
# srun -n 13 ./parallel_phantom -npx 1847 -lambda $lambda -sigma $sigma -np_row 1 -np_col 13
# srun -n 14 ./parallel_phantom -npx 1916 -lambda $lambda -sigma $sigma -np_row 2 -np_col 7
# srun -n 15 ./parallel_phantom -npx 1983 -lambda $lambda -sigma $sigma -np_row 3 -np_col 5
# srun -n 16 ./parallel_phantom -npx 2048 -lambda $lambda -sigma $sigma -np_row 4 -np_col 4
# srun -n 17 ./parallel_phantom -npx 2112 -lambda $lambda -sigma $sigma -np_row 1 -np_col 17
# srun -n 18 ./parallel_phantom -npx 2173 -lambda $lambda -sigma $sigma -np_row 3 -np_col 6
# srun -n 19 ./parallel_phantom -npx 2232 -lambda $lambda -sigma $sigma -np_row 1 -np_col 19
# srun -n 20 ./parallel_phantom -npx 2290 -lambda $lambda -sigma $sigma -np_row 4 -np_col 5
# srun -n 21 ./parallel_phantom -npx 2347 -lambda $lambda -sigma $sigma -np_row 3 -np_col 7
# srun -n 22 ./parallel_phantom -npx 2402 -lambda $lambda -sigma $sigma -np_row 2 -np_col 11
# srun -n 23 ./parallel_phantom -npx 2456 -lambda $lambda -sigma $sigma -np_row 1 -np_col 23
# srun -n 24 ./parallel_phantom -npx 2509 -lambda $lambda -sigma $sigma -np_row 6 -np_col 4
# srun -n 25 ./parallel_phantom -npx 2560 -lambda $lambda -sigma $sigma -np_row 5 -np_col 5
# srun -n 27 ./parallel_phantom -npx 2611 -lambda $lambda -sigma $sigma -np_row 3 -np_col 9
# srun -n 28 ./parallel_phantom -npx 2661 -lambda $lambda -sigma $sigma -np_row 4 -np_col 7
# srun -n 29 ./parallel_phantom -npx 2710 -lambda $lambda -sigma $sigma -np_row 1 -np_col 29
# srun -n 30 ./parallel_phantom -npx 2758 -lambda $lambda -sigma $sigma -np_row 5 -np_col 6
# srun -n 31 ./parallel_phantom -npx 2805 -lambda $lambda -sigma $sigma -np_row 1 -np_col 31
# srun -n 32 ./parallel_phantom -npx 2851 -lambda $lambda -sigma $sigma -np_row 4 -np_col 8
# srun -n 33 ./parallel_phantom -npx 2897 -lambda $lambda -sigma $sigma -np_row 3 -np_col 11
# srun -n 34 ./parallel_phantom -npx 2942 -lambda $lambda -sigma $sigma -np_row 2 -np_col 17
# srun -n 35 ./parallel_phantom -npx 2986 -lambda $lambda -sigma $sigma -np_row 5 -np_col 7
# srun -n 36 ./parallel_phantom -npx 3030 -lambda $lambda -sigma $sigma -np_row 6 -np_col 6
# srun -n 37 ./parallel_phantom -npx 3072 -lambda $lambda -sigma $sigma -np_row 1 -np_col 37
# srun -n 38 ./parallel_phantom -npx 3115 -lambda $lambda -sigma $sigma -np_row 2 -np_col 19
# srun -n 39 ./parallel_phantom -npx 3157 -lambda $lambda -sigma $sigma -np_row 3 -np_col 13
# srun -n 40 ./parallel_phantom -npx 3198 -lambda $lambda -sigma $sigma -np_row 5 -np_col 8
# srun -n 41 ./parallel_phantom -npx 3239 -lambda $lambda -sigma $sigma -np_row 1 -np_col 41
# srun -n 42 ./parallel_phantom -npx 3279 -lambda $lambda -sigma $sigma -np_row 6 -np_col 7
# srun -n 43 ./parallel_phantom -npx 3319 -lambda $lambda -sigma $sigma -np_row 1 -np_col 43
# srun -n 44 ./parallel_phantom -npx 3358 -lambda $lambda -sigma $sigma -np_row 4 -np_col 11
# srun -n 45 ./parallel_phantom -npx 3397 -lambda $lambda -sigma $sigma -np_row 5 -np_col 9
# srun -n 46 ./parallel_phantom -npx 3435 -lambda $lambda -sigma $sigma -np_row 2 -np_col 23
# srun -n 47 ./parallel_phantom -npx 3473 -lambda $lambda -sigma $sigma -np_row 1 -np_col 47
# srun -n 48 ./parallel_phantom -npx 3511 -lambda $lambda -sigma $sigma -np_row 6 -np_col 8

