sinteractive --time=04:00:00 --nodes=1 --ntasks=4 --mem=32000M
module load openmpi/4.0.3/gnu/9.2.0

g++ -I/share/apps/common/openmpi/4.0.3/gnu/9.2.0/include -L/share/apps/common/openmpi/4.0.3/gnu/9.2.0/lib -lmpi kmeans_MPI_regblock.cpp  -O3 -o op_kmeans_MPI
mpirun -np 4 op_kmeans_MPI <n> <m> <k> <itr_max>