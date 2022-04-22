#include <mpi.h>
//#include <math.h>
//#include <stdio.h>
#include <iostream> 
#include <string> 
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <vector>

using namespace std;

#define MIN(a,b) ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) ( BLOCK_LOW((id)+1,p,n)-1 )
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW( (id)+1, p, n) - BLOCK_LOW((id), p, n))
#define BLOCK_OWNER(index,p,n) ( ( ((p)*(index)+1)-1 ) / (n) )
#define MPI_CALL(call)                                                                \
    {                                                                                 \
        int mpi_status = call;                                                        \
        if (0 != mpi_status) {                                                        \
            char mpi_error_string[MPI_MAX_ERROR_STRING];                              \
            int mpi_error_string_length = 0;                                          \
            MPI_Error_string(mpi_status, mpi_error_string, &mpi_error_string_length); \
            if (NULL != mpi_error_string)                                             \
                fprintf(stderr,                                                       \
                        "ERROR: MPI call \"%s\" in line %d of file %s failed "        \
                        "with %s "                                                    \
                        "(%d).\n",                                                    \
                        #call, __LINE__, __FILE__, mpi_error_string, mpi_status);     \
            else                                                                      \
                fprintf(stderr,                                                       \
                        "ERROR: MPI call \"%s\" in line %d of file %s failed "        \
                        "with %d.\n",                                                 \
                        #call, __LINE__, __FILE__, mpi_status);                       \
        }                                                                             \
    }

int RandMatrixGen(double* M, long long n, int m) {

	for (long long i = 0; i < n * m; i++)
		//M[i] = 2.0 * rand() / RAND_MAX - 1.0;
		M[i] = rand() % 100;
	return 0;
}

void initializeCn(double* A, long long n, int m, double* Cn, int k) {
	for (int i = 0; i < k; i++) {
		long long e = rand() % n;
		for (int j = 0; j < m; j++) {
			Cn[i * m + j] = A[e * m + j];
		}
	}

}

void calculateDistance(double* A, double* Cn, double* D, int n, int m, int k) {

	for (int l = 0; l < k; l += 2) {
		for (int i = 0; i < n; i += 2) {
			register double dis11 = D[i * k + l];
			register double dis12 = D[i * k + l + 1];
			register double dis21 = D[(i + 1) * k + l];
			register double dis22 = D[(i + 1) * k + l + 1];
			register int r_A = i * m;
			register int r_Cn = l * m;
			for (int j = 0; j < m; j++) {
				register double A11 = A[r_A + j];
				register double A21 = A[r_A + m + j];
				register double Cn11 = Cn[r_Cn + j];
				register double Cn21 = Cn[r_Cn + m + j];
				dis11 += (A11 - Cn11) * (A11 - Cn11);
				dis12 += (A11 - Cn21) * (A11 - Cn21);
				dis21 += (A21 - Cn11) * (A21 - Cn11);
				dis22 += (A21 - Cn21) * (A21 - Cn21);
			}
			D[i * k + l] = sqrt(dis11);
			D[(i + 1) * k + l] = sqrt(dis21);
			D[i * k + l + 1] = sqrt(dis12);
			D[(i + 1) * k + l + 1] = sqrt(dis22);
		}
	}
}


void decideClusterID(double* D, int n, int k, int* Cl) {
	for (int i = 0; i < n; i++) {
		int min = 9999999;
		int clusterID = -1;
		for (int j = 0; j < k; j++) {
			if (D[i * k + j] < min) {
				min = D[i * k + j];
				clusterID = j;
			}
		}
		Cl[i] = clusterID;
	}
}


void computeCentroidSum(double* A, int n, int m, int* Cl, double* C_sum, int* Cl_n) {
	for (int i = 0; i < n; i++) {
		int clusterID = Cl[i];
		Cl_n[clusterID]++;
		for (int j = 0; j < m; j++) {
			C_sum[clusterID * m + j] += A[i * m + j];
		}

	}
}

void computeNewCentroid(int m, int k, double* Cn, int* Cl_n) {

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < m; j++) {
			if (Cl_n[i] != 0) {
				Cn[i * m + j] = Cn[i * m + j] / Cl_n[i];
			}
		}
	}
}


void printMatrix(double* matrix, int n, int m) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			cout << matrix[i * m + j] << "  ";
		}
		cout << endl;
	}
}

void printClusterID(int* Cl, int n) {
	for (int i = 0; i < n; i++) {
		cout << "point: " << i << "clusterID: " << Cl[i] << endl;
	}
}

void computeNewCentroid_validate(double* A, int n, int m, int* Cl, int k, double* Cn, int* Cl_n) {
	for (int i = 0; i < n; i++) {
		int clusterID = Cl[i];
		Cl_n[clusterID]++;
		for (int j = 0; j < m; j++) {
			Cn[clusterID * m + j] += A[i * m + j];
		}

	}
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < m; j++) {
			if (Cl_n[i] != 0) {
				Cn[i * m + j] = Cn[i * m + j] / Cl_n[i];
			}
		}
	}
}
void validate(double* A, double* Cn, int n, int m, int k, int itr) {
	int* Cl = (int*)calloc(n, sizeof(double)); //vector to store cluster id of each element
	int* Cl_n = (int*)calloc(k, sizeof(double)); //vector to store number of elements in cluster
	double* D = (double*)calloc(n * k, sizeof(double)); //matrix to store distances from centroids
	for (int i = 0; i < itr; i++) {
		calculateDistance(A, Cn, D, n, m, k);

		/*cout << "itr:" << i << " printing D: \n";
		printMatrix(D, n, k);*/


		decideClusterID(D, n, k, Cl);

		/*cout << "printing Cluster IDs: \n";
		printClusterID(Cl, n);*/

		Cn = (double*)calloc(k * m, sizeof(double)); //reinitialize all values of Cn to 0
		Cl_n = (int*)calloc(k, sizeof(double));
		computeNewCentroid_validate(A, n, m, Cl, k, Cn, Cl_n);


	}
	cout << "printing Cn from validation function: \n";
	printMatrix(Cn, k, m);
}



int main(int argc, char* argv[])
{
	int id, p;
	double elapsed_time;
	long long index, prime, first, i, count, global_count;
	MPI_CALL(MPI_Init(&argc, &argv));
	MPI_CALL(MPI_Barrier(MPI_COMM_WORLD));
	MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &id));
	MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &p));

	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = -MPI_Wtime();
	//Initialization Phase
	long long n = stoll(argv[1]); //number of elements
	int m = atoi(argv[2]); //number of features
	int k = atoi(argv[3]); //number of clusters
	int itr = atoi(argv[4]); //number of iterations

	long long low_value = BLOCK_LOW(id, p, n); //lowest Global id of the elements in the partition
	long long high_value = BLOCK_HIGH(id, p, n); //highest Global id of the elements in the partition
	long long size = BLOCK_SIZE(id, p, n); //number of the elements in the partition

	double* D = (double*)calloc(size * k, sizeof(double)); //matrix to store distances from centroids
	int* Cl = (int*)calloc(size, sizeof(double)); //vector to store cluster id of each element
	int* Cl_n = (int*)calloc(k, sizeof(int)); //vector to store number of elements in cluster in a processor
	int* Cl_n_total = (int*)calloc(k, sizeof(int)); //vector to store total number of elements in cluster for all processor
	double* Cn = (double*)calloc(k * m, sizeof(double)); //matrix to store k centroids
	double* A_local = (double*)malloc(size * m * sizeof(double)); //matrix to store all features for the local elements
	double* C_sum = (double*)calloc(k * m, sizeof(double));


	int* counts = (int*)calloc(p, sizeof(int)); // Declare the counts (= element*feature) in each partition
	int* displacements = (int*)calloc(p, sizeof(int)); // Declare the displacements in each partition
	for (int i = 0; i < p; i++) {
		counts[i] = BLOCK_SIZE(i, p, n) * m;
		displacements[i] = BLOCK_LOW(i, p, n) * m;
	}

	if (id == 0) {
		cout << "n " << n << " m: " << m << " k: " << k << endl;
		double* A = (double*)malloc(n * m * sizeof(double)); //matrix to store all features for all the elements
		srand(time(NULL));
		RandMatrixGen(A, n, m); // we create the data using a random generator here. If we want to use real dataset we need to modify this part
		initializeCn(A, n, m, Cn, k);
		//send the element blocks to proccessors using scatterv
		MPI_CALL(MPI_Scatterv(A, counts, displacements, MPI_DOUBLE, A_local, counts[id], MPI_DOUBLE, 0, MPI_COMM_WORLD));
		//free(A);
	}
	else {
		MPI_CALL(MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, A_local, counts[id], MPI_DOUBLE, 0, MPI_COMM_WORLD));
	}

	MPI_CALL(MPI_Bcast(Cn, k * m, MPI_DOUBLE, 0, MPI_COMM_WORLD)); //broadcast initial centroids to all processors


	//Computation Stage
	for (int i = 0; i < itr; i++) {
		calculateDistance(A_local, Cn, D, size, m, k);
		decideClusterID(D, size, k, Cl);
		/*free(C_sum);
		free(Cl_n);
		free(Cn);
		free(Cl_n_total);*/

		//Cn = (double*)calloc(k * m, sizeof(double)); //reinitialize all values of Cn to 0
		//Cl_n_total = (int*)calloc(k, sizeof(int));
		computeCentroidSum(A_local, size, m, Cl, C_sum, Cl_n);
		MPI_CALL(MPI_Allreduce(C_sum, Cn, k * m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
		MPI_CALL(MPI_Allreduce(Cl_n, Cl_n_total, k, MPI_INT, MPI_SUM, MPI_COMM_WORLD));
		computeNewCentroid(m, k, Cn, Cl_n_total);
		C_sum = (double*)calloc(k * m, sizeof(double)); //reinitialize all values of Cn to 0
		Cl_n = (int*)calloc(k, sizeof(int));
		D = (double*)calloc(size * k, sizeof(double)); //reinitializing distance
		//MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time += MPI_Wtime();
	if (!id) {
		printf("Total elapsed time: %10.6f\n", elapsed_time);
	}

	MPI_CALL(MPI_Finalize());
	return 0;
}