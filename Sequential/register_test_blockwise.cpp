#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <vector>
#include <cmath>

using namespace std;

class Timer {
public:
    void start() {
        err = clock_gettime(CLOCK_REALTIME, &start_time);
    }
    void end() {
        err = clock_gettime(CLOCK_REALTIME, &end_time);
    }
    double get() {
        return (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec) / (double)1000000000;
    }
private:
    int err = 0;
    struct timespec start_time, end_time;
};

int RandMatrixGen(double* M, int n, int m) {

    for (int i = 0; i < n * m; i++)
        //M[i] = 2.0 * rand() / RAND_MAX - 1.0;
        M[i] = rand() % 100;
    return 0;
}

//int RandCentroidGen(double* M, int n, int m) {
//    srand(time(NULL));
//
//}

void calculateDistance(double* A, double* Cn, double* D, int n, int m, int k) {

    for (int l = 0; l < k; l+=2) {
        for (int i = 0; i < n; i+=2) {
          register double dis11 = D[i*k+l];
	  register double dis12 = D[i*k+l+1];
	  register double dis21 = D[(i+1)*k+l];
	  register double dis22 = D[(i+1)*k+l+1];

          register int r_A = i*m;
          register int r_Cn = l*m;


            for (int j = 0; j < m; j++) {

		register double A11 = A[r_A+j];
          	register double A21= A[r_A+m+j];
//          	register double A13= A[r_A +j+ 2];
//          	register double A14= A[r_A +j+ 3];

          	register double Cn11=Cn[r_Cn+j];
         	register double Cn21=Cn[r_Cn+m+j];
//          	register double Cn13=Cn[r_Cn+j + 2];
//          	register double Cn14=Cn[r_Cn+j + 3];
    
                dis11 += (A11- Cn11)*(A11- Cn11);
	        dis12 += (A11- Cn21)*(A11- Cn21);
		dis21 += (A21- Cn11)*(A21- Cn11);
		dis22 += (A21- Cn21)*(A21- Cn21);
            }

            D[i * k + l] = sqrt(dis11);
	    D[(i+1)*k+l] = sqrt(dis21);
	    D[i*k+l+1] = sqrt (dis12);
	    D[(i+1)*k+l+1] = sqrt (dis22);
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

void decideClusterID(double* D, int n,  int k, int* Cl) {
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

void printClusterID(int* Cl, int n) {
    for (int i = 0; i < n; i++) {
            cout << "point: " <<i << "clusterID: " << Cl[i] << endl;
    }
}



void computeNewCentroid(double* A, int n, int m, int* Cl, int k, double* Cn, int* Cl_n) {
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

void initializeCn(double* A, int n, int m,  double* Cn, int k) {
    for (int i = 0; i < k; i++) {
        int e = rand() % n;
        for (int j = 0; j < m; j++) {
            Cn[i * m + j] = A[e * m + j];
        }
    }

}


/// <summary>
/// arg 1 : n
/// arg 2 : m
/// arg 3 : k
/// arg 4 : number of itr
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
/// <returns></returns>
int main(int argc, char** argv) {
    Timer timer;
    int n = atoi(argv[1]); //number of elements
    int m = atoi(argv[2]); //number of features
    int k = atoi(argv[3]); //number of clusters
    int itr = atoi(argv[4]); //number of iterations

    cout << "n " << n << " m: " << m << " k: " << k <<"\n";
    double* A = (double*)malloc(n * m * sizeof(double)); //matrix to store all features for all the elements
    srand(time(NULL));
    RandMatrixGen(A, n, m);
//    cout << "printing A: \n";
//    printMatrix(A, n, m);
    double* D = (double*)calloc(n * k, sizeof(double)); //matrix to store distances from centroids
    int* Cl = (int*)calloc(n, sizeof(double)); //vector to store cluster id of each element
    int* Cl_n = (int*)calloc(n, sizeof(double)); //vector to store number of elements in cluster
    double* Cn = (double*)calloc(k * m, sizeof(double)); //matrix to store k centroids

    //RandMatrixGen(Cn, k, m);
    initializeCn(A, n, m, Cn, k);
//    cout << "printing Cn: \n";
//    printMatrix(Cn, k, m);

    timer.start();
    for (int i = 0; i < itr; i++) {
        calculateDistance(A, Cn, D, n, m, k);

//	cout << "itr:" << i << "Printing D\n";
//        printMatrix(D, n, k);


        decideClusterID(D, n, k, Cl);

//        cout << "printing Cluster IDs: \n";
//        printClusterID(Cl, n);

        Cn = (double*)calloc(k * m, sizeof(double)); //reinitialize all values of Cn to 0
        Cl_n = (int*)calloc(n, sizeof(double));
	D = (double*)calloc(n * k, sizeof(double));
        computeNewCentroid(A, n, m, Cl, k, Cn, Cl_n);

//        cout << "printing Cn: \n";
//        printMatrix(Cn, k, m);
    }

    timer.end();


    //cout << "printing Cluster IDs: \n";
    /*printClusterID(Cl, n);*/





    std::cout << " execution time (in sec): " << timer.get() << ", ";
}

