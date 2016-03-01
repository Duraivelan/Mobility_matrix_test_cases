#include <iostream>
#include <random>

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

void inverse(double* A, int N)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete IPIV;
    delete WORK;
}

using namespace std;
int main(){

    double A [3*3] = {10,2,3,40,5,6,7,8,10};
    
    cout<< A[0]<<'\t'<<A[1]<<endl;
    cout<<A[2]<<'\t'<<A[3]<<endl;

    inverse(A, 3);

    cout<< A[0]<<'\t'<<A[1]<<endl;
    cout<<A[2]<<'\t'<<A[3]<<endl;


    return 0;
}
