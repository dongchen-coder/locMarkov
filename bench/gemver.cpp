#include "./utility/rt.h"

int N;

#define A_OFFSET 0
#define u1_OFFSET N * N  
#define v1_OFFSET N * N + N
#define u2_OFFSET N * N + N + N
#define v2_OFFSET N * N + N + N + N
#define W_OFFSET N * N + N + N + N + N 
#define X_OFFSET N * N + N + N + N + N + N
#define Y_OFFSET N * N + N + N + N + N + N + N
#define Z_OFFSET N * N + N + N + N + N + N + N + N

void gemver_trace(double alpha, double beta, double* A, double* u1, double* v1, double* u2, double* v2, double* w, double* x, double* y, double* z)
{
	int i,j;
    vector<int> idx;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            A[i * N + j] = A[i * N + j] + u1[i] * v1[j] + u2[i] * v2[j];
			rtTmpAccess(A_OFFSET + i * N + j, 0, 0, idx);
			rtTmpAccess(u1_OFFSET + i, 1, 1, idx);
			rtTmpAccess(v1_OFFSET + j, 2, 1, idx);
			rtTmpAccess(u2_OFFSET + i, 3, 2, idx);
			rtTmpAccess(v2_OFFSET + j, 4, 2, idx);
			rtTmpAccess(A_OFFSET + i * N + j, 5, 0, idx);
       }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            x[i] = x[i] + beta * A[j * N + i] * y[j];
			rtTmpAccess(X_OFFSET + i, 6, 3, idx);
			rtTmpAccess(A_OFFSET + j * N + i, 7, 0, idx);
			rtTmpAccess(Y_OFFSET + j, 8, 4, idx);
			rtTmpAccess(X_OFFSET + i, 9, 3, idx);
        }
    }

    for (i = 0; i < N; i++)
    {
        idx.clear(); idx.push_back(i);
        x[i] = x[i] + z[i];
		rtTmpAccess(X_OFFSET + i, 10, 3, idx);
		rtTmpAccess(Z_OFFSET + i, 11, 5, idx);
		rtTmpAccess(X_OFFSET + i, 12, 3, idx);
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            w[i] = w[i] +  alpha * A[i * N + j] * x[j];
			rtTmpAccess(W_OFFSET + i, 13, 6, idx);
			rtTmpAccess(A_OFFSET + i * N + j, 14, 0, idx);
			rtTmpAccess(X_OFFSET + j, 15, 3, idx);
			rtTmpAccess(W_OFFSET + i, 16, 6, idx);
        }
    }
}

int main(int argc, char const *argv[])
{
    if (argc != 2) {
        cout << "This benchmark needs 1 loop bounds" << endl;
        return 0;
    }
    for (int i = 1; i < argc; i++) {
        if (!isdigit(argv[i][0])) {
            cout << "arguments must be integer" << endl;
            return 0;
        }
    }
    
    N = stoi(argv[1]);
    
    double* u1 = (double*)malloc( N * sizeof(double));
    double* v1 = (double*)malloc( N * sizeof(double));
    double* u2 = (double*)malloc( N * sizeof(double));
    double* v2 = (double*)malloc( N * sizeof(double));
    double* A = (double*)malloc( (N*N) * sizeof(double));
    double* y = (double*)malloc( N * sizeof(double));
    double* x = (double*)malloc( N * sizeof(double));
    double* z = (double*)malloc( N * sizeof(double));
    double* w = (double*)malloc( N * sizeof(double));

    for (int i = 0; i < N; ++i)
    {
        u1[i] = i % 256;
        u2[i] = i % 32;
        v1[i] = i % 2;
        v2[i] = i % 512;
        y[i] = i % 4;
        z[i] = i % 64;
        x[i] = 0.0;
        w[i] = 0.0;
    }

    for (int i = 0; i < N*N; ++i)
    {
        A[i] = 0.0;
    }

    double alpha = 1.0;
    double beta = 1.5;

    gemver_trace(alpha, beta, A, u1, v1, u2, v2, w, x, y, z);
 	
    dumpRIHistogram();
    
    return 0;
}
