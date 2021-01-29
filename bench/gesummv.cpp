#include "./utility/rt.h"

int N;

#define A_OFFSET 0
#define B_OFFSET N * N
#define TMP_OFFSET N * N + N * N
#define X_OFFSET N * N + N * N + N
#define Y_OFFSET N * N + N * N + N + N

void gesummv_trace(double alpha, double beta, double* A, double* B, double* tmp, double* x, double* y)
{
	int i, j;    
    vector<int> idx;
    
	for (i = 0; i < N; i++)
    {
        idx.clear(); idx.push_back(i);
        tmp[i] = 0;
        y[i] = 0;

		rtTmpAccess(TMP_OFFSET + i, 0, 0, idx);
		rtTmpAccess(Y_OFFSET + i, 1, 1, idx);
	
        for (j = 0; j < N; j++)
        {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            tmp[i] = A[i * N + j] * x[j] + tmp[i];
            y[i] = B[i * N + j] * x[j] + y[i];
        	rtTmpAccess(A_OFFSET + i * N + j, 2, 2, idx);
			rtTmpAccess(X_OFFSET + j, 3, 3, idx);
			rtTmpAccess(TMP_OFFSET + i, 4, 0, idx);
			rtTmpAccess(TMP_OFFSET + i, 5, 0, idx);
			rtTmpAccess(B_OFFSET + i * N + j, 6, 4, idx);
            rtTmpAccess(X_OFFSET + j, 7, 3, idx);
            rtTmpAccess(Y_OFFSET + i, 8, 1, idx);
            rtTmpAccess(Y_OFFSET + i, 9, 1, idx);
		}

        idx.clear(); idx.push_back(i);
        y[i] = alpha * tmp[i] + beta * y[i];
        
		rtTmpAccess(TMP_OFFSET + i, 10, 0, idx);
		rtTmpAccess(Y_OFFSET + i, 11, 1, idx);
		rtTmpAccess(Y_OFFSET + i, 12, 1, idx);
    }

	return;
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
    
    double* A = (double*)malloc( (N*N) * sizeof(double));
    double* x = (double*)malloc( N * sizeof(double));
    double* tmp = (double*)malloc( N * sizeof(double));
    double* B = (double*)malloc( (N*N) * sizeof(double));
    double* y = (double*)malloc( N * sizeof(double));

    for (int i = 0; i < N; ++i)
    {
        x[i] = i % 256;
    }

    for (int i = 0; i < N*N; ++i)
    {
        A[i] = i / 10;
        B[i] = i / 25;
    }

    double alpha = 1.0;
    double beta = 1.5;

    gesummv_trace(alpha, beta, A, B, tmp, x, y);
  
    dumpRIHistogram();
    
    return 0;
}
