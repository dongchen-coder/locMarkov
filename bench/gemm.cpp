#include "./utility/rt.h"

int NI;
int NJ;
int NK;

#define A_OFFSET 0
#define B_OFFSET NI * NK
#define C_OFFSET NI * NK + NK * NJ

void gemm_trace(double alpha, double beta, double* A, double* B, double* C) {
    int i,j,k;
    vector<int> idx;
    
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NJ; j++)
        {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            C[i * NJ + j] *= beta;
			rtTmpAccess(C_OFFSET + i * NJ + j, 0, 0, idx);
			rtTmpAccess(C_OFFSET + i * NJ + j, 1, 0, idx);

            for (k = 0; k < NK; k++)
            {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
                C[i * NJ + j] += alpha * A[i * NK + k] * B[k * NJ + j];
                rtTmpAccess(A_OFFSET + i * NK + k, 2, 1, idx);
                rtTmpAccess(B_OFFSET + k * NJ + j, 3, 2, idx);
				rtTmpAccess(C_OFFSET + i * NJ + j, 4, 0, idx);
				rtTmpAccess(C_OFFSET + i * NJ + j, 5, 0, idx);
            }
        }
    }
    return;
}

int main(int argc, char* argv[]) {
    
    if (argc != 4) {
        cout << "This benchmark needs 3 loop bounds" << endl;
        return 0;
    }
    for (int i = 1; i < argc; i++) {
        if (!isdigit(argv[i][0])) {
            cout << "arguments must be integer" << endl;
            return 0;
        }
    }
    
    NI = stoi(argv[1]);
    NJ = stoi(argv[2]);
    NK = stoi(argv[3]);
    
    double* A = (double*)malloc(NI * NK * sizeof(double));
    double* B = (double*)malloc(NK * NJ * sizeof(double));
    double* C = (double*)malloc(NI * NJ * sizeof(double));

    for (int i = 0; i < NI*NK; ++i)
    {
        A[i] = i % 256;
    }

    for (int i = 0; i < NK * NJ; ++i)
    {
        B[i] = i % 48;
    }

    double alpha = 1.0;
    double beta = 1.5;

    gemm_trace(alpha, beta, A, B, C);

    dumpRIHistogram();
    
    return 0;
}
