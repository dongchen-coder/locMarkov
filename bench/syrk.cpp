#include "./utility/rt.h"

int NI;
int NJ;

#define A_OFFSET 0
#define C_OFFSET NI * NJ

void syrk_trace(double alpha, double beta, double* A, double* C)
{
    int i, j, k;
    vector<int> idx;
    
    /*  C := alpha*A*A' + beta*C */
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NI; j++)
        {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            C[i * NI + j] = C[i * NI + j] * beta;
			rtTmpAccess(C_OFFSET + i * NI + j, 0, 0, idx);
			rtTmpAccess(C_OFFSET + i * NI + j, 1, 0, idx);
        }
    }
    
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NI; j++)
        {
            for (k = 0; k < NJ; k++)
            {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
                C[i * NI + j] = C[i * NI + j] + alpha * A[i * NJ + k] * A[j * NJ + k];
				rtTmpAccess(A_OFFSET + i * NJ + k, 2, 1, idx);
				rtTmpAccess(A_OFFSET + j * NJ + k, 3, 1, idx);
				rtTmpAccess(C_OFFSET + i * NI + j, 4, 0, idx);
				rtTmpAccess(C_OFFSET + i * NI + j, 5, 0, idx);
            }
        }
    }
    return;
}

int main(int argc, char const *argv[])
{
    
    if (argc != 3) {
        cout << "This benchmark needs 2 loop bounds" << endl;
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
    
    double* A = (double*)malloc( (NI*NJ) * sizeof(double));
    double* C = (double*)malloc( (NI*NI) * sizeof(double));

    double alpha = 0.0;
    double beta = 1.5;

    syrk_trace(alpha, beta, A, C);

    dumpRIHistogram();
    
    return 0;
}
