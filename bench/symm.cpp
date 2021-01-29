#include "./utility/rt.h"

int M;
int N;

#define A_OFFSET 0
#define B_OFFSET M * M
#define C_OFFSET M * M + M * N

void symm_trace(double* A, double* B, double* C, double alpha, double beta) {
	
	int i, j, k;	
	double temp2;
    vector<int> idx;

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++ ) {
			temp2 = 0;
			for (k = 0; k < i; k++) {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
				C[k * N + j] += alpha*B[i * N + j] * A[i * M + k];
				temp2 += B[k * N + j] * A[i * M + k];
                rtTmpAccess(B_OFFSET + i * N + j, 0, 1, idx);
                rtTmpAccess(A_OFFSET + i * M + k, 1, 2, idx);
                rtTmpAccess(C_OFFSET + k * N + j, 2, 0, idx);
                rtTmpAccess(C_OFFSET + k * N + j, 3, 0, idx);
                rtTmpAccess(B_OFFSET + k * N + j, 4, 1, idx);
                rtTmpAccess(A_OFFSET + i * M + k, 5, 2, idx);
			}
            idx.clear(); idx.push_back(i); idx.push_back(j);
            
			C[i * N + j] = beta * C[i * N + j] + alpha*B[i * N + j] * A[i * M + i] + alpha * temp2;
            rtTmpAccess(C_OFFSET + i * N + j, 6, 0, idx);
            rtTmpAccess(B_OFFSET + i * N + j, 7, 1, idx);
            rtTmpAccess(A_OFFSET + i * M + i, 8, 2, idx);
            rtTmpAccess(C_OFFSET + i * N + j, 9, 0, idx);
		}
	}
	return;
}

int main(int argc, char* argv[]) {

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
    
    N = stoi(argv[1]);
    M = stoi(argv[2]);
    
	double* A = (double *)malloc(M * M * sizeof(double));
	double* B = (double *)malloc(M * N * sizeof(double));
	double* C = (double *)malloc(M * N * sizeof(double));
	double alpha = 0.2;
	double beta = 0.8;

	symm_trace(A, B, C, alpha, beta);

    dumpRIHistogram();
    
	return 0;
}
