#include "./utility/rt.h"

int N;
int M;

#define A_OFFSET 0
#define B_OFFSET M * M
#define C_OFFSET M * M + N * M

void syr2k_trace(double* A, double* B, double* C, double alpha, double beta) {

	int i, j, k;
    vector<int> idx;
    
	for (i = 0; i < N; i++) {
		for (j = 0; j <= i; j++) {
            idx.clear(); idx.push_back(i); idx.push_back(j);
			C[i * N + j] *= beta;
			rtTmpAccess(C_OFFSET + i * N + j, 0, 0, idx);
			rtTmpAccess(C_OFFSET + i * N + j, 1, 0, idx);
		}
		for (k = 0; k < M; k++) {
			for (j = 0; j <= i; j++) {
                idx.clear(); idx.push_back(i); idx.push_back(k); idx.push_back(j);
				C[i * N + j] += A[j * M + k]*alpha*B[i * M + k] + B[j * M + k]*alpha*A[i * N + k];
				rtTmpAccess(A_OFFSET + j * M + k, 2, 1, idx);
				rtTmpAccess(B_OFFSET + i * M + k, 3, 2, idx);
				rtTmpAccess(B_OFFSET + j * M + k, 4, 2, idx);
				rtTmpAccess(A_OFFSET + i * N + k, 5, 1, idx);
				rtTmpAccess(C_OFFSET + i * N + j, 6, 0, idx);
				rtTmpAccess(C_OFFSET + i * N + j, 7, 0, idx);
			}
		}
	}
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
	double* B = (double *)malloc(N * M * sizeof(double));
	double* C = (double *)malloc(N * N * sizeof(double));
	double alpha = 0.2;
	double beta = 0.8;

	syr2k_trace(A, B, C, alpha, beta);

    dumpRIHistogram();
    
	return 0;
}
