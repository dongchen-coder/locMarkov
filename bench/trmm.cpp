#include "./utility/rt.h"

int M;
int N;

#define A_OFFSET 0
#define B_OFFSET M * M

void trmm_trace(double* A, double* B, double alpha) {

	int i, j, k;
    vector<int> idx;

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			for (k = i+1; k < M; k++) {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
				B[i * N + j] += A[k * M + i] * B[k * N + j];
				rtTmpAccess(A_OFFSET + k * M + i, 0, 0, idx);
				rtTmpAccess(B_OFFSET + k * N + j, 1, 1, idx);
				rtTmpAccess(B_OFFSET + i * N + j, 2, 1, idx);
				rtTmpAccess(B_OFFSET + i * N + j, 3, 1, idx);
			}
            idx.clear(); idx.push_back(i); idx.push_back(j);
			B[i * N + j] = alpha * B[i * N + j];
     		rtTmpAccess(B_OFFSET + i * N + j, 4, 1, idx);
			rtTmpAccess(B_OFFSET + i * N + j, 5, 1, idx);
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
    
    M = stoi(argv[1]);
    N = stoi(argv[2]);
    
	double* A = (double *)malloc(M * M * sizeof(double));
	double* B = (double *)malloc(M * N * sizeof(double));
	double alpha = 0.2;

	trmm_trace(A, B, alpha);
 
    dumpRIHistogram();
    
    return 0;
}
