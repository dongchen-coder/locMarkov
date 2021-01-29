#include "./utility/rt.h"
#include <math.h>

int N;
int M;

#define A_OFFSET 0
#define R_OFFSET N * M
#define Q_OFFSET N * M + N * N

void gramschmidt_trace(double* A, double* R, double* Q) {

	int k, i, j;
	double nrm;
    vector<int> idx;

	for (k = 0; k < N; k++) {
		nrm = 0.0;
		for (i = 0; i < M; i++) {
            idx.clear(); idx.push_back(k); idx.push_back(i);
			nrm += A[i * N + k] * A[i * N + k];
			rtTmpAccess(A_OFFSET + i * N + k, 0, 0, idx);
			rtTmpAccess(A_OFFSET + i * N + k, 1, 0, idx);
		}
        idx.clear(); idx.push_back(k);
		R[k * N + k] = sqrt(nrm);
		rtTmpAccess(R_OFFSET + k * N + k, 2, 1, idx);
		for (i = 0; i < M; i++) {
            idx.clear(); idx.push_back(k); idx.push_back(i);
			Q[i * N + k] = A[i * N + k] / R[k * N + k];
			rtTmpAccess(A_OFFSET + i * N + k, 3, 0, idx);
			rtTmpAccess(R_OFFSET + k * N + k, 4, 1, idx);
			rtTmpAccess(Q_OFFSET + i * N + k, 5, 2, idx);
		}
		for (j = k + 1; j < N; j++) {
            idx.clear(); idx.push_back(k); idx.push_back(j);
			R[k * N + j] = 0.0;
			rtTmpAccess(R_OFFSET + k * N + j, 6, 1, idx);
			for (i = 0; i < M; i++) {
                idx.clear(); idx.push_back(k); idx.push_back(j); idx.push_back(i);
				R[k * N + j] += Q[i * N + k] * A[i * N + j];
				rtTmpAccess(Q_OFFSET + i * N + k, 7, 2, idx);
				rtTmpAccess(A_OFFSET + i * N + j, 8, 0, idx);
				rtTmpAccess(R_OFFSET + k * N + j, 9, 1, idx);
                rtTmpAccess(R_OFFSET + k * N + j, 10, 1, idx);
			}
			for (i = 0; i < M; i++) {
                idx.clear(); idx.push_back(k); idx.push_back(j); idx.push_back(i);
        		A[i * N + j] = A[i * N + j] - Q[i * N + k] * R[k * N + j];
				rtTmpAccess(A_OFFSET + i * N + j, 11, 0, idx);
				rtTmpAccess(Q_OFFSET + i * N + k, 12, 2, idx);
				rtTmpAccess(R_OFFSET + k * N + j, 13, 1, idx);
				rtTmpAccess(A_OFFSET + i * N + j, 14, 0, idx);
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
    
    M = stoi(argv[1]);
    N = stoi(argv[2]);
    
	double * A = (double *)malloc(N * M * sizeof(double));
	double * R = (double *)malloc(N * N * sizeof(double));
	double * Q = (double *)malloc(M * N * sizeof(double));

	gramschmidt_trace(A, R, Q);

    dumpRIHistogram();
    
	return 0;
}
