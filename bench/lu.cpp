#include "./utility/rt.h"

int N;

#define A_OFFSET 0

void lu_trace(double* A) {

	int i, j, k;
    vector<int> idx;

	for (i = 0; i < N; i++) {
		for (j = 0; j <i; j++) {
			for (k = 0; k < j; k++) {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
				A[i * N + j] -= A[i * N + k] * A[k * N + j];
				rtTmpAccess(A_OFFSET + i * N + k, 0, 0, idx);
				rtTmpAccess(A_OFFSET + k * N + j, 1, 0, idx);
				rtTmpAccess(A_OFFSET + i * N + j, 2, 0, idx);
				rtTmpAccess(A_OFFSET + i * N + j, 3, 0, idx);
			}
            idx.clear(); idx.push_back(i); idx.push_back(j);
			A[i * N + j] /= A[j * N + j];
			rtTmpAccess(A_OFFSET + j * N + j, 4, 0, idx);
			rtTmpAccess(A_OFFSET + i * N + j, 5, 0, idx);
			rtTmpAccess(A_OFFSET + i * N + j, 6, 0, idx);
		}
		for (j = i; j < N; j++) {
			for (k = 0; k < i; k++) {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
				A[i * N + j] -= A[i * N + k] * A[k * N + j];
				rtTmpAccess(A_OFFSET + i * N + k, 7, 0, idx);
				rtTmpAccess(A_OFFSET + k * N + j, 8, 0, idx);
				rtTmpAccess(A_OFFSET + i * N + j, 9, 0, idx);
				rtTmpAccess(A_OFFSET + i * N + j, 10, 0, idx);
			}
		}
	}
}

int main(int argc, char* argv[]) {

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
    
	double* A = (double *)malloc(N * N * sizeof(double));

	lu_trace(A);

    dumpRIHistogram();
    
    return 0;
}
