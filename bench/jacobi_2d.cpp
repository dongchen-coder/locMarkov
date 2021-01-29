#include "./utility/rt.h"

int TSTEPS;
int N;

#define A_OFFSET 0
#define B_OFFSET N * N

void jacobi_2d_trace(double* A, double* B) {

	int t, i, j;
    vector<int> idx;

	for (t = 0; t < TSTEPS; t++) {
		for (i = 1; i < N - 1; i++) {
			for (j = 1; j < N - 1; j++) {
                idx.clear(); idx.push_back(t); idx.push_back(i); idx.push_back(j);
				B[i * N + j] = 0.2 * (A[i * N + j] + A[i * N + j-1] + A[i * N + 1+j] + A[(1+i) * N + j] + A[(i-1) * N + j]);
				rtTmpAccess(A_OFFSET + i * N + j, 0, 0, idx);
				rtTmpAccess(A_OFFSET + i * N + j-1, 1, 0, idx);
				rtTmpAccess(A_OFFSET + i * N + 1+j, 2, 0, idx);
				rtTmpAccess(A_OFFSET + (1+i) * N + j, 3, 0, idx);
				rtTmpAccess(A_OFFSET + (i-1) * N + j, 4, 0, idx);
				rtTmpAccess(B_OFFSET + i * N + j, 5, 1, idx);
			}
		}
		for (i = 1; i < N - 1; i++) {
			for (j = 1; j < N - 1; j++) {
                idx.clear(); idx.push_back(t); idx.push_back(i); idx.push_back(j);
				A[i * N + j] = 0.2 * (B[i * N + j] + B[i * N + j-1] + B[i * N + 1+j] + B[(1+i) * N + j] + B[(i-1) * N + j]);
				rtTmpAccess(B_OFFSET + i * N + j, 6, 1, idx);
				rtTmpAccess(B_OFFSET + i * N + j-1, 7, 1, idx);
				rtTmpAccess(B_OFFSET + i * N + 1+j, 8, 1, idx);
				rtTmpAccess(B_OFFSET + (1+i) * N + j, 9, 1, idx);
				rtTmpAccess(B_OFFSET + (i-1) * N + j, 10, 1, idx);
				rtTmpAccess(A_OFFSET + i * N + j, 11, 0, idx);
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
    TSTEPS = stoi(argv[2]);
    
	double* A = (double *)malloc(N * N * sizeof(double));
	double* B = (double *)malloc(N * N * sizeof(double));

	jacobi_2d_trace(A, B);

    dumpRIHistogram();
    
    return 0;
}
