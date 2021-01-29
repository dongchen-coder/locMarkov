#include "./utility/rt.h"

int TSTEPS;
int N;

#define A_OFFSET 0
#define B_OFFSET N

void jacobi_1d_trace(double* A, double* B) {

	int t, i;
    vector<int> idx;

	for (t = 0; t < TSTEPS; t++) {
		for (i = 1; i < N - 1; i++) {
            idx.clear(); idx.push_back(t); idx.push_back(i);
			B[i] = 0.33333 * (A[i-1] + A[i] + A[i + 1]);
			rtTmpAccess(A_OFFSET + i - 1, 0, 0, idx);
			rtTmpAccess(A_OFFSET + i, 1, 0, idx);
			rtTmpAccess(A_OFFSET + i + 1, 2, 0, idx);
			rtTmpAccess(B_OFFSET + i, 3, 1, idx);
		}
		for (i = 1; i < N - 1; i++) {
            idx.clear(); idx.push_back(t); idx.push_back(i);
			A[i] = 0.33333 * (B[i-1] + B[i] + B[i + 1]);
			rtTmpAccess(B_OFFSET + i - 1, 4, 1, idx);
			rtTmpAccess(B_OFFSET + i, 5, 1, idx);
			rtTmpAccess(B_OFFSET + i + 1, 6, 1, idx);
			rtTmpAccess(A_OFFSET + i, 7, 0, idx);
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
    
	double* A = (double *)malloc(N * sizeof(double));
	double* B = (double *)malloc(N * sizeof(double));

	jacobi_1d_trace(A, B);

    dumpRIHistogram();
    
    return 0;
}
