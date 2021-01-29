#include "./utility/rt.h"

int TSTEPS;
int N;

#define A_OFFSET 0

void seidel_2d_trace(double* A) {

	int t, i, j;
    vector<int> idx;

	for (t = 0; t <= TSTEPS - 1; t++)
		for (i = 1; i<= N - 2; i++)
			for (j = 1; j <= N - 2; j++) {
                idx.clear(); idx.push_back(t); idx.push_back(i); idx.push_back(j);
				A[i * N + j] = (A[(i-1) * N + j-1] + A[(i-1) * N + j] + A[(i-1) * N + j+1] + A[i * N + j-1] + A[i * N + j] + A[i * N + j+1] + A[(i+1) * N + j-1] + A[(i+1) * N + j] + A[(i+1) * N + j+1]) / 9.0;
				rtTmpAccess(A_OFFSET + (i-1) * N + j-1, 0, 0, idx);
				rtTmpAccess(A_OFFSET + (i-1) * N + j, 1, 0, idx);
				rtTmpAccess(A_OFFSET + (i-1) * N + j+1, 2, 0, idx);
				rtTmpAccess(A_OFFSET + i * N + j-1, 3, 0, idx);
				rtTmpAccess(A_OFFSET + i * N + j, 4, 0, idx);
				rtTmpAccess(A_OFFSET + i * N + j+1, 5, 0, idx);
				rtTmpAccess(A_OFFSET + (i+1) * N + j-1, 6, 0, idx);
				rtTmpAccess(A_OFFSET + (i+1) * N + j, 7, 0, idx);
				rtTmpAccess(A_OFFSET + (i+1) * N + j+1, 8, 0, idx);
				rtTmpAccess(A_OFFSET + i * N + j, 9, 0, idx);
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

	seidel_2d_trace(A);

    dumpRIHistogram();
    
	return 0;
}
