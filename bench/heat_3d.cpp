#include "./utility/rt.h"

int TSTEPS;
int N;

#define A_OFFSET 0
#define B_OFFSET N * N * N

void heat_3d_trace(double* B, double* A) {

	int t, i, j, k;
    vector<int> idx;

	for (t = 1; t <= TSTEPS; t++) {
		for (i = 1; i < N-1; i++) {
			for (j = 1; j < N-1; j++) {
				for (k = 1; k < N-1; k++) {
                    idx.clear(); idx.push_back(t); idx.push_back(i); idx.push_back(j); idx.push_back(k);
					B[i * N * N + j * N + k] =   0.125 * (A[(i+1) * N * N + j * N + k] - 2.0 * A[i * N * N + j * N + k] + A[(i-1) * N * N + j * N + k])
                                 + 0.125 * (A[i * N * N + (j+1) * N + k] - 2.0 * A[i * N * N + j * N + k] + A[i * N * N + (j-1) * N + k])
                                 + 0.125 * (A[i * N * N + j * N + k+1] - 2.0 * A[i * N * N + j * N + k] + A[i * N * N + j * N + k-1])
                                 + A[i * N * N + j * N + k];
                	rtTmpAccess(A_OFFSET + (i+1) * N * N + j * N + k, 0, 0, idx);
					rtTmpAccess(A_OFFSET + i * N * N + j * N + k, 1, 0, idx);
					rtTmpAccess(A_OFFSET + (i-1) * N * N + j * N + k, 2, 0, idx);
					rtTmpAccess(A_OFFSET + i * N * N + (j+1) * N + k, 3, 0, idx);
					rtTmpAccess(A_OFFSET + i * N * N + j * N + k, 4, 0, idx);
					rtTmpAccess(A_OFFSET + i * N * N + (j-1) * N + k, 5, 0, idx);
					rtTmpAccess(A_OFFSET + i * N * N + j * N + k+1, 6, 0, idx);
					rtTmpAccess(A_OFFSET + i * N * N + j * N + k, 7, 0, idx);
					rtTmpAccess(A_OFFSET + i * N * N + j * N + k-1, 8, 0, idx);
					rtTmpAccess(A_OFFSET + i * N * N + j * N + k, 9, 0, idx);
					rtTmpAccess(B_OFFSET + i * N * N + j * N + k, 10, 1, idx);
				}
            }
        }
        for (i = 1; i < N-1; i++) {
           for (j = 1; j < N-1; j++) {
               for (k = 1; k < N-1; k++) {
                    idx.clear(); idx.push_back(t); idx.push_back(i); idx.push_back(j); idx.push_back(k);
					A[i * N * N + j * N + k] =   0.125 * (B[(i+1) * N * N + j * N + k] - 2.0 * B[i * N * N + j * N + k] + B[(i-1) * N * N + j * N + k])
                                + 0.125 * (B[i * N * N + (j+1) * N + k] - 2.0 * B[i * N * N + j * N + k] + B[i * N * N + (j-1) * N + k])
                                + 0.125 * (B[i * N * N + j * N + k+1] - 2.0 * B[i * N * N + j * N + k] + B[i * N * N + j * N + k-1])
                                + B[i * N * N + j * N + k];
					rtTmpAccess(B_OFFSET + (i+1) * N * N + j * N + k, 11, 1, idx);
					rtTmpAccess(B_OFFSET + i * N * N + j * N + k, 12, 1, idx);
					rtTmpAccess(B_OFFSET + (i-1) * N * N + j * N + k, 13, 1, idx);
					rtTmpAccess(B_OFFSET + i * N * N + (j+1) * N + k, 14, 1, idx);
					rtTmpAccess(B_OFFSET + i * N * N + j * N + k, 15, 1, idx);
					rtTmpAccess(B_OFFSET + i * N * N + (j-1) * N + k, 16, 1, idx);
					rtTmpAccess(B_OFFSET + i * N * N + j * N + k+1, 17, 1, idx);
					rtTmpAccess(B_OFFSET + i * N * N + j * N + k, 18, 1, idx);
					rtTmpAccess(B_OFFSET + i * N * N + j * N + k-1, 19, 1, idx);
					rtTmpAccess(B_OFFSET + i * N * N + j * N + k, 20, 1, idx);
					rtTmpAccess(A_OFFSET + i * N * N + j * N + k, 21, 0, idx);
               }
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

	double* A = (double *)malloc(N * N * N * sizeof(double));
	double* B = (double *)malloc(N * N * N * sizeof(double));

	heat_3d_trace(B, A);

    dumpRIHistogram();
    
	return 0;
}
