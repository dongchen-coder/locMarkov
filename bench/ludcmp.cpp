#include "./utility/rt.h"

int N;

#define A_OFFSET 0
#define B_OFFSET N * N
#define Y_OFFSET N * N + N
#define X_OFFSET N * N + N + N

void ludcmp_trace(double* A, double* b, double* y, double* x) {

	int i, j, k;
	double w;
    vector<int> idx;

	for (i = 0; i < N; i++) {
		for (j = 0; j <i; j++) {
            idx.clear(); idx.push_back(i); idx.push_back(j);
			w = A[i * N + j];
			rtTmpAccess(A_OFFSET + i * N + j, 0, 0, idx);
			for (k = 0; k < j; k++) {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
				w -= A[i * N + k] * A[k * N + j];
				rtTmpAccess(A_OFFSET + i * N + k, 1, 0, idx);
				rtTmpAccess(A_OFFSET + k * N + j, 2, 0, idx);
			}
            
            idx.clear(); idx.push_back(i); idx.push_back(j);
			A[i * N + j] = w / A[j * N + j];
			rtTmpAccess(A_OFFSET + j * N + j, 3, 0, idx);
			rtTmpAccess(A_OFFSET + i * N + j, 4, 0, idx);
		}
		for (j = i; j < N; j++) {
            idx.clear(); idx.push_back(i); idx.push_back(j);
			w = A[i * N + j];
			rtTmpAccess(A_OFFSET + i * N + j, 5, 0, idx);
			for (k = 0; k < i; k++) {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
				w -= A[i * N + k] * A[k * N + j];
				rtTmpAccess(A_OFFSET + i * N + k, 6, 0, idx);
				rtTmpAccess(A_OFFSET + k * N + j, 7, 0, idx);
			}
            idx.clear(); idx.push_back(i); idx.push_back(j);
			A[i * N + j] = w;
			rtTmpAccess(A_OFFSET + i * N + j, 8, 0, idx);
		}
	}

	for (i = 0; i < N; i++) {
        idx.clear(); idx.push_back(i);
		w = b[i];
		rtTmpAccess(B_OFFSET + i, 9, 1, idx);
		for (j = 0; j < i; j++) {
            idx.clear(); idx.push_back(i); idx.push_back(j);
			w -= A[i * N + j] * y[j];
			rtTmpAccess(A_OFFSET + i * N + j, 10, 0, idx);
			rtTmpAccess(Y_OFFSET + j, 11, 2, idx);
		}
        idx.clear(); idx.push_back(i);
		y[i] = w;
		rtTmpAccess(Y_OFFSET + i, 12, 2, idx);
	}
	for (i = N-1; i >=0; i--) {
        idx.clear(); idx.push_back(i);
		w = y[i];
		rtTmpAccess(Y_OFFSET + i, 13, 2, idx);
		for (j = i+1; j < N; j++) {
            idx.clear(); idx.push_back(i); idx.push_back(j);
			w -= A[i * N + j] * x[j];
			rtTmpAccess(A_OFFSET + i * N + j, 14, 0, idx);
			rtTmpAccess(X_OFFSET + j, 15, 3, idx);
		}
        idx.clear(); idx.push_back(i);
		x[i] = w / A[i * N + i];
		rtTmpAccess(A_OFFSET + i * N + i, 16, 0, idx);
		rtTmpAccess(X_OFFSET + i, 17, 3, idx);
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
    
	double* A = (double *)malloc(N * N *sizeof(double));
	double* b = (double *)malloc(N * sizeof(double));
	double* y = (double *)malloc(N * sizeof(double));
	double* x = (double *)malloc(N * sizeof(double));

	ludcmp_trace(A, b, y, x);

    dumpRIHistogram();
    
    return 0;
}
