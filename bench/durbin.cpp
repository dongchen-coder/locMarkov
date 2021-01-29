#include "./utility/rt.h"

int N;

#define Y_OFFSET 0
#define R_OFFSET N
#define Z_OFFSET N + N

void durbin_trace(double *y, double* r, double* z) {

	int k, i;
	double alpha, beta, sum;
    vector<int> idx;

    /*
	y[0] = -r[0];
	beta = 1.0;
	alpha = -r[0];
    */
	for (k = 1; k < N; k++) {
		beta = (1-alpha*alpha)*beta;
		sum = 0.0;
		for (i=0; i<k; i++) {
            idx.clear(); idx.push_back(k); idx.push_back(i);
			sum += r[k-i-1]*y[i];
			rtTmpAccess(R_OFFSET + k-i-1, 0, 0, idx);
			rtTmpAccess(Y_OFFSET + i, 1, 1, idx);
		}
        idx.clear(); idx.push_back(k);
		alpha = - (r[k] + sum)/beta;
		rtTmpAccess(R_OFFSET + k, 2, 0, idx);

		for (i=0; i<k; i++) {
            idx.clear(); idx.push_back(k); idx.push_back(i);
			z[i] = y[i] + alpha*y[k-i-1];
			rtTmpAccess(Y_OFFSET + i, 3, 1, idx);
			rtTmpAccess(Y_OFFSET + k-i-1, 4, 1, idx);
			rtTmpAccess(Z_OFFSET + i, 5, 2, idx);
		}
		for (i=0; i<k; i++) {
            idx.clear(); idx.push_back(k); idx.push_back(i);
			y[i] = z[i];
			rtTmpAccess(Z_OFFSET + i, 6, 2, idx);
			rtTmpAccess(Y_OFFSET + i, 7, 1, idx);
		}
        
        idx.clear(); idx.push_back(k);
		y[k] = alpha;
		rtTmpAccess(Y_OFFSET + k, 8, 1, idx);
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
    
	double * y = (double *) malloc(N * sizeof(double));
	double * r = (double *) malloc(N * sizeof(double));
	double * z = (double *) malloc(N * sizeof(double));

	durbin_trace(y, r, z);

    dumpRIHistogram();
    
	return 0;
}

