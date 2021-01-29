#include "./utility/rt.h"

int NR;
int NQ;
int NP;

#define SUM_OFFSET 0
#define A_OFFSET NP
#define C4_OFFSET NP + NR * NQ * NP

void doitgen_trace(double* sum, double* A, double* C4) {

	int r, q, p, s;
    vector<int> idx;

	for (r = 0; r < NR; r++) {
		for (q = 0; q < NQ; q++)  {
			for (p = 0; p < NP; p++)  {
                idx.clear(); idx.push_back(r); idx.push_back(q); idx.push_back(p);
				sum[p] = 0.0;
				rtTmpAccess(SUM_OFFSET + p, 0, 0, idx);
				for (s = 0; s < NP; s++) {
                    idx.clear(); idx.push_back(r); idx.push_back(q); idx.push_back(p); idx.push_back(s);
					sum[p] += A[r * NQ * NP + q * NP + s] * C4[s * NP + p];
					rtTmpAccess(A_OFFSET + r * NQ * NP + q * NP + s, 1, 1, idx);
					rtTmpAccess(C4_OFFSET + s * NP + p, 2, 2, idx);
					rtTmpAccess(SUM_OFFSET + p, 3, 0, idx);
					rtTmpAccess(SUM_OFFSET + p, 4, 0, idx);
				}
			}
			for (p = 0; p < NP; p++) {
                idx.clear(); idx.push_back(r); idx.push_back(q); idx.push_back(p);
				A[r * NQ * NP + q * NP + p] = sum[p];
				rtTmpAccess(SUM_OFFSET + p, 5, 0, idx);
				rtTmpAccess(A_OFFSET + r * NQ * NP + q * NP + p, 6, 1, idx);
			}
    	}
	}

}

int main(int argc, char* argv[]) {
	
    if (argc != 4) {
        cout << "This benchmark needs 3 loop bounds" << endl;
        return 0;
    }
    for (int i = 1; i < argc; i++) {
        if (!isdigit(argv[i][0])) {
            cout << "arguments must be integer" << endl;
            return 0;
        }
    }
    
    NR = stoi(argv[1]);
    NQ = stoi(argv[2]);
    NP = stoi(argv[3]);
    
	double * sum = (double *) malloc(NP * sizeof(double));
	double * A = (double *) malloc(NR * NQ * NP * sizeof(double));
	double * C4 = (double *) malloc(NP * NP);
	
	doitgen_trace(sum, A, C4);

    dumpRIHistogram();
    
	return 0;
}


