#include "./utility/rt.h"

int NI;
int NJ;
int NK;
int NL;

#define TMP_OFFSET 0
#define A_OFFSET NI * NJ
#define B_OFFSET NI * NJ + NI * NK
#define D_OFFSET NI * NJ + NI * NK + NJ * NK
#define C_OFFSET NI * NJ + NI * NK + NJ * NK + NI * NL

void mm2_trace(double* tmp, double* A, double* B, double* C, double* D, double alpha, double beta) {

    int i, j, k;
    vector<int> idx;
    
    for (i = 0; i < NI; i++) {
        for (j = 0; j < NJ; j++) {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            
            tmp[i * NJ + j] = 0.0;
            rtTmpAccess(TMP_OFFSET + i * NJ + j, 0, 0, idx);
            for (k = 0; k < NK; ++k) {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
                
                tmp[i * NJ + j] += alpha * A[i * NK + k] * B[k * NJ + j];
				rtTmpAccess(A_OFFSET + i * NK + k, 1, 1, idx);
				rtTmpAccess(B_OFFSET + k * NJ + j, 2, 2, idx);
				rtTmpAccess(TMP_OFFSET + i * NJ + j, 3, 0, idx);
				rtTmpAccess(TMP_OFFSET + i * NJ + j, 4, 0, idx);
            }
        }
    }
    for (i = 0; i < NI; i++) {
        for (j = 0; j < NL; j++) {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            
            D[i * NL + j] *= beta;
			rtTmpAccess(D_OFFSET + i * NL + j, 5, 3, idx);
            rtTmpAccess(D_OFFSET + i * NL + j, 6, 3, idx);
            for (k = 0; k < NJ; ++k) {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
                
                D[i * NL + j] += tmp[i * NJ + k] * C[k * NL + j];
            	rtTmpAccess(TMP_OFFSET + i * NJ + k, 7, 0, idx);
				rtTmpAccess(C_OFFSET + k * NL + j, 8, 4, idx);
				rtTmpAccess(D_OFFSET + i * NL + j, 9, 3, idx);
				rtTmpAccess(D_OFFSET + i * NL + j, 10, 3, idx);
			}
        }
    }
}


int main(int argc, char* argv[]) {
	
    if (argc != 5) {
        cout << "This benchmark needs 4 loop bounds" << endl;
        return 0;
    }
    for (int i = 1; i < argc; i++) {
        if (!isdigit(argv[i][0])) {
            cout << "arguments must be integer" << endl;
            return 0;
        }
    }
    
    NI = stoi(argv[1]);
    NJ = stoi(argv[2]);
    NK = stoi(argv[3]);
    NL = stoi(argv[4]);
    
	double* tmp = (double*)malloc( NI * NJ * sizeof(double));
	double* A = (double*)malloc( NI * NK * sizeof(double));
	double* B = (double*)malloc( NK * NJ * sizeof(double));
	double* C = (double*)malloc( NJ * NL * sizeof(double));
	double* D = (double*)malloc( NI * NL * sizeof(double));
	double alpha = 0.1;
	double beta = 0.5;
	
	mm2_trace(tmp, A, B, C, D, alpha, beta);

    string name(argv[0]);
    size_t found = name.find_last_of("/\\") + 1;
    string conf = name.substr(found, name.size()-found) + "_" + to_string(NI) + "_" + to_string(NJ) + "_" + to_string(NK) + "_" + to_string(NL);
    
    dumpRIHistogram(conf);
    predictionWithBmc(conf);
    predictionWithRefBmc(conf);   
	return 0;
}

