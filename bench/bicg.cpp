#include "./utility/rt.h"

int NX;
int NY;

#define A_OFFSET 0
#define S_OFFSET NX * NY
#define Q_OFFSET NX * NY + NY
#define R_OFFSET NX * NY + NY + NX
#define P_OFFSET NX * NY + NY + NX + NX


void bicg_cpu_trace(double* A, double* r, double* s, double* p, double* q, unsigned int nx, int ny) {
 
	int i,j;
    vector<int> idx;

    for (i = 0; i < NY; i++)
    {
        idx.clear(); idx.push_back(i);
        s[i] = 0.0;
    	rtTmpAccess(S_OFFSET + i, 0, 0, idx);
	}

    for (i = 0; i < NX; i++)
    {
        idx.clear(); idx.push_back(i);
        q[i] = 0.0;
		rtTmpAccess(Q_OFFSET + i, 1, 1, idx);
        for (j = 0; j < NY; j++)
        {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            s[j] = s[j] + r[i] * A[i * NY + j];
            q[i] = q[i] + A[i * NY + j] * p[j];
        	rtTmpAccess(S_OFFSET + j, 2, 0, idx);
			rtTmpAccess(R_OFFSET + i, 3, 2, idx);
			rtTmpAccess(A_OFFSET + i * NY + j, 4, 3, idx);
			rtTmpAccess(S_OFFSET + j, 5, 0, idx);
			rtTmpAccess(Q_OFFSET + i, 6, 1, idx);
			rtTmpAccess(A_OFFSET + i * NY + j, 7, 3, idx);
			rtTmpAccess(P_OFFSET + j, 8, 4, idx);
			rtTmpAccess(Q_OFFSET + i, 9, 1, idx);
		}
    }

    return;
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
    
    NX = stoi(argv[1]);
    NY = stoi(argv[2]);
    
    double* A = (double*)malloc( (NX*NY)*sizeof(double));
    double* r = (double*)malloc( NX*sizeof(double));
    double* s = (double*)malloc( NY*sizeof(double));
    double* q = (double*)malloc( NX*sizeof(double));
    double* p = (double*)malloc( NY*sizeof(double));

    for (int i = 0; i < NX; i++) {
        r[i] = i % 256;
    }

    for (int i = 0; i< NY; i++) {
        p[i] = i % 256;
    }

    for (int i = 0; i< NY*NX; i++) {
        A[i] = i % 128;
    }

    bicg_cpu_trace(A, r, s, p ,q, NX, NY);

    dumpRIHistogram();
    
    return 0;
}
