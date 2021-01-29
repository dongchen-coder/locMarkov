#include "./utility/rt.h"

int NI;
int NJ;
int NL;
int NK;
int NM;

#define A_OFFSET 0
#define B_OFFSET NK * NI
#define C_OFFSET NK * NI + NK * NJ
#define D_OFFSET NK * NI + NK * NJ + NJ * NM
#define E_OFFSET NK * NI + NK * NJ + NJ * NM + NM * NL
#define F_OFFSET NK * NI + NK * NJ + NJ * NM + NM * NL + NI * NJ
#define G_OFFSET NK * NI + NK * NJ + NJ * NM + NM * NL + NI * NJ + NJ * NL

void mm3_cpu_trace(int *A, int *B, int *C, int *D, int *E, int *F, int *G) {
    
	int i, j, k;
    vector<int> idx;

    /* E := A*B */
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NJ; j++)
        {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            E[i * NJ + j] = 0;
			
			rtTmpAccess(E_OFFSET + i * NJ + j, 0, 0, idx);
	
            for (k = 0; k < NK; ++k)
            {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
                E[i * NJ + j] += A[i * NK + k] * B[k * NJ + j];
				rtTmpAccess(A_OFFSET + i * NK + k, 1, 1, idx);
				rtTmpAccess(B_OFFSET + k * NJ + j, 2, 2, idx);
				rtTmpAccess(E_OFFSET + i * NJ + j, 3, 0, idx);
				rtTmpAccess(E_OFFSET + i * NJ + j, 4, 0, idx);
            }
        }
    }

    /* F := C*D */
    for (i = 0; i < NJ; i++)
    {
        for (j = 0; j < NL; j++)
        {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            F[i * NL + j] = 0;

			rtTmpAccess(F_OFFSET + i * NL + j, 5, 3, idx);
			
            for (k = 0; k < NM; ++k)
            {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
                F[i * NL + j] += C[i * NM + k] * D[k * NL + j];
				rtTmpAccess(C_OFFSET + i * NM + k, 6, 4, idx);
				rtTmpAccess(D_OFFSET + k * NL + j, 7, 5, idx);
				rtTmpAccess(F_OFFSET + i * NL + j, 8, 3, idx);
				rtTmpAccess(F_OFFSET + i * NL + j, 9, 3, idx);
			}
        }
    }

    /* G := E*F */
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NL; j++)
        {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            G[i * NL + j] = 0;

			rtTmpAccess(G_OFFSET + i * NL + j, 10, 6, idx);

            for (k = 0; k < NJ; ++k)
            {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
                G[i * NL + j] += E[i * NJ + k] * F[k * NL + j];
				rtTmpAccess(E_OFFSET + i * NJ + k, 11, 0, idx);
				rtTmpAccess(F_OFFSET + k * NL + j, 12, 3, idx);
				rtTmpAccess(G_OFFSET + i * NL + j, 13, 6, idx);
				rtTmpAccess(G_OFFSET + i * NL + j, 14, 6, idx);
			}
        }
    }

    return;
}



int main(int argc, char* argv[]) {
    
    if (argc != 6) {
        cout << "This benchmark needs 5 loop bounds" << endl;
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
    NM = stoi(argv[5]);
    
    int* a = (int*)malloc( NK * NI * sizeof(int) );
    int* b = (int*)malloc( NK * NJ * sizeof(int) );
    int* c = (int*)malloc( NJ * NM * sizeof(int) );
    int* d = (int*)malloc( NM * NL * sizeof(int) );
    int* e = (int*)malloc( NI * NJ * sizeof(int) );
    int* f = (int*)malloc( NJ * NL * sizeof(int) );
    int* g = (int*)malloc( NI * NL * sizeof(int) );

    for (int i = 0; i < NK * NI; i++) {
            a[i] = i % 256;
    }

	for (int i = 0; i < NK * NJ; i++) {
            b[i] = i % 256;
    }

	for (int i = 0; i < NJ * NM; i++) {
            c[i] = i % 256;
    }

	for (int i = 0; i < NM * NL; i++) {
            d[i] = i % 256;
    }

    mm3_cpu_trace(a, b, c, d, e, f, g);
    
    dumpRIHistogram();
    
    return 0;
}
