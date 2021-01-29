#include "./utility/rt.h"

int NX;
int NY;


#define A_OFFSET 0
#define X_OFFSET NX * NY
#define Y_OFFSET NX * NY + NY
#define TMP_OFFSET NX * NY + NY + NY


void atax_cpu_trace(double* A, double* x, double* y, double* tmp, unsigned int dim_size_x, unsigned int dim_size_y) {
    
	int i,j;
    vector<int> idx;

    for (i= 0; i < NY; i++)
    {
        idx.clear(); idx.push_back(i);
        y[i] = 0;
    	rtTmpAccess(Y_OFFSET + i, 0, 0, idx);
	}

    for (i = 0; i < NX; i++)
    {
        idx.clear(); idx.push_back(i);
        tmp[i] = 0;
        rtTmpAccess(TMP_OFFSET + i, 1, 1, idx);
            
        for (j = 0; j < NY; j++) {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            tmp[i] = tmp[i] + A[i * NY + j] * x[j];
            rtTmpAccess(TMP_OFFSET + i, 2, 1, idx);
            rtTmpAccess(A_OFFSET + i * NY + j, 3, 2, idx);
            rtTmpAccess(X_OFFSET + j, 4, 3, idx);
            rtTmpAccess(TMP_OFFSET + i, 5, 1, idx);
        }

        for (j = 0; j < NY; j++) {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            y[j] = y[j] + A[i * NY + j] * tmp[i];
            rtTmpAccess(Y_OFFSET + j, 6, 0, idx);
            rtTmpAccess(A_OFFSET + i * NY + j, 7, 2, idx);
            rtTmpAccess(TMP_OFFSET + i, 8, 1, idx);
            rtTmpAccess(Y_OFFSET + j, 9, 0, idx);
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
    
    double* A = (double*)malloc( (NX * NY) * sizeof(double) );
    double* x = (double*)malloc( NX * sizeof(double) );
    double* y = (double*)malloc( NY * sizeof(double) );
    double* tmp = (double*)malloc( NX * sizeof(double) );

    for (int i = 0; i < (NX * NY); i++) {
        A[i] = i % 256;
    }

    for (int i = 0; i < NX; i++) {
        x[i] = i % 256;
    }

    atax_cpu_trace(A, x, y, tmp, NX, NY);

    dumpRIHistogram();
    
    return 0;
}
