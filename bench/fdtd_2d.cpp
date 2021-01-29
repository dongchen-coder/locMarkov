#include "./utility/rt.h"

int _PB_TMAX;
int _PB_NY;
int _PB_NX;

#define FICT_OFFSET 0
#define EY_OFFSET _PB_TMAX
#define EX_OFFSET _PB_TMAX + _PB_NX * _PB_NY
#define HZ_OFFSET _PB_TMAX + _PB_NX * _PB_NY + _PB_NX * _PB_NY


void fdtd_2d_trace(double* _fict_, double* ey, double* ex, double* hz) {

	int t, i, j;
    vector<int> idx;
    
	for(t = 0; t < _PB_TMAX; t++) {
		for (j = 0; j < _PB_NY; j++) {
            idx.clear(); idx.push_back(i); idx.push_back(j);
			ey[0 * _PB_NY + j] = _fict_[t];
			rtTmpAccess(FICT_OFFSET + t, 0, 0, idx);
			rtTmpAccess(EY_OFFSET + 0 * _PB_NY + j, 1, 1, idx);
		}

		for (i = 1; i < _PB_NX; i++)
			for (j = 0; j < _PB_NY; j++) {
                idx.clear(); idx.push_back(i); idx.push_back(j);
				ey[i * _PB_NY + j] = ey[i * _PB_NY + j] - 0.5 * (hz[i * _PB_NY + j] - hz[(i-1) * _PB_NY + j]);
				rtTmpAccess(EY_OFFSET + i * _PB_NY + j, 2, 1, idx);
				rtTmpAccess(HZ_OFFSET + i * _PB_NY + j, 3, 2, idx);
				rtTmpAccess(HZ_OFFSET + (i-1) * _PB_NY + j, 4, 2, idx);
				rtTmpAccess(EY_OFFSET + i * _PB_NY + j, 5, 1, idx);
			}

		for (i = 0; i < _PB_NX; i++)
			for (j = 1; j < _PB_NY; j++) {
                idx.clear(); idx.push_back(i); idx.push_back(j);
				ex[i * _PB_NY + j] = ex[i * _PB_NY + j] - 0.5 * (hz[i * _PB_NY + j]-hz[i * _PB_NY + j - 1]);
				rtTmpAccess(EX_OFFSET + i * _PB_NY + j, 6, 3, idx);
				rtTmpAccess(HZ_OFFSET + i * _PB_NY + j, 7, 2, idx);
				rtTmpAccess(HZ_OFFSET + i * _PB_NY + j - 1, 8, 2, idx);
				rtTmpAccess(EX_OFFSET + i * _PB_NY + j, 9, 3, idx);
			}

		for (i = 0; i < _PB_NX - 1; i++)
			for (j = 0; j < _PB_NY - 1; j++) {
                idx.clear(); idx.push_back(i); idx.push_back(j);
				hz[i * _PB_NY + j] = hz[i * _PB_NY + j] - 0.7 * (ex[i * _PB_NY + j + 1] - ex[i * _PB_NY + j] + ey[(i+1) * _PB_NY + j] - ey[i * _PB_NY + j]);
    			rtTmpAccess(HZ_OFFSET + i * _PB_NY + j, 10, 2, idx);
				rtTmpAccess(EX_OFFSET + i * _PB_NY + j + 1, 11, 3, idx);
				rtTmpAccess(EX_OFFSET + i * _PB_NY + j, 12, 3, idx);
				rtTmpAccess(EY_OFFSET + (i+1) * _PB_NY + j, 13, 1, idx);
				rtTmpAccess(EY_OFFSET + i * _PB_NY + j, 14, 1, idx);
				rtTmpAccess(HZ_OFFSET + i * _PB_NY + j, 15, 2, idx);
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
    
    _PB_TMAX = stoi(argv[1]);
    _PB_NX = stoi(argv[2]);
    _PB_NY = stoi(argv[3]);
    
	double* _fict_ = (double *)malloc(_PB_TMAX * sizeof(double));
	double* ey = (double *)malloc(_PB_NX * _PB_NY * sizeof(double));
	double* ex = (double *)malloc(_PB_NX * _PB_NY * sizeof(double));
	double* hz = (double *)malloc(_PB_NX * _PB_NY * sizeof(double));	

	fdtd_2d_trace(_fict_, ey, ex, hz);

    dumpRIHistogram();
    
	return 0;
}
