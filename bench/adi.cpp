#include "./utility/rt.h"

int N;
int TSTEPS;

#define P_OFFSET 0
#define Q_OFFSET 0 + N * N
#define V_OFFSET 0 + N * N + N * N 
#define U_OFFSET 0 + N * N + N * N + N * N

void adi_trace(double* p, double* q, double* v, double* u) {

	int t, i, j;
	double DX, DY, DT, B1, B2, mul1, mul2, a, b, c, d, e, f;
    vector<int> idx;

	DX = 1.0 / N;
	DY = 1.0 / N;
	DT = 1.0 / TSTEPS;
	B1 = 2.0;
	B2 = 1.0;
	mul1 = B1 * DT / (DX * DX);
	mul2 = B2 * DT / (DY * DY);

	a = -mul1 / 2.0;
	b = 1.0 + mul1;
	c = a;
	d = -mul2 / 2.0;
	e = 1.0 + mul2;
	f = d;

	for (t=1; t <= TSTEPS; t++) {
	//Column Sweep
		for (i=1; i< N-1; i++) {
            idx.clear(); idx.push_back(t); idx.push_back(i);
			v[0 * N + i] = 1.0;
			p[i * N + 0] = 0.0;
			q[i * N + 0] = v[0 * N + i];
			rtTmpAccess(V_OFFSET + 0 * N + i, 0, 0, idx);
			rtTmpAccess(P_OFFSET + i * N + 0, 1, 1, idx);
			rtTmpAccess(V_OFFSET + 0 * N + i, 2, 0, idx);
			rtTmpAccess(Q_OFFSET + i * N + 0, 3, 2, idx);
			for (j=1; j< N-1; j++) {
                idx.clear(); idx.push_back(t); idx.push_back(i); idx.push_back(j);
				p[i * N + j] = -c / (a*p[i * N + j-1]+b);
				q[i * N + j] = (-d*u[j * N + i-1]+(1.0 + 2.0 * d)*u[j * N + i] - f*u[j * N + i+1]-a*q[i * N + j-1])/(a*p[i * N + j-1]+b);
				rtTmpAccess(P_OFFSET + i * N + j-1, 4, 1, idx);
				rtTmpAccess(P_OFFSET + i * N + j, 5, 1, idx);
				rtTmpAccess(U_OFFSET + j * N + i-1, 6, 3, idx);
				rtTmpAccess(U_OFFSET + j * N + i, 7, 3, idx);
				rtTmpAccess(U_OFFSET + j * N + i+1, 8, 3, idx);
				rtTmpAccess(Q_OFFSET + i * N + j-1, 9, 2, idx);
				rtTmpAccess(P_OFFSET + i * N + j-1, 10, 1, idx);
				rtTmpAccess(Q_OFFSET + i * N + j, 11, 2, idx);
			}
            idx.clear(); idx.push_back(t); idx.push_back(i);
			v[(N - 1) * N + i] = 1.0;
			rtTmpAccess(V_OFFSET + (N - 1) * N + i, 12, 0, idx);
			for (j= N-2; j>=1; j--) {
                idx.clear(); idx.push_back(t); idx.push_back(i); idx.push_back(j);
				v[j * N + i] = p[i * N + j] * v[(j+1) * N + i] + q[i * N + j];
				rtTmpAccess(P_OFFSET + i * N + j, 13, 1, idx);
				rtTmpAccess(V_OFFSET + (j+1) * N + i, 14, 0, idx);
				rtTmpAccess(Q_OFFSET + i * N + j, 15, 2, idx);
				rtTmpAccess(V_OFFSET + j * N + i, 16, 0, idx);
			}
		}
		//Row Sweep
		for (i=1; i < N - 1; i++) {
            idx.clear(); idx.push_back(t); idx.push_back(i);
			u[i * N + 0] = 1.0;
			p[i * N + 0] = 0.0;
			q[i * N + 0] = u[i * N + 0];
			rtTmpAccess(U_OFFSET + i * N + 0, 17, 3, idx);
			rtTmpAccess(P_OFFSET + i * N + 0, 18, 1, idx);
			rtTmpAccess(U_OFFSET + i * N + 0, 19, 3, idx);
			rtTmpAccess(Q_OFFSET + i * N + 0, 20, 2, idx);
			for (j=1; j< N - 1; j++) {
                idx.clear(); idx.push_back(t); idx.push_back(i); idx.push_back(j);
				p[i * N + j] = -f / (d*p[i * N + j-1]+e);
				q[i * N + j] = (-a*v[(i-1) * N + j]+(1.0 + 2.0 * a)*v[i * N + j] - c*v[(i+1) * N + j]-d*q[i * N + j-1])/(d*p[i * N + j-1]+e);
				rtTmpAccess(P_OFFSET + i * N + j-1, 21, 1, idx);
				rtTmpAccess(P_OFFSET + i * N + j, 22, 1, idx);
				rtTmpAccess(V_OFFSET + (i-1) * N + j, 23, 0, idx);
				rtTmpAccess(V_OFFSET + i * N + j, 24, 0, idx);
				rtTmpAccess(V_OFFSET + (i+1) * N + j, 25, 0, idx);
				rtTmpAccess(Q_OFFSET + i * N + j-1, 26, 2, idx);
				rtTmpAccess(P_OFFSET + i * N + j-1, 27, 1, idx);
				rtTmpAccess(Q_OFFSET + i * N + j, 28, 2, idx);
			}
            idx.clear(); idx.push_back(t); idx.push_back(i);
			u[i * N + N - 1 ] = 1.0;
			rtTmpAccess(U_OFFSET + i * N + N - 1, 29, 3, idx);
			for (j= N - 2; j>=1; j--) {
                idx.clear(); idx.push_back(t); idx.push_back(i); idx.push_back(j);
				u[i * N + j] = p[i * N + j] * u[i * N + j+1] + q[i * N + j];
				rtTmpAccess(P_OFFSET + i * N + j, 30, 1, idx);
				rtTmpAccess(U_OFFSET + i * N + j+1, 31, 3, idx);
				rtTmpAccess(Q_OFFSET + i * N + j, 32, 2, idx);
				rtTmpAccess(U_OFFSET + i * N + j, 33, 3, idx);
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
    
	double * p = (double *) malloc(N * N * sizeof(double));
	double * q = (double *) malloc(N * N * sizeof(double));
	double * v = (double *) malloc(N * N * sizeof(double));
	double * u = (double *) malloc(N * N * sizeof(double));

	adi_trace(p, q, v, u);

    dumpRIHistogram();
    
	return 0;
}


