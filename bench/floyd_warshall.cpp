#include "./utility/rt.h"

int N;

#define PATH_OFFSET 0

void floyd_warshall_trace(double* path) {

	int k, i, j;
    vector<int> idx;

	for (k = 0; k < N; k++) 
		for(i = 0; i < N; i++)
			for (j = 0; j < N; j++) {
                idx.clear(); idx.push_back(k); idx.push_back(i); idx.push_back(j);
				//path[i * N + j] = path[i * N + j] < path[i * N + k] + path[k * N + j] ? path[i * N + j] : path[i * N + k] + path[k * N + j];
				rtTmpAccess(PATH_OFFSET + i * N + j, 0, 0, idx);
				rtTmpAccess(PATH_OFFSET + i * N + k, 1, 0, idx);
				rtTmpAccess(PATH_OFFSET + k * N + j, 2, 0, idx);
				if (path[i * N + j] < path[i * N + k] + path[k * N + j]) {
					path[i * N + j] = path[i * N + j];
					rtTmpAccess(PATH_OFFSET + i * N + j, 3, 0, idx);
				} else {
					path[i * N + j] = path[i * N + k] + path[k * N + j];
					rtTmpAccess(PATH_OFFSET + i * N + k, 4, 0, idx);
					rtTmpAccess(PATH_OFFSET + k * N + j, 5, 0, idx);
					rtTmpAccess(PATH_OFFSET + i * N + j, 6, 0, idx);
				}
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
    
	double* path = (double *)malloc(N * N * sizeof(double));

	floyd_warshall_trace(path);

    dumpRIHistogram();
    
	return 0;
}
