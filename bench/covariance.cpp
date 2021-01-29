#include "./utility/rt.h"

int M;
int N;

#define DATA_OFFSET 0
#define MEAN_OFFSET M * N
#define SYMMAT_OFFSET M * N + M

#define FLOAT_N 3214212.01
#define EPS 0.005

void covariance_trace(double* data, double* mean, double* symmat, unsigned int m, unsigned int n) {

	int j, i, j1, j2;
    vector<int> idx;
    
    /* Determine mean of column vectors of input data matrix */
    for (j = 0; j < M; j++)
    {
        idx.clear(); idx.push_back(j);
        mean[j] = 0.0;
        rtTmpAccess(MEAN_OFFSET + j, 0, 0, idx);
        
        for (i = 0; i < N; i++)
        {
            idx.clear(); idx.push_back(j); idx.push_back(i);
            mean[j] = mean[j] + data[i * M + j];
            rtTmpAccess(DATA_OFFSET + i * M + j, 1, 1, idx);
            rtTmpAccess(MEAN_OFFSET + j, 2, 0, idx);
            rtTmpAccess(MEAN_OFFSET + j, 3, 0, idx);
        }
        idx.clear(); idx.push_back(j);
        mean[j] = mean[j] / FLOAT_N;
        rtTmpAccess(MEAN_OFFSET + j, 4, 0, idx);
        rtTmpAccess(MEAN_OFFSET + j, 5, 0, idx);
    }
    
    /* Center the column vectors. */
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            data[i * M + j] = data[i * M + j] - mean[j];
            rtTmpAccess(MEAN_OFFSET + j, 6, 0, idx);
            rtTmpAccess(DATA_OFFSET + i * M + j, 7, 1, idx);
            rtTmpAccess(DATA_OFFSET + i * M + j, 8, 1, idx);
        }
    }
    
    /* Calculate the m * m covariance matrix. */
    for (j1 = 0; j1 < M; j1++)
    {
        for (j2 = j1; j2 < M; j2++)
        {
            idx.clear(); idx.push_back(j1); idx.push_back(j2);
            symmat[j1 * M + j2] = 0.0;
            rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 9, 2, idx);
            
            for (i = 0; i < N; i++)
            {
                idx.clear(); idx.push_back(j1); idx.push_back(j2); idx.push_back(i);
                symmat[j1 * M + j2] = symmat[j1 * M + j2] + data[i * M + j1] * data[i * M + j2];
                rtTmpAccess(DATA_OFFSET + i * M + j1, 10, 1, idx);
                rtTmpAccess(DATA_OFFSET + i * M + j2, 11, 1, idx);
                rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 12, 2, idx);
                rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 13, 2, idx);
            }
            
            idx.clear(); idx.push_back(j1); idx.push_back(j2);
            symmat[j2 * M + j1] = symmat[j1 * M + j2];
            rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 14, 2, idx);
            rtTmpAccess(SYMMAT_OFFSET + j2 * M + j1, 15, 2, idx);
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
    
    M = stoi(argv[1]);
    N = stoi(argv[2]);
    
    double* data = (double*)malloc( (M*N)*sizeof(double));
    double* mean = (double*)malloc( M*sizeof(double));
    double* symmat = (double*)malloc( (M*M)*sizeof(double));

    for (int i = 0; i < M*N; i++) {
        data[i] = i % 256;
    }

    covariance_trace(data, mean, symmat, M, N);
    
    dumpRIHistogram();
    
    return 0;
}
