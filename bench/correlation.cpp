#include "./utility/rt.h"
#include<math.h>

int M;
int N;

#define FLOAT_N 3214212.01
#define EPS 0.005

#define DATA_OFFSET 0
#define MEAN_OFFSET M * N
#define STDDEV_OFFSET M * N + M
#define SYMMAT_OFFSET M * N + M + M

void correlation_trace(double* data, double* mean, double* stddev, double* symmat) {

	int i, j, j1, j2;
    vector<int> idx;

    // Determine mean of column vectors of input data matrix
    for (j = 0; j < M; j++)
    {
        idx.clear(); idx.push_back(j);
        mean[j] = 0.0;
		rtTmpAccess(MEAN_OFFSET + j, 0, 0, idx);

        for (i = 0; i < N; i++)
        {
            idx.clear(); idx.push_back(j); idx.push_back(i);
            mean[j] += data[i * M + j];
            rtTmpAccess(DATA_OFFSET + i * M + j, 1, 1, idx);
        	rtTmpAccess(MEAN_OFFSET + j, 2, 0, idx);
			rtTmpAccess(MEAN_OFFSET + j, 3, 0, idx);
		}

        idx.clear(); idx.push_back(j);
        mean[j] /= (double)FLOAT_N;
    	rtTmpAccess(MEAN_OFFSET + j, 4, 0, idx);
		rtTmpAccess(MEAN_OFFSET + j, 5, 0, idx);
	}

    // Determine standard deviations of column vectors of data matrix.
    for (j = 0; j < M; j++)
    {
        idx.clear(); idx.push_back(j);
        stddev[j] = 0.0;
		rtTmpAccess(STDDEV_OFFSET + j, 6, 2, idx);

        for (i = 0; i < N; i++)
        {
            idx.clear(); idx.push_back(j); idx.push_back(i);
            stddev[j] += (data[i * M + j] - mean[j]) * (data[i * M + j] - mean[j]);
        	
			rtTmpAccess(DATA_OFFSET + i * M + j, 7, 1, idx);
			rtTmpAccess(MEAN_OFFSET + j, 8, 0, idx);
			rtTmpAccess(DATA_OFFSET + i * M + j, 9, 1, idx);
			rtTmpAccess(MEAN_OFFSET + j, 10, 0, idx);
            rtTmpAccess(STDDEV_OFFSET + j, 11, 2, idx);
			rtTmpAccess(STDDEV_OFFSET + j, 12, 2, idx);
		}

        idx.clear(); idx.push_back(j);
        stddev[j] /= FLOAT_N;
//        stddev[j] = sqrt_of_array_cell(stddev, j);
		stddev[j] = sqrt(stddev[j]);

        stddev[j] = stddev[j] <= EPS ? 1.0 : stddev[j];
    	rtTmpAccess(STDDEV_OFFSET + j, 13, 2, idx);
		rtTmpAccess(STDDEV_OFFSET + j, 14, 2, idx);
		rtTmpAccess(STDDEV_OFFSET + j, 15, 2, idx);
		rtTmpAccess(STDDEV_OFFSET + j, 16, 2, idx);
        
        if (stddev[j] > EPS) {
            rtTmpAccess(STDDEV_OFFSET + j, 17, 2, idx);
        }
		rtTmpAccess(STDDEV_OFFSET + j, 18, 2, idx);
	}

	// Center and reduce the column vectors.
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            idx.clear(); idx.push_back(i); idx.push_back(j);
            data[i * M + j] -= mean[j];
            data[i * M + j] /= (sqrt(FLOAT_N)*stddev[j]) ;
            
            rtTmpAccess(MEAN_OFFSET + j, 19, 0, idx);
			rtTmpAccess(DATA_OFFSET + i * M + j, 20, 1, idx);
			rtTmpAccess(DATA_OFFSET + i * M + j, 21, 1, idx);
            rtTmpAccess(STDDEV_OFFSET + j, 22, 2, idx);
			rtTmpAccess(DATA_OFFSET + i * M + j, 23, 1, idx);
			rtTmpAccess(DATA_OFFSET + i * M + j, 24, 1, idx);
        }
    }

    // Calculate the m * m correlation matrix.
    for (j1 = 0; j1 < M-1; j1++)
    {
        idx.clear(); idx.push_back(j1);
        symmat[j1 * M + j1] = 1.0;
		rtTmpAccess(SYMMAT_OFFSET + j1 * M + j1, 25, 3, idx);

        for (j2 = j1+1; j2 < M; j2++)
        {
            idx.clear(); idx.push_back(j1); idx.push_back(j2);
            symmat[j1 * M + j2] = 0.0;
			rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 26, 3, idx);

            for (i = 0; i < N; i++)
            {
                idx.clear(); idx.push_back(j1); idx.push_back(j2); idx.push_back(i);
                symmat[j1 * M + j2] += (data[i * M + j1] * data[i * M + j2]);
                rtTmpAccess(DATA_OFFSET + i * M + j1, 27, 1, idx);
                rtTmpAccess(DATA_OFFSET + i * M + j2, 28, 1, idx);
                rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 29, 3, idx);
				rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 30, 3, idx);
            }

            idx.clear(); idx.push_back(j1); idx.push_back(j2);
            symmat[j2 * M + j1] = symmat[j1 * M + j2];
			rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 31, 3, idx);
			rtTmpAccess(SYMMAT_OFFSET + j2 * M + j1, 32, 3, idx);
        }
    }
    
    idx.clear();
    symmat[(M-1) * M + M-1] = 1.0;	
	rtTmpAccess(SYMMAT_OFFSET + (M-1) * M + M-1, 33, 3, idx);

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
    double* stddev = (double*)malloc( M*sizeof(double));
    double* symmat = (double*)malloc( (M*M)*sizeof(double));

    for (int i = 0; i < M; i++) {
        data[i] = i % 256;
    }

    correlation_trace(data, mean, stddev, symmat);

    dumpRIHistogram();
    
    return 0;
}
