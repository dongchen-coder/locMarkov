#include "./utility/rt.h"

#include <math.h>

int N;

#define match(b1, b2) (((b1)+(b2)) == 3 ? 1 : 0)
#define max_score(s1, s2) ((s1 >= s2) ? s1 : s2)


#define TABLE_OFFSET 0
#define SEQ_OFFSET N * N

void nussinov_trace(double* table, double* seq) {

	int i, j, k;
    vector<int> idx;

	for (i = N - 1; i >= 0; i--) {
		for (j=i+1; j < N; j++) {
            idx.clear(); idx.push_back(i); idx.push_back(j);
			if (j-1>=0) {
				table[i * N + j] = max_score(table[i * N + j], table[i * N + j-1]);
				rtTmpAccess(TABLE_OFFSET + i * N + j, 0, 0, idx);
				rtTmpAccess(TABLE_OFFSET + i * N + j-1, 1, 0, idx);
				//if (table[i * N + j] > table[i * N + j-1])
				//	rtTmpAccess(TABLE_OFFSET + i * N + j, 2, 0);
				//else
				//	rtTmpAccess(TABLE_OFFSET + i * N + j-1, 3, 0);
				rtTmpAccess(TABLE_OFFSET + i * N + j, 2, 0, idx);
			}
			if (i+1 < N) {
				table[i * N +j] = max_score(table[i * N + j], table[(i+1) * N + j]);
				rtTmpAccess(TABLE_OFFSET + i * N + j, 3, 0, idx);
				rtTmpAccess(TABLE_OFFSET + (i+1) * N + j, 4, 0, idx);
				//if (table[i * N + j] >= table[(i+1) * N + j])
				//	rtTmpAccess(TABLE_OFFSET + i * N + j, 7, 0);
				//else
				//	rtTmpAccess(TABLE_OFFSET + i * N + j, 8, 0);
				rtTmpAccess(TABLE_OFFSET + i * N + j, 5, 0, idx);
			}
			if (j-1>=0 && i+1 < N) {
				/* don't allow adjacent elements to bond */
				if (i<j-1) {
					table[i * N + j] = max_score(table[i * N + j], table[(i+1) * N + j-1]+match(seq[i], seq[j]));
					rtTmpAccess(TABLE_OFFSET + i * N + j, 6, 0, idx);
					rtTmpAccess(TABLE_OFFSET + (i+1) * N + j-1, 7, 0, idx);
					rtTmpAccess(SEQ_OFFSET + i, 8, 1, idx);
					rtTmpAccess(SEQ_OFFSET + j, 9, 1, idx);
                    
					//if (table[i * N + j] >= table[(i+1) * N + j-1]+match(seq[i], seq[j])) {
					//	rtTmpAccess(TABLE_OFFSET + i * N + j, 14, 0);
					//} else {
					//	rtTmpAccess(TABLE_OFFSET + (i+1) * N + j-1, 15, 0);
                    //	rtTmpAccess(SEQ_OFFSET + i, 16, 1);
                    //	rtTmpAccess(SEQ_OFFSET + j, 17, 1);
					//}
					rtTmpAccess(TABLE_OFFSET + i * N + j, 10, 0, idx);
				} else {
					table[i * N + j] = max_score(table[i * N + j], table[(i+1) * N + j-1]);
					rtTmpAccess(TABLE_OFFSET + i * N + j, 11, 0, idx);
					rtTmpAccess(TABLE_OFFSET + (i+1) * N + j-1, 12, 0, idx);
					//if (table[i * N + j] >= table[(i+1) * N + j-1]) {
					//	rtTmpAccess(TABLE_OFFSET + i * N + j, 21, 0);
					//} else {
					//	rtTmpAccess(TABLE_OFFSET + i * N + j, 22, 0);
					//}
					rtTmpAccess(TABLE_OFFSET + i * N + j, 13, 0, idx);
				}
			}
			for (k=i+1; k<j; k++) {
				table[i * N + j] = max_score(table[i * N + j], table[i * N + k] + table[(k+1) * N + j]);
				rtTmpAccess(TABLE_OFFSET + i * N + j, 14, 0, idx);
				rtTmpAccess(TABLE_OFFSET + i * N + k, 15, 0, idx);
				rtTmpAccess(TABLE_OFFSET + (k+1) * N + j, 16, 0, idx);
				//if (table[i * N + j] >= table[i * N + k] + table[(k+1) * N + j]) {
				//	rtTmpAccess(TABLE_OFFSET + i * N + j, 27, 0);
				//} else {
				//	rtTmpAccess(TABLE_OFFSET + i * N + k, 28, 0);
                //	rtTmpAccess(TABLE_OFFSET + (k+1) * N + j, 29, 0);
				//}
				rtTmpAccess(TABLE_OFFSET + i * N + j, 17, 0, idx);
			}
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
    
	double* table = (double *)malloc(N * N * sizeof(double));
	double* seq = (double *)malloc(N * sizeof(double));

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			table[i * N + j] = (i * N + j) % 128;
		}
		seq[i] = i % 64;
	}

	nussinov_trace(table, seq);

    dumpRIHistogram();
    
	return 0;
}
