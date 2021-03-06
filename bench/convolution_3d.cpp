#include "./utility/mc_kth.h"

int NI;
int NJ;
int NK;

#define A_OFFSET 0
#define B_OFFSET NI * NJ * NK

void conv3D_trace(double* A, double* B)
{
    int i, j, k;
    double c11, c12, c13, c21, c22, c23, c31, c32, c33;
    vector<int> idx;
    
    c11 = +2;  c21 = +5;  c31 = -8;
    c12 = -3;  c22 = +6;  c32 = -9;
    c13 = +4;  c23 = +7;  c33 = +10;

    for (int i = 1; i < NI - 1; ++i) // 0
    {
        for (int j = 1; j < NJ - 1; ++j) // 1
        {
            for (int k = 1; k < NK -1; ++k) // 2
            {
                idx.clear(); idx.push_back(i); idx.push_back(j); idx.push_back(k);
                /*
                B[i][j][k] = c11 * A[(i - 1)][(j - 1)][(k - 1)]  +  c13 * A[(i + 1)][(j - 1)][(k - 1)]
                +   c21 * A[(i - 1)][(j - 1)][(k - 1)]  +  c23 * A[(i + 1)][(j - 1)][(k - 1)]
                +   c31 * A[(i - 1)][(j - 1)][(k - 1)]  +  c33 * A[(i + 1)][(j - 1)][(k - 1)]
                +   c12 * A[(i + 0)][(j - 1)][(k + 0)]  +  c22 * A[(i + 0)][(j + 0)][(k + 0)]
                +   c32 * A[(i + 0)][(j + 1)][(k + 0)]  +  c11 * A[(i - 1)][(j - 1)][(k + 1)]
                +   c13 * A[(i + 1)][(j - 1)][(k + 1)]  +  c21 * A[(i - 1)][(j + 0)][(k + 1)]
                +   c23 * A[(i + 1)][(j + 0)][(k + 1)]  +  c31 * A[(i - 1)][(j + 1)][(k + 1)]
                +   c33 * A[(i + 1)][(j + 1)][(k + 1)];
                 */
                B[i * NJ * NK + j * NK + k] = c11 * A[(i - 1) * NJ * NK + (j - 1) * NK + (k - 1)]
                                            +   c13 * A[(i + 1) * NJ * NK + (j - 1) * NK + (k - 1)]
                                            +   c21 * A[(i - 1) * NJ * NK + (j - 1) * NK + (k - 1)]
                                            +   c23 * A[(i + 1) * NJ * NK + (j - 1) * NK + (k - 1)]
                                            +   c31 * A[(i - 1) * NJ * NK + (j - 1) * NK + (k - 1)]
                                            +   c33 * A[(i + 1) * NJ * NK + (j - 1) * NK + (k - 1)]
                                            +   c12 * A[(i + 0) * NJ * NK + (j - 1) * NK + (k + 0)]
                                            +   c22 * A[(i + 0) * NJ * NK + (j + 0) * NK + (k + 0)]
                                            +   c32 * A[(i + 0) * NJ * NK + (j + 1) * NK + (k + 0)]
                                            +   c11 * A[(i - 1) * NJ * NK + (j - 1) * NK + (k + 1)]
                                            +   c13 * A[(i + 1) * NJ * NK + (j - 1) * NK + (k + 1)]
                                            +   c21 * A[(i - 1) * NJ * NK + (j + 0) * NK + (k + 1)]
                                            +   c23 * A[(i + 1) * NJ * NK + (j + 0) * NK + (k + 1)]
                                            +   c31 * A[(i - 1) * NJ * NK + (j + 1) * NK + (k + 1)]
                                            +   c33 * A[(i + 1) * NJ * NK + (j + 1) * NK + (k + 1)];

                rtTmpAccess(A_OFFSET + (i - 1) * NJ * NK + (j - 1) * NK + (k - 1), 0, 0, idx);
                rtTmpAccess(A_OFFSET + (i + 1) * NJ * NK + (j - 1) * NK + (k - 1), 1, 0, idx);
                rtTmpAccess(A_OFFSET + (i - 1) * NJ * NK + (j - 1) * NK + (k - 1), 2, 0, idx);
                rtTmpAccess(A_OFFSET + (i + 1) * NJ * NK + (j - 1) * NK + (k - 1), 3, 0, idx);
                rtTmpAccess(A_OFFSET + (i - 1) * NJ * NK + (j - 1) * NK + (k - 1), 4, 0, idx);
                rtTmpAccess(A_OFFSET + (i + 1) * NJ * NK + (j - 1) * NK + (k - 1), 5, 0, idx);
                rtTmpAccess(A_OFFSET + (i + 0) * NJ * NK + (j - 1) * NK + (k + 0), 6, 0, idx);
                rtTmpAccess(A_OFFSET + (i + 0) * NJ * NK + (j + 0) * NK + (k + 0), 7, 0, idx);
                rtTmpAccess(A_OFFSET + (i + 0) * NJ * NK + (j + 1) * NK + (k + 0), 8, 0, idx);
                rtTmpAccess(A_OFFSET + (i - 1) * NJ * NK + (j - 1) * NK + (k + 1), 9, 0, idx);
                rtTmpAccess(A_OFFSET + (i + 1) * NJ * NK + (j - 1) * NK + (k + 1), 10, 0, idx);
                rtTmpAccess(A_OFFSET + (i - 1) * NJ * NK + (j + 0) * NK + (k + 1), 11, 0, idx);
                rtTmpAccess(A_OFFSET + (i + 1) * NJ * NK + (j + 0) * NK + (k + 1), 12, 0, idx);
                rtTmpAccess(A_OFFSET + (i - 1) * NJ * NK + (j + 1) * NK + (k + 1), 13, 0, idx);
                rtTmpAccess(A_OFFSET + (i + 1) * NJ * NK + (j + 1) * NK + (k + 1), 14, 0, idx);
                rtTmpAccess(B_OFFSET + (i * NJ * NK + j * NK + k), 15, 1, idx);
            }
        }
    }

}

int main(int argc, char* argv[]) 
{
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
    
    NI = stoi(argv[1]);
    NJ = stoi(argv[2]);
    NK = stoi(argv[3]);
    
    double* A = (double*)malloc( (NI * NJ * NK)*sizeof(double));
    double* B = (double*)malloc( (NI * NJ * NK)*sizeof(double));

    for (int i = 0; i < (NI * NJ * NK); i++) {
            A[i] = i % 256;
    }

    conv3D_trace(A, B);

    string name(argv[0]);
    size_t found = name.find_last_of("/\\") + 1;
    string conf = name.substr(found, name.size()-found) + "_" + to_string(NI) + "_" + to_string(NJ) + "_" + to_string(NK);
    
    dumpRIHistogram(conf);
    predictionWithBmc(conf);
    
    return 0;
}
