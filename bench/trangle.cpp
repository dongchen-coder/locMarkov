#include "./utility/rt.h"

void trangle_trace(int *a, int *b, unsigned int i_size, unsigned int j_size) {
    
    vector<int> idx;
    
    for (int i = 1; i < i_size+1; i++) {
        for (int j = i; j < j_size+1; j++) {
            idx.clear(), idx.push_back(i); idx.push_back(j);
            //b[i* (DIM_SIZE + 2) +j] =  a[i* (DIM_SIZE + 2)+j] + a[i* (DIM_SIZE + 2)+j + 1] + a[i* (DIM_SIZE + 2)+j - 1] + a[(i-1)* (DIM_SIZE + 2) +j] + a[(i+1)* (DIM_SIZE + 2) +j];
            rtTmpAccess(i * (i_size + 2) + j, 0, 0, idx);
            rtTmpAccess(i * (i_size + 2) + j + 1, 1, 0, idx);
            rtTmpAccess(i * (i_size + 2) + j - 1, 2, 0, idx);
            rtTmpAccess( (i-1) * (i_size + 2) + j, 3, 0, idx);
            rtTmpAccess( (i+1) * (i_size + 2) + j, 4, 0, idx);
            rtTmpAccess(i * (i_size + 2) + j + (i_size + 2)* (j_size + 2), 5, 0, idx);
		}
    }
    return;
}



int main(int argc, char* argv[]) {

	int ISIZE=0;
	int JSIZE=0;
	
	if (argc != 3) {
		cout << "need 2 parameters for i, j loops" << endl;
	} else {
		ISIZE = stoi(argv[1]);
		JSIZE = stoi(argv[2]);
	}

    int* a = (int*)malloc( (ISIZE+2)* (JSIZE+2)*sizeof(int));
    int* b = (int*)malloc( (ISIZE+2)* (JSIZE+2)*sizeof(int));

    trangle_trace(a, b, ISIZE, JSIZE);

    dumpRIHistogram();
    
    return 0;
}
