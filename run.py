import os
import itertools

bound_values = ["16"]

numOfSymBounds = {"cholesky" : 1, "durbin" : 1, "floyd_warshall" : 1, "gemver" : 1, "gesummv" : 1, "lu" : 1, "ludcmp" : 1, "mvt" : 1, "nussinov" : 1, "stencil" : 1, "trisolv" : 1, "trangle" : 2, "adi" : 2, "atax" : 2, "bicg" : 2, "convolution_2d" : 2, "correlation" : 2, "covariance" : 2, "deriche" : 2, "gramschmidt" : 2, "heat_3d" : 2, "jacobi_1d" : 2, "jacobi_2d" : 2, "seidel_2d" : 2, "symm" : 2, "syr2d" : 2, "syrk" : 2, "trmm" : 2, "convolution_3d" : 3, "doitgen" : 3, "fdtd_2d" : 3, "gemm" : 3, "2mm" : 4, "3mm" : 5}


for bench in numOfSymBounds.keys():
    
    all_permutations = list(itertools.product(bound_values, repeat=numOfSymBounds[bench])) 
    all_permutations = [list(x) for x in all_permutations]
   
    for klevel in [1,2,3,4,5,6,7,8]:
        for tlevel in [1,2,3,4,5,6,7,8]:
        
            print bench, klevel, tlevel    
        
            os.system("echo \"#define CLS 32\" > ./bench/utility/rt_conf.h")
            os.system("echo \"#define DS 8\" >> ./bench/utility/rt_conf.h")
            os.system("echo \"#define KLEVEL " + str(klevel) +  "\" >> ./bench/utility/rt_conf.h")
            os.system("echo \"#define TLEVEL " + str(tlevel) +  "\" >> ./bench/utility/rt_conf.h")
            os.system("make bench_gen " + " BENCH=" + bench)
            
            for bound in all_permutations:
                os.system("./bin/bench_bin/" + bench + "_CLS32_DS8 " + " ".join(bound))


