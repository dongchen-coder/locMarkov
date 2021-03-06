LIB_DIR=./lib
SRC_DIR=./src
OBJ_DIR=./obj
BIN_DIR=./bin

BENCH_DIR=./bench
BENCH_BIN_DIR=./bin/bench_bin

BENCH?=stencil 2mm 3mm adi atax bicg cholesky correlation covariance deriche doitgen durbin fdtd_2d floyd_warshall gemm gemver gesummv gramschmidt heat_3d jacobi_1d jacobi_2d lu ludcmp mvt nussinov seidel_2d symm syr2d syrk trisolv trmm convolution_2d convolution_3d trangle

CC=g++
CCFLAG= -std=c++11 -O2

bench_gen:
	mkdir -p $(BENCH_BIN_DIR)
	$(foreach name, $(BENCH), $(CC) -std=c++11 -O -o $(BENCH_BIN_DIR)/$(name)_CLS32_DS8 $(BENCH_DIR)/$(name).cpp ;)

