LIB_DIR=./lib
SRC_DIR=./src
OBJ_DIR=./obj
BIN_DIR=./bin

BENCH_DIR=./bench
BENCH_BIN_DIR=./bin/bench_bin

#bench=stencil 2mm 3mm adi atax bicg cholesky correlation covariance deriche doitgen durbin fdtd_2d floyd_warshall gemm gemver gesummv gramschmidt heat_3d jacobi_1d jacobi_2d lu ludcmp mvt nussinov seidel_2d symm syr2d syrk trisolv trmm convolution_2d convolution_3d trangle

bench=2mm

bench_1para= stencil cholesky durbin floyd_warshall gemver gesummv lu ludcmp mvt nussinov stencil trisolv
bench_2para= heat_3d trangle adi atax bicg convolution_2d correlation covariance deriche gramschmidt heat_3d jacobi_1d jacobi_2d seidel_2d symm syr2d syrk trmm
bench_3para= doitgen convolution_3d doitgen fdtd_2d gemm  
bench_4para= 2mm
bench_5para= 3mm

train_size=4 8 12 16 20 24


CC=g++
CCFLAG= -std=c++11 -O2

bench_gen:
	mkdir -p $(BENCH_BIN_DIR)
	echo "#define CLS 32" >> $(BENCH_DIR)/utility/rt_conf.h
	echo "#define DS 8" >> $(BENCH_DIR)/utility/rt_conf.h
	$(foreach name, $(bench), $(CC) -std=c++11 -O -o $(BENCH_BIN_DIR)/$(name)_CLS32_DS8 $(BENCH_DIR)/$(name).cpp ;)

ris_raw_gen:
	#mkdir -p $(RIS_RAW_DIR)_ELM
	#$(foreach name, $(bench), rm -r -f $(RIS_RAW_DIR)_ELM/$(name) ;)
	#$(foreach name, $(bench), mkdir -p $(RIS_RAW_DIR)_ELM/$(name) ;)
	#mkdir -p $(RIS_RAW_DIR)_CLS32_DS8
	#$(foreach name, $(bench), rm -r -f $(RIS_RAW_DIR)_CLS32_DS8/$(name) ;)
	#$(foreach name, $(bench), mkdir -p $(RIS_RAW_DIR)_CLS32_DS8/$(name) ;)
	mkdir -p $(IBOUND_RAW_DIR)
	$(foreach name, $(bench), rm -r -f $(IBOUND_RAW_DIR)/$(name) ;)
	$(foreach name, $(bench), mkdir -p $(IBOUND_RAW_DIR)/$(name) ;)
	$(foreach size, $(train_size), \
		$(foreach name, $(bench_1para), $(BENCH_BIN_DIR)/$(name)_IBOUND $(size) > $(IBOUND_RAW_DIR)/$(name)/$(name)_$(size).txt ;) \
	)
	$(foreach size1, $(train_size), \
		$(foreach size2, $(train_size), \
			$(foreach name, $(bench_2para), $(BENCH_BIN_DIR)/$(name)_IBOUND $(size1) $(size2) > $(IBOUND_RAW_DIR)/$(name)/$(name)_$(size1)_$(size2).txt ;) \
		) \
	)
	$(foreach size1, $(train_size), \
		$(foreach size2, $(train_size), \
			$(foreach size3, $(train_size), \
				$(foreach name, $(bench_3para), $(BENCH_BIN_DIR)/$(name)_IBOUND $(size1) $(size2) $(size3) > $(IBOUND_RAW_DIR)/$(name)/$(name)_$(size1)_$(size2)_$(size3).txt ;) \
			) \
		) \
	)
	$(foreach size1, $(train_size), \
		$(foreach size2, $(train_size), \
			$(foreach size3, $(train_size), \
				$(foreach size4, $(verify_size), \
					$(foreach name, $(bench_4para), $(BENCH_BIN_DIR)/$(name)_IBOUND $(size1) $(size2) $(size3) $(size4) > $(IBOUND_RAW_DIR)/$(name)/$(name)_$(size1)_$(size2)_$(size3)_$(size4).txt ;) \
				) \
			) \
		) \
	)
	$(foreach size1, $(train_size), \
		$(foreach size2, $(train_size), \
			$(foreach size3, $(train_size), \
				$(foreach size4, $(train_size), \
					$(foreach size5, $(train_size), \
						$(foreach name, $(bench_5para), $(BENCH_BIN_DIR)/$(name)_IBOUND $(size1) $(size2) $(size3) $(size4) $(size5) > $(IBOUND_RAW_DIR)/$(name)/$(name)_$(size1)_$(size2)_$(size3)_$(size4)_$(size5).txt ;) \
					) \
				) \
			) \
		) \
	)

