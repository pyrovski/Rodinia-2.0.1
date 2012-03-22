include ../../common/make.config
targets=euler3d  euler3d_double  pre_euler3d  pre_euler3d_double 

all: $(targets)

euler3d: euler3d.cu
	nvcc -Xptxas -v -O3 --gpu-architecture=compute_13 --gpu-code=compute_13 euler3d.cu -o euler3d -I$(CUDA_SDK_PATH)/common/inc  -L$(CUDA_SDK_PATH)/lib $(CUTIL_LIB)

euler3d_double: euler3d_double.cu
	nvcc -Xptxas -v -O3 --gpu-architecture=compute_13 --gpu-code=compute_13 euler3d_double.cu -o euler3d_double -I$(CUDA_SDK_PATH)/common/inc  -L$(CUDA_SDK_PATH)/lib $(CUTIL_LIB)


pre_euler3d: pre_euler3d.cu
	nvcc -Xptxas -v -O3 --gpu-architecture=compute_13 --gpu-code=compute_13 pre_euler3d.cu -o pre_euler3d -I$(CUDA_SDK_PATH)/common/inc  -L$(CUDA_SDK_PATH)/lib $(CUTIL_LIB)

pre_euler3d_double: pre_euler3d_double.cu
	nvcc -Xptxas -v -O3 --gpu-architecture=compute_13 --gpu-code=compute_13 pre_euler3d_double.cu -o pre_euler3d_double -I$(CUDA_SDK_PATH)/common/inc  -L$(CUDA_SDK_PATH)/lib $(CUTIL_LIB)


clean:
	rm -f euler3d euler3d_double pre_euler3d pre_euler3d_double *.linkinfo

install: $(targets)
	install -m $(install_perms) -t $(CUDA_BIN_DIR) $(targets)
