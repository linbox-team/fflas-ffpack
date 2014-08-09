#include <stdio.h>
#include <cuda_runtime.h>
#include <cusparse.h>

int main() {
	cusparseHandle_t handle = 0;
	cusparseCreate( &handle );
	return 0 ;
}
