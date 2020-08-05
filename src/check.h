#pragma once
#ifndef CHECK_CUDA
static void checked_cuda(cudaError_t err, const char *file, int line) {
    if (err != cudaSuccess) {
        fprintf(stderr, "%s in %s at line %d\n", cudaGetErrorString(err), file, line );
        exit( EXIT_FAILURE );
    }
}
#define CHECK_CUDA( err ) (checked_cuda( err, __FILE__, __LINE__ ))
#endif

#ifndef CHECK_FREAD
static void checked_fread(void *ptr, size_t size, size_t count, FILE *fp,
        const char *file, int line) {
    size_t actual = fread(ptr, size, count, fp);
    if ( actual != count ) {
        fprintf(stderr, "Read %" PRIu64 " elements but expected %" PRIu64 " (%zu * %zu) in %s at line %d\n",
                actual, count, size, count, file, line);
        exit( EXIT_FAILURE );
    }
}
#define CHECK_FREAD(ptr,size,count,fp) (checked_fread(ptr, size, count, fp, __FILE__, __LINE__)) 
#endif

