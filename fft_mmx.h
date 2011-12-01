#ifndef _QCC_TEST_FFT_H
#define _QCC_TEST_FFT_H

#ifndef LOG2SIZE
#define LOG2SIZE 10	// only even number is allowed. 
#endif
#define FFTSIZE (1<<LOG2SIZE)
#define SQRTSIZE (1<<(LOG2SIZE/2))

#define USING_PREFETCH 1
#define PREFETCH_SIZE 1024
#define MINIMUM_SEGMENT 4096

struct icpx_t {
	signed short re[4],im[4];
};

typedef struct icpx_t icpx_vector[FFTSIZE/4];
typedef struct icpx_t icpx_matrix[SQRTSIZE][SQRTSIZE/4];

struct vcpx_t {
	float re[4],im[4];
};

typedef struct vcpx_t vcpx_vector[FFTSIZE/4];
typedef struct vcpx_t vcpx_matrix[SQRTSIZE][SQRTSIZE/4];

void fft_init(void);
void fft(short f[],icpx_vector g);
void rfft(icpx_vector f,short g[]);
void fft2d(unsigned char f[SQRTSIZE][SQRTSIZE],icpx_matrix g);
void rfft2d(icpx_matrix f,unsigned char g[SQRTSIZE][SQRTSIZE]);

void ffft_init(void);
void ffft(float f[],vcpx_vector g);
void frfft(vcpx_vector f,float g[]);
void ffft2d(unsigned char f[SQRTSIZE][SQRTSIZE],vcpx_matrix g);
void frfft2d(vcpx_matrix f,unsigned char g[SQRTSIZE][SQRTSIZE]);

void differential(vcpx_vector f);
void integral(vcpx_vector f);
void laplace2d(vcpx_matrix f);

#endif

