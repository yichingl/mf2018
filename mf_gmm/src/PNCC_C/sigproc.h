#ifndef _SIGPROC_H_
#define _SIGPROC_H_

/* Windowing functions */
void Hanning(double *data, int windowSize);
void Hamming(double *data, int windowSize);

/* FFT */
int fft(int windowSize, double *realData, double *imaginaryData);

/* DCT */
int dct(int in_dim, double *in_data, int out_dim, double *out_data);

#endif /* !_SIGPROC_H_ */
