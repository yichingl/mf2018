#ifndef _PNCC_H_
#define _PNCC_H_

#ifdef WIDEBAND
    #define SAMPRATE    16000.0 /* Sampling frequency */

    #define FFTSIZE  1024 /* FFT point */
    #define NUMFILTS 40   /* Number of filters in the filter bank */

    #define LOWFREQ  100.0  // 200.0
    #define HIGHFREQ 6800.0 // (SAMPRATE / 2.0)
#else
    #define SAMPRATE    8000.0 /* Sampling frequency */
    #define FFTSIZE  512 /* FFT point */
    #define NUMFILTS 28   /* Number of filters in the filter bank */

    #define LOWFREQ  100.0  // 200.0
    #define HIGHFREQ 3700.0 // (SAMPRATE / 2.0)
#endif

#define FRAMELEN    0.0256  /* 25.6 ms window length, */
#define FRAMEPERIOD 0.010   /* 10 ms frame period */

#define PREEMPHASIS 0.97 /* pre-emphasis coefficient */

#define POWERCOEF (1.0 / 15.0)

#define FILTTYPE    2     /* Filter type (2: Gammatone Filter, 1: HTK MEL, 0: Snaley MEL) */
#define GAMMATHRESH 0.005

#define SMTYPE 0 /* Smoothing Type */
#define M      2 /* window length (2M+1) for median-time power */ 
#define N      4 /* window length (2N+1) for spectral smoothing */
    
#define LAMBDA   0.999 /* lambda_a for Q_{le} and Q_f */
#define LAMBDA2  0.5   /* lambda_b for Q_{le} and Q_f */
#define LAMBDA3  0.85  /* lambda_t */
#define DELTA    1
#define DELTA2   0.2   /* mu_t */
#define FACTOR   2.0   /* c */
#define LAMBDA_L 0.999 /* lambda_mu */
#define LAMBDA_S 0.999 /* lambda_mu */

#define CEPDIM 13

#define EPS 2.22e-16

#define STRLEN 256

#endif /* !_PNCC_H_ */
