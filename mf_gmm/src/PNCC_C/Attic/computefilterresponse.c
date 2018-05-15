/*
 * Compute the Gammatone Filter Response
 *
 * Programmed by Kazumasa Yamamoto, May 08, 2012.
 *
 * Based on Chanwoo Kim's MATLAB code
 * for the IEEETran Speech, Audio, and Langauge Processing and ICASSP 2012
 * (chanwook@cs.cmu.edu)
 */

/*
 * This code is based on Snaley's auditory toolbox. Modified by Chanwoo Kim
 */

#include <matrix.h>

matrix ComputeFilterResponse(int iNumFilts, int iFFTSize,
			     double dLowFreq, double dHighFreq, double dSampRate)
{
  matrix aadFilterCoeff = ComputeFilterCoeff(dSampRate, iNumFilts, dLowFreq, dHighFreq);

  for (i = 0; i < iNumFilts; i++) {
    A0[i]   = aadFilterCoeff[i][1];
    A11[i]  = aadFilterCoeff[i][2];
    A12[i]  = aadFilterCoeff[i][3];
    A13[i]  = aadFilterCoeff[i][4];
    A14[i]  = aadFilterCoeff[i][5];
    A2[i]   = aadFilterCoeff[i][6];
    B0[i]   = aadFilterCoeff[i][7];
    B1[i]   = aadFilterCoeff[i][8];
    B2[i]   = aadFilterCoeff[i][9];
    gain[i] = aadFilterCoeff[i][10];
  }

  matrix H = new_matrix(FFT_SIZE / 2, iNumFilts);

  for (chan = 0; chan < iNumFilts; chan++) {

    H1Num[0] = A0[chan] / gain[chan];
    H1Num[1] = A11[chan] / gain[chan];
    H1Num[2] = A2[chan] / gain[chan];
    H1Den[0] = B0[chan];
    H1Den[1] = B1[chan];
    H1Den[2] = B2[chan];

    [H1, W] = freqz(H1Num, H1Den, FFT_SIZE / 2);

    H2Num = [A0(chan) A12(chan) A2(chan)];
    H2Den = [B0(chan) B1(chan) B2(chan)];

    [H2, W] = freqz(H2Num, H2Den, FFT_SIZE / 2);

    H3Num = [A0(chan) A13(chan) A2(chan)];
    H3Den = [B0(chan) B1(chan) B2(chan)];

    [H3, W] = freqz(H3Num, H3Den, FFT_SIZE / 2);

    H4Num = [A0(chan) A14(chan) A2(chan)];
    H4Den = [B0(chan) B1(chan) B2(chan)];

    [H4, W] = freqz(H4Num, H4Den, FFT_SIZE / 2);

    HNum = conv(H1Num, H2Num);
    HNum = conv(HNum,  H3Num);
    HNum = conv(HNum,  H4Num);

    HDen = conv(H1Den, H2Den);
    HDen = conv(H1Den, H3Den);
    HDen = conv(H1Den, H4Den);

    iFiltIndex = length(gain) + 1 - chan;

    H(:, iFiltIndex)  = H1 .* H2 .* H3 .* H4;
  }
}


/*
 * This function computes the filter coefficients for a bank of
 * Gammatone filters.  These filters were defined by Patterson and
 * Holdworth for simulating the cochlea.
 *
 * The result is returned as an array of filter coefficients.  Each row
 * of the filter arrays contains the coefficients for four second order
 * filters.  The transfer function for these four filters share the same
 * denominator (poles) but have different numerators (zeros).  All of these
 * coefficients are assembled into one vector that the ERBFilterBank
 * can take apart to implement the filter.
 *
 * The filter bank contains "numChannels" channels that extend from
 * half the sampling rate (fs) to "lowFreq".  Alternatively, if the numChannels
 * input argument is a vector, then the values of this vector are taken to
 * be the center frequency of each desired filter.  (The lowFreq argument is
 * ignored in this case.)
 *
 * Note this implementation fixes a problem in the original code by
 * computing four separate second order filters.  This avoids a big
 * problem with round off errors in cases of very small cfs (100Hz) and
 * large sample rates (44kHz).  The problem is caused by roundoff error
 * when a number of poles are combined, all very close to the unit
 * circle.  Small errors in the eigth order coefficient, are multiplied
 * when the eigth root is taken to give the pole location.  These small
 * errors lead to poles outside the unit circle and instability.  Thanks
 * to Julius Smith for leading me to the proper explanation.
 *
 * Execute the following code to evaluate the frequency
 * response of a 10 channel filterbank.
 *       fcoefs = MakeERBFilters(16000,10,100);
 *       y = ERBFilterBank([1 zeros(1,511)], fcoefs);
 *       resp = 20*log10(abs(fft(y')));
 *       freqScale = (0:511)/512*16000;
 *       semilogx(freqScale(1:255),resp(1:255,:));
 *       axis([100 16000 -60 0])
 *       xlabel('Frequency (Hz)'); ylabel('Filter Response (dB)');
 *
 * Rewritten by Malcolm Slaney@Interval.  June 11, 1998.
 * (c) 1998 Interval Research Corporation
 */
matrix ComputeFilterCoeff(double dSampRate, int iNumChannels, double dLowFreq, double dHighFreq)
{
  double T = 1.0 / dSampRate;

  // This is just exactly the same as ERBSpace in Snaley's toolbox
  vector cf = ComputeERBSpace(dLowFreq, dHighFreq, iNumChannels);

  /*
   * Change the followFreqing three parameters if you wish to use a different
   * ERB scale.  Must change in ERBSpace too.
   */
  double EarQ  = 9.26449; // Glasberg and Moore Parameters
  double minBW = 24.7;
  int    oder  = 1;

  matrix aadFilterCoeff = new_matrix(iNumChannels, 10);

  A0 = T;
  A2 = 0.0;
  B0 = 1.0;
  for (i = 0; i < iNumChannels; i++) {
    ERB[i] = pow(pow((cf[i] / EarQ), order) + pow(minBW, order), 1.0 / order);
    B[i]   = 1.019 * 2.0 * M_PI * ERB[i];
    B1[i] = -2.0 * cos(2.0 * cf[i] * pi * T) / exp(B[i] * T);
    B2[i] = exp(-2.0 * B[i] * T);

    A11[i] = -(2.0 * T * cos(2.0 * cf[i] * M_PI * T) / exp(B[i] * T)
	       + 2.0 * sqrt(3.0 + pow(2.0, 1.5)) * T
	       * sin(2.0 * cf[i] * M_PI * T) / exp(B[i] * T)) /2.0;
    A12[i] = -(2.0 * T * cos(2.0 * cf[i] * M_PI * T) / exp(B[i] * T)
	       - 2.0 * sqrt(3.0 - pow(2.0, 1.5)) * T
	       * sin(2.0 * cf[i] * M_PI * T) / exp(B[i] * T)) /2.0;
    A13[i] = -(2.0 * T * cos(2.0 * cf[i] * M_PI * T) / exp(B[i] * T)
	       + 2.0 * sqrt(3.0 + pow(2.0, 1.5)) * T
	       * sin(2.0 * cf[i] * M_PI * T) / exp(B[i] * T)) /2.0;
    A14[i] = -(2.0 * T * cos(2.0 * cf[i] * M_PI * T) / exp(B[i] * T)
	       - 2.0 * sqrt(3.0 - pow(2.0, 1.5)) * T
	       * sin(2.0 * cf[i] * M_PI * T) / exp(B[i] * T)) /2.0;

    gain[i] = ;
    /* "i" is sqrt(-1)
    abs((-2.0 * exp(4.0 * i * cf * pi * T) * T +
		2 * exp(-(B * T) + 2 * i * cf * pi * T).*T.*
		(cos(2*cf*pi*T) - sqrt(3 - 2^(3/2)) *
		 sin(2*cf*pi*T))) .*
	       (-2*exp(4*i*cf*pi*T)*T +
		2*exp(-(B*T) + 2*i*cf*pi*T).*T.*
		(cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) *
		 sin(2*cf*pi*T))).*
	       (-2*exp(4*i*cf*pi*T)*T +
		2*exp(-(B*T) + 2*i*cf*pi*T).*T.*
		(cos(2*cf*pi*T) -
		 sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .*
	       (-2*exp(4*i*cf*pi*T)*T + 2*exp(-(B*T) + 2*i*cf*pi*T).*T.*
		(cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./
	       (-2 ./ exp(2*B*T) - 2*exp(4*i*cf*pi*T) +
		2*(1 + exp(4*i*cf*pi*T))./exp(B*T)).^4);
    */

    aadFilterCoeff[i][0] = A0;
    aadFilterCoeff[i][1] = A11[i];
    aadFilterCoeff[i][2] = A12[i];
    aadFilterCoeff[i][3] = A13[i];
    aadFilterCoeff[i][4] = A14[i];
    aadFilterCoeff[i][5] = A2;
    aadFilterCoeff[i][6] = B0;
    aadFilterCoeff[i][7] = B1[i];
    aadFilterCoeff[i][8] = B2[i];
    aadFilterCoeff[i][9] = gain[i];
  }

  return aadFilterCoeff;
}


/*
 * This function computes an array of N frequencies uniformly spaced between
 * highFreq and lowFreq on an ERB scale.  N is set to 100 if not specified.
 *
 * See also linspace, logspace, MakeERBCoeffs, MakeERBFilters.
 *
 * For a definition of ERB, see Moore, B. C. J., and Glasberg, B. R. (1983).
 * "Suggested formulae for calculating auditory-filter bandwidths and
 * excitation patterns," J. Acoust. Soc. Am. 74, 750-753.
 */
vector ComputeERBSpace(double lowFreq, double highFreq, int N)
{
  /*
   * Change the following three parameters if you wish to use a different
   * ERB scale.  Must change in MakeERBCoeffs too.
   */
  double EarQ  = 9.26449; // Glasberg and Moore Parameters
  double minBW = 24.7;
  int    order = 1;

  vector cfArray = new_vector(N);

  /*
   * All of the followFreqing expressions are derived in Apple TR #35, "An
   * Efficient Implementation of the Patterson-Holdsworth Cochlear
   * Filter Bank."  See pages 33-34.
   */
  for (i = 0; i < N; i++) {
    cfArray[i] = -(EarQ * minBW)
      + exp((i + 1) * (-log(highFreq + EarQ * minBW) + log(lowFreq + EarQ * minBW)) / N)
      * (highFreq + EarQ * minBW);
  }

  return cfArray;
}


/*
function [H] = ComputeFilterResponse(iNumFilts, iFFTSize, dLowFreq, dHighFreq, dSampRate)

    iNumChannels = iNumFilts;
    FFT_SIZE     = iFFTSize;

    aadFilterCoeff = ComputerFilterCoeff(dSampRate,iNumChannels, dLowFreq, dHighFreq);

    A0  = aadFilterCoeff(:,1);
    A11 = aadFilterCoeff(:,2);
    A12 = aadFilterCoeff(:,3);
    A13 = aadFilterCoeff(:,4);
    A14 = aadFilterCoeff(:,5);
    A2  = aadFilterCoeff(:,6);
    B0  = aadFilterCoeff(:,7);
    B1  = aadFilterCoeff(:,8);
    B2  = aadFilterCoeff(:,9);
    gain= aadFilterCoeff(:,10);


    H = zeros(FFT_SIZE / 2, length(gain));

    for chan = 1 : length(gain)

        H1Num = [A0(chan)/gain(chan) A11(chan)/gain(chan) ...
               A2(chan)/gain(chan)];
        H1Den = [B0(chan) B1(chan) B2(chan)];

        [H1, W] = freqz(H1Num, H1Den, FFT_SIZE / 2);

        H2Num = [A0(chan) A12(chan) A2(chan)];
        H2Den = [B0(chan) B1(chan) B2(chan)];

        [H2, W] = freqz(H2Num, H2Den, FFT_SIZE / 2);

        H3Num = [A0(chan) A13(chan) A2(chan)];
        H3Den = [B0(chan) B1(chan) B2(chan)];

        [H3, W] = freqz(H3Num, H3Den, FFT_SIZE / 2);

        H4Num = [A0(chan) A14(chan) A2(chan)];
        H4Den = [B0(chan) B1(chan) B2(chan)];

        [H4, W] = freqz(H4Num, H4Den, FFT_SIZE / 2);

        HNum = conv(H1Num, H2Num);
        HNum = conv(HNum,  H3Num);
        HNum = conv(HNum,  H4Num);

        HDen = conv(H1Den, H2Den);
        HDen = conv(H1Den, H3Den);
        HDen = conv(H1Den, H4Den);


       % hold on
      %  freqz(HNum, HDen);
        iFiltIndex = length(gain) + 1 - chan;

        H(:, iFiltIndex)  = H1 .* H2 .* H3 .* H4;

     % plot(linspace(0, 8000, length(H)), 20 * log10(abs(H)));
     % plot(linspace(0, 8000, length(H)), (abs(H)));



     % hold on

    end

    %axis([0 8000 -100 20])

 %   xlabel('Frequency (Hz)');
 %   ylabel('Frequency Response (dB)')

 %   grid on

end

















function [aadFilterCoeff] = ComputerFilterCoeff(dSampRate,iNumChannels, dLowFreq, dHighFreq)
% function [fcoefs]=MakeERBFilters(fs,numChannels,lowFreq)
% This function computes the filter coefficients for a bank of
% Gammatone filters.  These filters were defined by Patterson and
% Holdworth for simulating the cochlea.
%
% The result is returned as an array of filter coefficients.  Each row
% of the filter arrays contains the coefficients for four second order
% filters.  The transfer function for these four filters share the same
% denominator (poles) but have different numerators (zeros).  All of these
% coefficients are assembled into one vector that the ERBFilterBank
% can take apart to implement the filter.
%
% The filter bank contains "numChannels" channels that extend from
% half the sampling rate (fs) to "lowFreq".  Alternatively, if the numChannels
% input argument is a vector, then the values of this vector are taken to
% be the center frequency of each desired filter.  (The lowFreq argument is
% ignored in this case.)

% Note this implementation fixes a problem in the original code by
% computing four separate second order filters.  This avoids a big
% problem with round off errors in cases of very small cfs (100Hz) and
% large sample rates (44kHz).  The problem is caused by roundoff error
% when a number of poles are combined, all very close to the unit
% circle.  Small errors in the eigth order coefficient, are multiplied
% when the eigth root is taken to give the pole location.  These small
% errors lead to poles outside the unit circle and instability.  Thanks
% to Julius Smith for leading me to the proper explanation.

% Execute the following code to evaluate the frequency
% response of a 10 channel filterbank.
%       fcoefs = MakeERBFilters(16000,10,100);
%       y = ERBFilterBank([1 zeros(1,511)], fcoefs);
%       resp = 20*log10(abs(fft(y')));
%       freqScale = (0:511)/512*16000;
%       semilogx(freqScale(1:255),resp(1:255,:));
%       axis([100 16000 -60 0])
%       xlabel('Frequency (Hz)'); ylabel('Filter Response (dB)');

% Rewritten by Malcolm Slaney@Interval.  June 11, 1998.
% (c) 1998 Interval Research Corporation

    T = 1/dSampRate;
    if length(iNumChannels) == 1
          % This is just exactly the same as ERBSpace in Snaley's toolbox
        cfArray = ComputeERBSpace(dLowFreq, dHighFreq, iNumChannels);
        cf      = cfArray;
    else
        cf = numChannels(1:end);
        if size(cf,2) > size(cf,1)
            cf = cf';
        end
    end

    % Change the followFreqing three parameters if you wish to use a different
    % ERB scale.  Must change in ERBSpace too.
    EarQ = 9.26449;                             %  Glasberg and Moore Parameters
    minBW = 24.7;
    order = 1;

    ERB = ((cf/EarQ).^order + minBW^order).^(1/order);
    B=1.019*2*pi*ERB;

    A0 = T;
    A2 = 0;
    B0 = 1;
    B1 = -2*cos(2*cf*pi*T)./exp(B*T);
    B2 = exp(-2*B*T);

    A11 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./ ...
            exp(B*T))/2;
    A12 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./ ...
            exp(B*T))/2;
    A13 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./ ...
            exp(B*T))/2;
    A14 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./ ...
            exp(B*T))/2;

    gain = abs((-2*exp(4*i*cf*pi*T)*T + ...
                     2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
                             (cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))* ...
                              sin(2*cf*pi*T))) .* ...
               (-2*exp(4*i*cf*pi*T)*T + ...
                 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
                  (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) * ...
                   sin(2*cf*pi*T))).* ...
               (-2*exp(4*i*cf*pi*T)*T + ...
                 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
                  (cos(2*cf*pi*T) - ...
                   sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .* ...
               (-2*exp(4*i*cf*pi*T)*T + 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
               (cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./ ...
              (-2 ./ exp(2*B*T) - 2*exp(4*i*cf*pi*T) +  ...
               2*(1 + exp(4*i*cf*pi*T))./exp(B*T)).^4);

    allfilts = ones(length(cf),1);
    aadFilterCoeff = [A0*allfilts A11 A12 A13 A14 A2*allfilts B0*allfilts B1 B2 gain];

end

function cfArray = ComputeERBSpace(lowFreq, highFreq, N)
% function cfArray = ERBSpace(lowFreq, highFreq, N)
% This function computes an array of N frequencies uniformly spaced between
% highFreq and lowFreq on an ERB scale.  N is set to 100 if not specified.
%
% See also linspace, logspace, MakeERBCoeffs, MakeERBFilters.
%
% For a definition of ERB, see Moore, B. C. J., and Glasberg, B. R. (1983).
% "Suggested formulae for calculating auditory-filter bandwidths and
% excitation patterns," J. Acoust. Soc. Am. 74, 750-753.

    if nargin < 1
        lowFreq = 100;
    end

    if nargin < 2
        highFreq = 44100/4;
    end

    if nargin < 3
        N = 100;
    end

    % Change the following three parameters if you wish to use a different
    % ERB scale.  Must change in MakeERBCoeffs too.
    EarQ = 9.26449;                             %  Glasberg and Moore Parameters
    minBW = 24.7;
    order = 1;

    % All of the followFreqing expressions are derived in Apple TR #35, "An
    % Efficient Implementation of the Patterson-Holdsworth Cochlear
    % Filter Bank."  See pages 33-34.
    cfArray = -(EarQ*minBW) + exp((1:N)'*(-log(highFreq + EarQ*minBW) + ...
            log(lowFreq + EarQ*minBW))/N) * (highFreq + EarQ*minBW);

end
*/
