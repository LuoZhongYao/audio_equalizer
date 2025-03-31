function [ b a ] = parmeq( dB , cf , bw , sr  )

% Second-order parametric EQ filter design. Computes the b and a
% coefficients to be passed to the Matlab "filter" function. This
% function is intended to be used as a stage in a parametric EQ.
%
%   [ b a ] = parmeq( dB , cf , bw , sr )
%
% Inputs:
%   dB = boost/cut (decibels)
%   cf = center frequency in Hz.
%   bw = bandwidth in Hz.
%   sr = sample rate in samples per second
%
% Outputs:
%   b = numerator coefficients
%   a = denominator coefficients

  G0 = 1 ;
  GB = 1 / sqrt( 2 ) ;
  G  = 10 ^ ( dB / 20 ) ;
  w0 = 2 * pi * cf / sr ;
  Dw = 2 * pi * bw / sr ;
  beta = tan(Dw/2) * sqrt(abs(GB^2 - G0^2)) / sqrt(abs(G^2 - GB^2));
  b = [(G0 + G*beta), -2*G0*cos(w0), (G0 - G*beta)] / (1+beta);
  a = [1, -2*cos(w0)/(1+beta), (1-beta)/(1+beta)];

end

% "Digital Filters for Audio Processing"
% www.osd.rutgers.edu/gs/07papers/dsp.pdf
%
% see also:
% Sophocles J. Orfanidis. "Introduction to Signal Processing".
% Prentice-Hall, 1996, p. 766.
