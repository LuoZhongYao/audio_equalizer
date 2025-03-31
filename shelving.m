function [ b , a ]  = shelving( type , g , fc , Q , sr )
%
% Derive coefficients for a shelving filter with a given amplitude and
% cutoff frequency. Use MATLAB filter() to apply filter to signal. See
% ShelvingDemo for examples of the use of shelving.
%
%   [ b , a ]  = shelving( type , g , fc , Q , sr )
%
% Inputs:
%   type = a character string defining filter type.
%          Choices are: 'Bass_Shelf' or 'Treble_Shelf'
%   g    = the logarithmic gain in dB; positive for boost, negative for cut
%   fc   = the center frequency
%   Q    = adjusts the slope of shelving
%   sr   = sample rate in samples per second
%
% Outputs:
%  [ b , a ] = filter coefficients to be used with MATLAB filter()

% Error Check
  if((strcmp(type,'Bass_Shelf') ~= 1) && (strcmp(type,'Treble_Shelf') ~= 1))
    error(['Unsupported Filter Type: ' type]);
  end
%
  K = tan((pi * fc)/sr);
  V0 = 10^(g/20);
  root2 = 1/Q; %sqrt(2)
%
% Invert gain if a cut
  if(V0 < 1)
    V0 = 1/V0;
  end
%
%%%%%%%%%%%%%%%%%%%%
%    BASS BOOST
%%%%%%%%%%%%%%%%%%%%
%
  if(( g > 0 ) && (strcmp(type,'Bass_Shelf')))

    b0 = (1 + sqrt(V0)*root2*K + V0*K^2) / (1 + root2*K + K^2);
    b1 =             (2 * (V0*K^2 - 1) ) / (1 + root2*K + K^2);
    b2 = (1 - sqrt(V0)*root2*K + V0*K^2) / (1 + root2*K + K^2);
    a1 =                (2 * (K^2 - 1) ) / (1 + root2*K + K^2);
    a2 =             (1 - root2*K + K^2) / (1 + root2*K + K^2);
%
%%%%%%%%%%%%%%%%%%%%
%    BASS CUT
%%%%%%%%%%%%%%%%%%%%
%
  elseif (( g < 0 ) && (strcmp(type,'Bass_Shelf')))

    b0 =             (1 + root2*K + K^2) / (1 + root2*sqrt(V0)*K + V0*K^2);
    b1 =                (2 * (K^2 - 1) ) / (1 + root2*sqrt(V0)*K + V0*K^2);
    b2 =             (1 - root2*K + K^2) / (1 + root2*sqrt(V0)*K + V0*K^2);
    a1 =             (2 * (V0*K^2 - 1) ) / (1 + root2*sqrt(V0)*K + V0*K^2);
    a2 = (1 - root2*sqrt(V0)*K + V0*K^2) / (1 + root2*sqrt(V0)*K + V0*K^2);
%
%%%%%%%%%%%%%%%%%%%%
%   TREBLE BOOST
%%%%%%%%%%%%%%%%%%%%
%
  elseif (( g > 0 ) && (strcmp(type,'Treble_Shelf')))

    b0 = (V0 + root2*sqrt(V0)*K + K^2) / (1 + root2*K + K^2);
    b1 =             (2 * (K^2 - V0) ) / (1 + root2*K + K^2);
    b2 = (V0 - root2*sqrt(V0)*K + K^2) / (1 + root2*K + K^2);
    a1 =              (2 * (K^2 - 1) ) / (1 + root2*K + K^2);
    a2 =           (1 - root2*K + K^2) / (1 + root2*K + K^2);
%
%%%%%%%%%%%%%%%%%%%%
%   TREBLE CUT
%%%%%%%%%%%%%%%%%%%%
%
  elseif (( g < 0 ) && (strcmp(type,'Treble_Shelf')))

    b0 =               (1 + root2*K + K^2) / (V0 + root2*sqrt(V0)*K + K^2);
    b1 =                  (2 * (K^2 - 1) ) / (V0 + root2*sqrt(V0)*K + K^2);
    b2 =               (1 - root2*K + K^2) / (V0 + root2*sqrt(V0)*K + K^2);
    a1 =             (2 * ((K^2)/V0 - 1) ) / (1 + root2/sqrt(V0)*K + (K^2)/V0);
    a2 = (1 - root2/sqrt(V0)*K + (K^2)/V0) / (1 + root2/sqrt(V0)*K + (K^2)/V0);

%%%%%%%%%%%%%%%%%%%%
%   All-Pass
%%%%%%%%%%%%%%%%%%%%
  else
    b0 = V0;
    b1 = 0;
    b2 = 0;
    a1 = 0;
    a2 = 0;
  end
%
% return values
  a = [  1, a1, a2];
  b = [ b0, b1, b2];

end
%
%
% All coefficients are calculated as described in DAFX (p. 50-55).
% Author:    Jeff Tackett 08/22/05