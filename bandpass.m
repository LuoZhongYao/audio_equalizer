function [ b , a ]  = bandpass( Fc, Fbw, Fs )
%
% Derive coefficients for bandpass filter
%
% Usage: [ b, a ] = bandpass( Fc, Fbw, Fs )
%
% Inputs:
%   Fc  = filter center frequency
%   Fbw = filter bandwidth
%   Fs  = Sample rate
%
% Outputs:
%   [ b, a ] filter coefficients in transfer function format
%
% Remarks:
%   None.
%
F0 = Fs/2; %Folding frequency
Wb = Fbw/F0; % Normalized bandwidth
Wc = Fc/F0;  % Normalized center frequency

c = (tan(pi*Wb/2)-1) / (tan(pi*Wb/2)+1); % Eqn 2.27 of DAFX
d = -cos(pi*Wc);                         % Eqn 2.28 of DAFX

b0 = (1+c)/2;
b1 = 0;
b2 = -b0;
b = [ b0, b1, b2 ];

a0 = 1;
a1 = d*(1-c);
a2 = -c;
a = [ a0, a1, a2 ];

return;
%
% Based on M-File 2.2 (apbandpass.m) in DAFX: Digital Audio Effects, 2nd Edition, Udo Zölzer editor
%
