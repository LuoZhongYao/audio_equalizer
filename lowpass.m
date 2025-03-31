function [ b , a ]  = lowpass( fcutoff , Q )
%
% Derive coefficients for lowpass filter
%
% [ b, a ] = lowpass( fcutoff, Q );
%
% Inputs:
%   fcutoff = normalized cutoff frequency, where fcutoff = 1 at f0 = fs/2;
%   Q  = filter quality
%
% Outputs:
%   [ b, a ] filter coefficients in transfer function format
%
% Remarks:
%   Q = 1/sqrt(2) provides maximally flat pass band up to the cutoff frequency.
%   Q < 1/sqrt(2) provides higher pass band attenuation
%   Q > 1/sqrt(2) provides additional gain around the cutoff frequency

    K = tan(pi*fcutoff/2);

    b0 = ( K.^2 * Q )/( K.^2*Q + K + Q );
    b1 = ( 2*K.^2*Q )/( K.^2*Q + K + Q );
    b2 = b0;

    a0 = 1;
    a1 = 2*Q*(K.^2 - 1)/( K.^2*Q + K + Q );
    a2 = ( K.^2*Q - K + Q )/( K.^2*Q + K + Q );

    b = [ b0, b1, b2 ];
    a = [ a0, a1, a2 ];

end
% All coefficients based on Table 2.2 of DAFX: Digital Audio Effects, 2nd Edition, Udo Zölzer editor
