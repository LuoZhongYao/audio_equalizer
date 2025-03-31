function y = highshelving ( x , Wc , G )

% y = highshelving ( x , Wc , G )
%
% Applies a high-frequency shelving filter to the input signal.
%
% Inputs:
%   x  = input vector
%   Wc = normalized cut-off frequency, 0<Wc<1, i.e., fc/(fS/2)
%   G  = gain in dB
%
%   y  = output vector

  V0 = 10 ^ ( G / 20 ) ;
  H0 = V0 - 1 ;
  if G >= 0
    c = ( tan( pi * Wc / 2 ) - 1 ) / ( tan( pi * Wc / 2 ) + 1 ) ;   % boost
  else
    c = ( tan( pi * Wc / 2 ) - 1 ) / ( tan( pi * Wc / 2 ) + 1 ) ; % cut
  end
  xh = 0 ;
  y = zeros( 1 , length( x ) ) ;
  for n = 1:length( x )
    xh_new = x( n ) - c * xh ;
    ap_y = c * xh_new + xh ;
    xh = xh_new ;
    y( n ) = 0.5 * H0 * ( x( n ) - ap_y ) + x( n ) ;
  end

end

% Author: M. Holters
% Modified by Stu Smith
%--------------------------------------------------------------------------
% This source code is provided without any warranties as published in
% DAFX book 2nd edition, copyright Wiley & Sons 2011, available at
% http://www.dafx.de. It may be used for educational purposes and not
% for commercial applications without further permission.
%--------------------------------------------------------------------------