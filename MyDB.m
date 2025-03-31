function [ Y ] = MyDB( X, Units )
%
% Matlab function to compute amplitude, in dB, of input signal X.
%
%  Usage: 
%   db_amp = MyDB( X, Units );
%
%  Input:
%   X = complex transfer function or time series
%   Units = 'power' or 'voltage', case insensitive
%
%  Output:
%   db_amp = 10.0*log10(abs(X)), if Units == 'power'
%   db_amp = 10.0*log10(abs(X)), if Units == 'voltage'
%
if ( nargin == 0 ) %then
    error('Wrong number of arguments! See help MyDB for more info.')
elseif ( nargin==1) %then
   Units = 'voltage';
end%if ( nargin = 0 ) %then

switch ( lower(Units) )
    case {'voltage'}
        Y = ( 20.*log10(abs(X)) + 300 ) - 300;
         
    case {'power'}
        Y = ( 10.*log10(abs(X)) + 300 ) - 300;
        
    otherwise
        error(['Unrecognized Units: ',Units]);
        
end%switch ( lower(Units) )

end%function 

