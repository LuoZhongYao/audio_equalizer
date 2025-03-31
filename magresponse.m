function magresponse( b , a , fs )

% Plots the frequency response of any filter defined in terms of b and a
% coefficient vectors.
%
%   magresponse( b , a , fs )
%
% Inputs:
%   b , a = vectors of the coefficients that define the filter
%   fs    = sample rate in samples per second

   im = [ 1 zeros( 1 , fs - 1 ) ] ;
   yim=filter(b,a,im);
   yimf=fft(yim);
   lf=ceil(length(yimf)/2);
   yimf=yimf(1:lf);
   fplot=linspace(0,round(fs/2),lf);
   semilogx(fplot,20*log10(abs(yimf)))
   ylabel('dB')
   xlabel('Frequency (Hz)');
   axis([10 22050 -9 9])
   grid on
   
end

% "Digital Filters for Audio Processing"
% www.osd.rutgers.edu/gs/07papers/dsp.pdf