%
% Matlab script to generate graphic equalizer filters
%
% Inputs:
%   Filter bit width, 16 or 32
%   Number of frequency bands, 4,5,6,7, or 8
%   Number of biquad IIRs per band, minimum of 2
%   Filename of filter definition .h file
%   Filename of filter definition workspace .mat file
%
% Outputs:
%   Filter definition .h file
%   Filter definition workspace .mat file
%   Plots of filter response, as magnitude and in dB
%
% Filter design parameters in the code:
%    peakAmp = design amplitude for each filter band, specifying the filter
%              response at the center of each band.  Typically the first
%              band and last band must have amplitudes bigger than 1.
%
% Remarks:
%   The total filter response (sum of all band outputs) is calculated and
%   the peak amplitude of this response is normalized to be unity.  This
%   prevents overflows when filtering a signal at 0 dBFS.
%
%   This routine uses the yulewalk filter design function (part of the
%   Signal Processing Toolbox).
clear all;
clc;

% Line specifications for plotting band responses *****************************
plotLineSpecs = { '-r', '-g', 'b', '-c', '-m', '-k', '--r', '--g', '--b', '--c', '--m', '--k' };

% Configuration of dialog window
DialogTitle     = 'Equalization Filter Design';
nDialogLines    = 1;
DialogPrompts   = { 'Bit Width:', '# Bands:', '# Filters/Band:' };
DialogDefaults  = { '16'        , '6'       , '2'               };

% Get user choices
DialogAnswers   = inputdlg(DialogPrompts,DialogTitle,nDialogLines,DialogDefaults);
nFilterBits     = str2num(DialogAnswers{1});
nBands          = str2num(DialogAnswers{2});
nFiltersPerBand = str2num(DialogAnswers{3});

nFilterOrder = 2*nFiltersPerBand;  % Filter order for each band filter

% Band edges, frequency normalized by F0 = FS/2;
fBandEdges = (0:nBands)/nBands;
nFreqs = 512; % Number of frequencies used in [0,Fs/2) to evaluate filter responses
              % Note: nFreqs = nFFT/2


% Filter specification setups *************************************************
switch ( nBands )
    case {4}
        % Frequencies, including band center and band edges
        %          Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge
        fBands = [    0  .125  .250  .375  .500  .625  .750  .875  1.00 ];

        % Band amplitude desired, with 6 dB at band edges
        mBand1 = [    1     1    .5     0     0     0     0     0     0 ];
        mBand2 = [    0     0    .5     1    .5     0     0     0     0 ];
        mBand3 = [    0     0     0     0    .5     1    .5     0     0 ];
        mBand4 = [    0     0     0     0     0     0    .5     1     1 ];

        % Peak amplitude for each band
       %peakAmp = [        +1          +1          +1          +1       ];
        if ( nFiltersPerBand == 3 )
            peakAmp = [    +1.1        +1          +1         +1.1      ];
        elseif ( nFiltersPerBand == 2 )
            peakAmp = [    +1.5        +1          +1         +1.5      ];
        end%if ( nFiltersPerBand == 3 )

        % Sign used for each band output, alternating signs prevents band edges from
        % being 180 out of phase and producing notch in overall filter response.
       %sBands =  [        +1          +1          +1          +1        ];
        sBands =  [        +1          -1          +1          -1        ];

        mBands = [ mBand1; mBand2; mBand3; mBand4 ];

    case {5}
        % Frequencies, including band center and band edges
        %           Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge
        fBands = [    0    .1    .2    .3    .4    .5    .6    .7    .8    .9     1 ];

        % Band amplitude desired, with 6 dB at band edges
        mBand1 = [    1     1    .5     0     0     0     0     0     0     0     0 ];
        mBand2 = [    0     0    .5     1    .5     0     0     0     0     0     0 ];
        mBand3 = [    0     0     0     0    .5     1    .5     0     0     0     0 ];
        mBand4 = [    0     0     0     0     0     0    .5     1    .5     0     0 ];
        mBand5 = [    0     0     0     0     0     0     0     0    .5     1     1 ];

        % Peak amplitude for each band
       %peakAmp =  [       +1          +1          +1          +1          +1       ];
        peakAmp =  [       +1.5        +1          +1          +1          +1.5     ];

        % Sign used for each band output, alternating signs prevents band edges from
        % being 180 out of phase and producing notch in overall filter response.
       %sBands =  [         +1          +1          +1          +1          +1       ];
        sBands =  [         +1          -1          +1          -1          +1       ];

        mBands = [ mBand1; mBand2; mBand3; mBand4; mBand5 ];

    case {6}
        % Frequencies, including band center and band edges
        %            Edge Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge
        fBands = [    0     1     2     3     4     5     6     7     8     9    10    11    12 ]/12;

        % Band amplitude desired, with 6 dB at band edges
        mBand1 = [    1     1    .5     0     0     0     0     0     0     0     0     0     0 ];
        mBand2 = [    0     0    .5     1    .5     0     0     0     0     0     0     0     0 ];
        mBand3 = [    0     0     0     0    .5     1    .5     0     0     0     0     0     0 ];
        mBand4 = [    0     0     0     0     0     0    .5     1    .5     0     0     0     0 ];
        mBand5 = [    0     0     0     0     0     0     0     0    .5     1    .5     0     0 ];
        mBand6 = [    0     0     0     0     0     0     0     0     0     0    .5     1     1 ];

        % Peak amplitude for each band
       %peakAmp =  [       +            +1         +1          +1          +1          +1        ];
        peakAmp =  [      +1.5          +1         +1          +1          +1          +1.5      ];

        % Sign used for each band output, alternating signs prevents band edges from
        % being 180 out of phase and producing notch in overall filter response.
       %sBands =  [        +1          +1          +1          +1          +1          +1       ];
        sBands =  [        +1          -1          +1          -1          +1          -1       ];

        mBands = [ mBand1; mBand2; mBand3; mBand4; mBand5; mBand6 ];

    case {7}
        % Frequencies, including band center and band edges
        %            Edge Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge
        fBands = [    0     1     2     3     4     5     6     7     8     9    10    11    12    13    14 ]/14;

        % Band amplitude desired, with 6 dB at band edges
        mBand1 = [    1     1    .5     0     0     0     0     0     0     0     0     0     0     0     0 ];
        mBand2 = [    0     0    .5     1    .5     0     0     0     0     0     0     0     0     0     0 ];
        mBand3 = [    0     0     0     0    .5     1    .5     0     0     0     0     0     0     0     0 ];
        mBand4 = [    0     0     0     0     0     0    .5     1    .5     0     0     0     0     0     0 ];
        mBand5 = [    0     0     0     0     0     0     0     0    .5     1    .5     0     0     0     0 ];
        mBand6 = [    0     0     0     0     0     0     0     0     0     0    .5     1    .5     0     0 ];
        mBand7 = [    0     0     0     0     0     0     0     0     0     0     0     0    .5     1     1 ];

        % Peak amplitude for each band
       %peakAmp =  [       +1          +1          +1          +1          +1          +1           +1      ];
        peakAmp =  [      +1.5         +1          +1          +1          +1          +1           +1.5    ];

        % Sign used for each band output, alternating signs prevents band edges from
        % being 180 out of phase and producing notch in overall filter response.
       %sBands =      [    +1          +1          +1          +1          +1          +1          +1       ];
        sBands =      [    +1          -1          +1          -1          +1          -1          +1       ];

        mBands = [ mBand1; mBand2; mBand3; mBand4; mBand5; mBand6; mBand7 ];

    case {8}
        % Frequencies, including band center and band edges
        %           Edge Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge  Cntr  Edge
        fBands = [    0 .0625 .1250 .1875 .2500 .3125 .3750 .4375 .5000 .5625 .6250 .6875 .7500 .8125 .8750 .9375 1.000 ];

        % Band amplitude desired, with 6 dB at band edges
        mBand1 = [    1     1    .5     0     0     0     0     0     0     0     0     0     0     0     0     0     0 ];
        mBand2 = [    0     0    .5     1    .5     0     0     0     0     0     0     0     0     0     0     0     0 ];
        mBand3 = [    0     0     0     0    .5     1    .5     0     0     0     0     0     0     0     0     0     0 ];
        mBand4 = [    0     0     0     0     0     0    .5     1    .5     0     0     0     0     0     0     0     0 ];
        mBand5 = [    0     0     0     0     0     0     0     0    .5     1    .5     0     0     0     0     0     0 ];
        mBand6 = [    0     0     0     0     0     0     0     0     0     0    .5     1    .5     0     0     0     0 ];
        mBand7 = [    0     0     0     0     0     0     0     0     0     0     0     0    .5     1    .5     0     0 ];
        mBand8 = [    0     0     0     0     0     0     0     0     0     0     0     0     0     0    .5     1     1 ];

        % Peak amplitude for each band
       %peakAmp =     [    +1           +1          +1          +1          +1          +1          +1          +1       ];
        peakAmp =     [   +1.5          +1          +1          +1          +1          +1          +1         1.5       ];

        % Sign used for each band output, alternating signs prevents band edges from
        % being 180 out of phase and producing notch in overall filter response.
       %sBands =      [    +1           +1          +1          +1          +1          +1           +1         +1       ];
        sBands =      [    +1           -1          +1          -1          +1          -1           +1         -1       ];

        mBands = [ mBand1; mBand2; mBand3; mBand4; mBand5; mBand6; mBand7; mBand8 ];

    otherwise
        error(['Unknown number of bands: ', num2str(nBands)]);

end%switch( nBands )


% Setup variables ************************************************************
tempSOS = tf2sos( ones(1,nFilterOrder+1), ones(1,nFilterOrder+1) );
tempSOS = zeros(size(tempSOS));

SOSBands = cell(nBands,1);
for iBand = 1:nBands
    SOSBands{iBand} = tempSOS;
end%for
GainBands = zeros(nBands,1);

hBandsRaw    = cell(nBands,1);
hBandsScaled = hBandsRaw;


% Calculate coefficients for each band's filter *******************************
for iBand = 1:nBands
    % Use Yule Walker algorithm to calculate filter
    [b,a] = yulewalk(nFilterOrder,fBands,mBands(iBand,:));
    [hBand,wBand] = freqz(b,a,nFreqs); % Get filter response

    [SOS,gSOS] = tf2sos(b,a); % Convert to Biquad coefficients

    % Round to nFilterBits bit coefficients
    SOS = round( (2^(nFilterBits-1))*SOS ) / ( 2^(nFilterBits-1) );

    gSOS = peakAmp(iBand)*gSOS/max(abs(hBand));      % Correct gain
    expGain  = ceil( log2(gSOS) );                   % Calculate gain exponent
    fracGain = round( (2^15) * gSOS / (2^expGain) ); % Gain mantissa is always a Q15 number
    gSOS = (fracGain * (2^expGain))/(2^15);          % Actual gain for this filter

    if ( gSOS < 1 ) % filter gain should be applied at start rather than end to prevent overflows
        SOS(1,1:3) = SOS(1,1:3)*gSOS;
        gSOS = 1;
    end%if ( gSOS < 1 )

    SOSBands{iBand} = SOS;  GainBands(iBand) = gSOS; % Store filter coefficients

    [b,a] = sos2tf(SOSBands{iBand},GainBands(iBand)); % Convert from biquad to transfer function coefficients
    [hBand,wBand] = freqz(b,a,nFreqs); % Calculate filter response
    hBandsRaw{iBand} = sBands(iBand)*hBand; % Store unscaled filter response for later use

    % Accumulate total filter response over all bands
    if ( iBand == 1 )
        hTotal = hBand;
    else
        hTotal = hTotal + sBands(iBand)*hBand;
    end%if ( iBand == 1 )

end

% Scale all band filters to produce peak of unity for total response
AvgResponse  = mean(abs(hTotal));
PeakResponse = max(abs(hTotal));
if ( PeakResponse > 1 ) % Peak Response too big, decrease gain in last biquad by adjusting b's
    for iBand = 1:nBands
        SOS = SOSBands{iBand};
        SOS(nFiltersPerBand,1:3) = SOS(nFiltersPerBand,1:3)/PeakResponse;
        SOSBands{iBand} = SOS;
    end
else % Peak Reponse too small, increase band gains
    for iBand = 1:nBands
        GainBands(iBand) = GainBands(iBand)/PeakResponse;
    end
end%if ( AvgResponse > 1 )

% Calculate scaled filter response and total response
for iBand = 1:nBands
    [b,a] = sos2tf(SOSBands{iBand},GainBands(iBand)); % Convert from biquad to transfer function coefficients
    [hBand,wBand] = freqz(b,a,nFreqs); % Calculate filter response
    hBandsScaled{iBand} = sBands(iBand)*hBand;

    % Accumulate total filter response over all bands
    if ( iBand == 1 )
        hTotal = hBand;
    else
        hTotal = hTotal + sBands(iBand)*hBand;
    end%if ( iBand == 1 )

end%for iBand = 1:nBands

PeakScaledResponse = max(abs(hTotal));
AvgScaledResponse = mean(abs(hTotal));

% Start plot of filter response as filter magnitude ***************************
hFigMag = figure; hold on;
xlabel('Freq/F0');
ylabel('Filter Magnitude');
title([num2str(nFilterOrder), 'th Order Filters over ',num2str(nBands),' bands']);
grid on;
axis([0,1,0,+inf]);
aBandEdges = 0.5*ones(size(fBandEdges));
aBandEdges(1) = 1; aBandEdges(end) = 1;
plot(fBandEdges,aBandEdges,'+','LineWidth',2,'MarkerSize',10);

% Start plot of filter response as dB *****************************************
hFigdB  = figure; hold on;
xlabel('Freq/F0');
ylabel('Filter Response [dB]');
title([num2str(nFilterOrder), 'th Order Filters over ',num2str(nBands),' bands']);
grid on;
axis([0,1,-40,+5]);
aBandEdges = db(0.5,'voltage')*ones(size(fBandEdges));
aBandEdges(1) = 0; aBandEdges(end) = 0;
plot(fBandEdges,aBandEdges,'+','LineWidth',2,'MarkerSize',10);

% Plot filter response for each band ******************************************
LegendText = cell(1,nBands+2);
LegendText{1} = 'Edges';
for iBand = 1:nBands
    LegendText{iBand+1} = ['Band ',num2str(iBand)];
    hBand = hBandsScaled{iBand};
   %hBand = hBandsRaw{iBand};
    iLineSpec = mod(iBand,length(plotLineSpecs)) + 1;
    figure(hFigMag);   plot(wBand/pi,    abs(hBand) , plotLineSpecs{iLineSpec} );
    figure(hFigdB);    plot(wBand/pi, db(abs(hBand)), plotLineSpecs{iLineSpec} );
end %for iBand = 1:nBands

% Finish plot with overall (total) filter response
LegendText{nBands+2} = 'Total';
figure(hFigMag); plot(wBand/pi,   abs(hTotal), '--r','LineWidth',2);
hLegendMag = legend(get(gcf,'CurrentAxes'),LegendText,'Location','Best');
set(hLegendMag,'FontSize',8);
figure(hFigdB);  plot(wBand/pi,db(abs(hTotal)),'--r','LineWidth',2);
hLegenddB = legend(get(gcf,'CurrentAxes'),LegendText,'Location','Best');
set(hLegenddB,'FontSize',8);


% Save workspace for later use ************************************************
DefaultFileName = ['myFilters',num2str(nBands),'x',num2str(nFiltersPerBand),'_Q',num2str(nFilterBits-1),'.mat'];
[matfile, matpath] =  uiputfile({'*.mat'},'Filter Definition Workspace',DefaultFileName);
if ( matfile ~= 0 ) %then - output workspace
    pathFilterWorkspace = [matpath filesep() matfile];
    save(pathFilterWorkspace,...
        'AvgResponse','AvgScaledResponse','GainBands','PeakResponse',...
        'PeakScaledResponse','SOSBands','fBands','fBandEdges',...
        'hBandsRaw','hBandsScaled','hTotal','mBands','nBands','nFilterBits',...
        'nFilterOrder','nFiltersPerBand','nFreqs','peakAmp','sBands','wBand');
end%if ( matfile ~= 0 )


% Output header file to initialize filter array *******************************
DefaultFileName = ['myFilters',num2str(nBands),'x',num2str(nFiltersPerBand),'_Q',num2str(nFilterBits-1),'.h'];
[filterfile, filterpath] =  uiputfile({'*.h'}, 'Filter Definition Header', DefaultFileName);

if ( filterfile ~= 0 ) %then - output filter definition
    pathFilterDefn = [filterpath filesep() filterfile];
    fidFilterDefn = fopen(pathFilterDefn,'w');

    % Print out start of header
    fprintf(fidFilterDefn,'/*******************************************************************************\n');
    fprintf(fidFilterDefn,'  Audio Equalizer Library Filter Definition\n   (Created by GraphicEqualizerFilterDesignScript.m)\n\n');
    fprintf(fidFilterDefn,'  File Name:\n');
    fprintf(fidFilterDefn,'     %s\n\n',filterfile);
    fprintf(fidFilterDefn,'  Summary:\n');
    fprintf(fidFilterDefn,'    %d Bit Filter definition for %d Bands, with %d Filters/Band.\n\n',nFilterBits,nBands,nFiltersPerBand);
    fprintf(fidFilterDefn,'  Description:\n');
    fprintf(fidFilterDefn,'    %d Bit Filter definitions for %d Bands, with %d Filters/Band.\n',nFilterBits,nBands,nFiltersPerBand);
    fprintf(fidFilterDefn,'    See the file %s for details of the filter design.\n',matfile);
    fprintf(fidFilterDefn,'  Remarks:\n<pre>\n');
    fprintf(fidFilterDefn,'    Unscaled Filter Response Peak = %g (%g dB)\n', PeakResponse, 20*log10(PeakResponse));
    fprintf(fidFilterDefn,'    Unscaled Filter Response Average = %g (%g dB)\n', AvgResponse, 20*log10(AvgResponse));
    fprintf(fidFilterDefn,'    Scaled Filter Response Peak = %g (%g dB)\n', PeakScaledResponse, 20*log10(PeakScaledResponse));
    fprintf(fidFilterDefn,'    Scaled Filter Response Average = %g (%g dB)\n', AvgScaledResponse, 20*log10(AvgScaledResponse));
    fprintf(fidFilterDefn,'</pre>\n\n');
    fprintf(fidFilterDefn,'*******************************************************************************/\n');

    % Write out EQUALIZER_FILTER or EQUALIZER_FILTER_32 structure initialization to define filters
    if ( nFilterBits == 16 )
        fprintf(fidFilterDefn,'EQUALIZER_FILTER myFilters[] =\n');
    else
        fprintf(fidFilterDefn,'EQUALIZER_FILTER_32 myFilters[] =\n');
    end%if ( nFilterBits == 16 )
    fprintf(fidFilterDefn,'\t{\n');

    % Write out filters for each band
    for iBand = 1:nBands % for each band...
        % Start definitions for band filters
        fprintf(fidFilterDefn,'\t  //Band %d\n',iBand);

        log2Alpha = 1; % Always one
        SOS = SOSBands{iBand}/(2^log2Alpha); % Scale band filter coefficients
        SOS = round(SOS*2^(nFilterBits-1));  % Convert coefficients to fixed point integers

        % Assign band gain to last biquad in the filter cascade
        gFilter = GainBands(iBand);
        if ( nFilterBits == 16 )
            if ( gFilter ~= 1 ) % Filter gain not unity
                expGain  = max([ceil(log2(gFilter)),0]);
                fracGain = round( 2^15 * gFilter / (2^expGain) ); % Always Q15 number
            else
                expGain = 1;
                fracGain = 2^15/2; % 0.5 * 2^1 = Unity Gain
            end%if ( gFilter ~= 1 )

        else % Filter bits == 32
            if ( gFilter ~= 1 ) % Filter gain not unity
                expGain  = max([ceil(log2(gFilter)),0]);
                fracGain = round( 2^31 * gFilter / (2^expGain) ); % Always Q31 number
            else
                expGain = 1;
                fracGain = 2^31/2;  % 0.5 * 2^1 = Unity Gain
            end%if ( gFilter ~= 1 )

        end%if ( nFilterBits == 16 )

        % Write out filter initialization values ******************************
        %
        %    sBands sign applied to fracGain value written out so that
        %       Yout = Yband1 + Yband2 + Yband3... + YbandN
        %    rather than
        %       Yout = Yband1 - Yband2 + Yband3... +/- YbandN
        %
        for iFilter = 1:nFiltersPerBand % for each filter...

            % Start new filter
            fprintf(fidFilterDefn,'\t\t{ //Filter %d\n',iFilter);

            % Write out {Fractional Gain, Gain Exponent}, Log2Alpha
            if ( iFilter == nFiltersPerBand ) % If last filter in the cascade, use gain
                fprintf(fidFilterDefn,'\t\t\t{%d,%d},%d,\n',sBands(iBand)*fracGain,expGain,log2Alpha);
            else % Otherwise, just set filter's gain to unity, won't actually get used in code
                if ( nFilterBits == 16 )
                    fprintf(fidFilterDefn,'\t\t\t{%d,%d},%d,\n',2^15/2,1,log2Alpha);
                else
                    fprintf(fidFilterDefn,'\t\t\t{%d,%d},%d,\n',2^31/2,1,log2Alpha);
                end%if ( nFilterBits == 16 )
            end%if ( iFilter == nFiltersPerBand )

            % Write out filter taps: {b0,b1,b2}, then {a1,a2} (a1 assumed to always be 1)
            fprintf(fidFilterDefn,'\t\t\t{%d,%d,%d},\n',SOS(iFilter,1),SOS(iFilter,2),SOS(iFilter,3));
            fprintf(fidFilterDefn,'\t\t\t{%d,%d},\n',SOS(iFilter,5),SOS(iFilter,6));

            % Write out initial values for Z's
            fprintf(fidFilterDefn,'\t\t\t{0L,0L}\n');

            % End filter definition
            fprintf(fidFilterDefn,'\t\t}, //end Filter %d\n',iFilter);

        end%for iFilter = 1:nFiltersPerBand

        % End definitions for band filters
        fprintf(fidFilterDefn,'\t  //end Band %d\n\n',iBand);

    end%for iBand = 1:nBands

    fprintf(fidFilterDefn,'\t};\n'); % Last line of initialization.

    fclose(fidFilterDefn); % Close header file, we're done.

end%if ( filename ~= 0 )
