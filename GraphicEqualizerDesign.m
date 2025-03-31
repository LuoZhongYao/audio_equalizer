%
% Matlab script to generate graphic equalizer filters
%
% Inputs:
%   Filter bit width, 16 or 32
%   Sample Rate, in Hz
%   Number of IIR Biquad Filters Per Band
%   Filter Center Frequencies, in Hz
%   Filter Gains, in dB, at center frequencies
%   Filter Bandwidths, in Hz
%   Shelving Fitler Qs, for low and high filters
%   Plot command, plot for linear frequency axis, semilogx for logarithmic
%   Filename of filter definition .h file
%   Filename of filter definition workspace .mat file
%
% Outputs:
%   Filter definition .h file
%   Filter definition workspace .mat file
%   Plot of overall filter responses, in dB
%
% Remarks:
%   Running this script will generate a input dialog window, with pre-wired
%   defaults if the script is run from an empty workspace.  Loading a
%   previously generated filter design workspace .mat file will restore the
%   previously used settings.
%
%   After design filters a question dialog window appears asking if the
%   filter design is done.  "Yes" moves the script onto saving the results,
%   "No" brings up another input dialog window, allowing the user to
%   iterate on filter design parameters.  The previously generated plot
%   remains on the screen, allowing the user to compare old versus new.
%
%   When filter design is done, and after saving the results, a question
%   dialog window appears asking whether to erase the current workspace.
%   Select "Yes" to clear the current workspace, "No" if you wish to rerun
%   the script and continue designing filters.
%

% Determine tool being used.
VerStructure = ver;
nVersions = length(VerStructure);
GotMatlab = false;
for ( iVer = 1:nVersions )
    if ( strcmp( upper(VerStructure(iVer).Name), 'MATLAB') )
        GotMatlab = true;
        break;
    end%if
end%for

% Line specifications for plotting filter responses *****************************
plotLineSpecs = { '-r', '-g', 'b', '-c', '-m', '--r', '--g', '--b', '--c', '--m' };


% Parameters ******************************************************************
nFreqs = 2048; % Number of frequencies used in [0,Fs/2] to evaluate filter responses

% Setup defaults for diaglog window *******************************************
Q=(sqrt(2)/2)
if ( ~exist('FilterDesignID','var') ) %then
    FilterDesignID = 0;
end%if
if ( ~exist('nFilterBits','var') ) %then
    nFilterBits = 16;
end%if
if ( ~exist('Fs','var') ) %then
    Fs = 48000; % Default sample rate, in Hz
end%if
if ( ~exist('nFiltersPerBand','var') )
    nFiltersPerBand = 2; % Two biquads per band
end%if
if ( ~exist('FcFilters','var') ) %then
    nFilters = 10;
    DeltaFc = 1/(nFilters-1);
    FcFilters = [    32    64    125    250    500    1000    2000    4000    8000    16000]; % Default filter center frequencies
end%if
if ( ~exist('DesiredGains','var') ) %then
    DesiredGains = [7    7    4    1     0    4    1    5    7    5]; % Default filter gains, in dB;
end%if
if ( ~exist('BWidthFilters','var') )%then
    BWidthFilters = FcFilters/Q; % Bandwidth for filters.
end%if
if ( ~exist('SignBands','var') ) %then
    SignBands = (2*mod(1:nFilters,2)-1);
end%if
if ( ~exist('Qs','var') ) %then
    Qs = (sqrt(2)/2)*[1,1]; % Default Shelving Q factors
end%if
if ( ~exist('PlotCommand','var') ) %then
    PlotCommand = 'semilogx'; % Default: logarithmic frequency axis
end%if


% Iterate filter design until satisfied ***************************************
bDone = false;
while (~bDone)
    FilterDesignID = FilterDesignID + 1; % Starting new design iteration

    % Setup configuration dialog window
    DialogTitle     = ['Graphic Equalization Filter Design ',num2str(FilterDesignID)];
    DialogPrompts   = { 'Filter_Bit_Width:', ...
                        'Sample_Rate_[Hz]:', ...
                        '#_Filters/Band:', ...
                        'Center_Frequencies_[Hz]:', ...
                        'Gains_[dB]:', ...
                        'Bandwidths_[Hz]:',...
                        'Signs [+1/-1]:', ...
                        'Qs_[LowPass,HighPass]:', ...
                        'Plot_Command_[plot,semilogx]:'
                      };
    FcFilterDefaults   = sprintf('%10g\t',FcFilters);
    GainFilterDefaults = sprintf('%10g\t',DesiredGains);
    BandWidthDefaults  = sprintf('%10g\t',BWidthFilters);
    SignBandsDefaults =sprintf('%10g\t',SignBands);
    QsDefault   = sprintf('%10g',Qs(1));
    for iBand = 2:(length(FcFilters)-1)
        QsDefault = [QsDefault,sprintf('\t        NA')];
    end%for
    QsDefault = [QsDefault,sprintf('\t%10g',Qs(2))];
    nColumns = max( [length(FcFilterDefaults), ...
                     length(GainFilterDefaults),  ...
                     length(BandWidthDefaults), ...
                     length(QsDefault)] );
    DialogDefaults = { num2str(nFilterBits), ...
                       num2str(Fs), ...
                       num2str(nFiltersPerBand), ...
                       FcFilterDefaults,...
                       GainFilterDefaults,...
                       BandWidthDefaults,...
                       SignBandsDefaults, ...
                       QsDefault,...
                       PlotCommand ...
                     };
    options.WindowStyle='normal';
    %             Row, Columns
    rowcolumns = [  1, length('Filter Bit Width:'); ... %nFilterbits
                    1, length('Sample_Rate_[Hz]:'); ... %Fs
                    1, length('#_Filters/Band:');   ... %# Filters/Band:
                    1, nColumns; ... % Center frequencies
                    1, nColumns; ... % Gains
                    1, nColumns; ... % Bandwidths
                    1, nColumns; ... % Filter Signs
                    1, nColumns; ... % Shelving Qs
                    1, length('Plot Command [plot,semilogx]:') % Plot command: linear or log frequency axis
                 ];

    % Get user choices
    DialogAnswers = inputdlg(DialogPrompts,DialogTitle,rowcolumns,DialogDefaults,options);
    assert( length(DialogAnswers) == 9,'Dialog terminated! Dialog Length: %d',length(DialogAnswers));
    nFilterBits   = str2num(DialogAnswers{1});
    Fs            = str2num(DialogAnswers{2});
    nFiltersPerBand = str2num(DialogAnswers{3});
    FcFilters     = str2num(DialogAnswers{4});
    DesiredGains  = str2num(DialogAnswers{5});
    BWidthFilters = str2num(DialogAnswers{6});
    SignBands     = str2num(DialogAnswers{7});
    Qs            = str2num(strrep(DialogAnswers{8},'NA',''));
    PlotCommand   = DialogAnswers{9};

    % Validate user input
    assert( Fs > 0,'Sampling is not positive: %g',Fs)

    nFcFilters   = length(FcFilters);
    nGainFilters = length(DesiredGains);
    nBandWidths = length(BWidthFilters);
    if ( nFcFilters ~= nGainFilters ) %then
        error(['Number of frequencies (',num2str(nFcFilters),') doesn''t match gains specified (',num2str(nGainFilters),')!']);
    end%if
    if ( nFcFilters ~= nBandWidths ) %then
        error(['Number of interior frequencies (',num2str(nFcFilters-2),') doesn''t match bandwidths specified (',num2str(nBandWidths),')!']);
    end%if

    switch ( lower(PlotCommand) )
        case 'plot'
            %Do Nothing
        case 'semilogx'
            if ( ~GotMatlab ) %then - Disable warnings about plotting 0 Hz on log freq axis
                warning('off','all');
            end%if
        otherwise
            error(['Unknown plot command: ',PlotCommand]);
    end%switch


    % Setup variables for later use *******************************************
    F0 = Fs/2;
    nBands = nFcFilters;
    nFilterOrder = 2*nFiltersPerBand; % All filters are single biquad IIRs

    peakAmp = (10*ones(1,nBands)).^(DesiredGains/20);
    tempSOS = tf2sos( ones(1,nFilterOrder+1), ones(1,nFilterOrder+1) );
    tempSOS = zeros(size(tempSOS));

    SOSBands = cell(nBands,1);
    for iBand = 1:nBands
        SOSBands{iBand} = tempSOS;
    end%for
    GainBands = zeros(1,nBands);

    hBandsRaw    = cell(1,nBands);
    hBandsScaled = hBandsRaw;


    % Start plot of filter response as dB *************************************
    hFigdB  = figure;
    switch ( lower(PlotCommand) )
        case {'plot'}
            plot(FcFilters(1),DesiredGains(1),'k+','LineWidth',2,'MarkerSize',10);
        case {'semilogx'}
            semilogx(FcFilters(1),DesiredGains(1),'k+','LineWidth',2,'MarkerSize',10);
        otherwise
            error(['Unknown plot command: ',PlotCommand]);
    end%switch
    grid on; hold on;
    LegendText = cell(1,nBands+2);
    xlabel('Frequency [Hz]');
    ylabel('Filter Response [dB]');
    LegendText{1} = 'Desired';
    legend ("location", "northeastoutside");


    % Calculate coefficients for each filter ***************************
    for iBand = 1:nBands
        if     ( iBand == 1 ) %then
            [ b, a ] = lowpass( (FcFilters(iBand)+BWidthFilters(iBand)/2)/F0, Qs(1) );
        elseif ( iBand == nBands ) %then
            [ b, a ] = highpass((FcFilters(iBand)-BWidthFilters(iBand)/2)/F0, Qs(2) );
        else
            [ b, a ] = bandpass(FcFilters(iBand) , BWidthFilters(iBand) , Fs ) ;
        end%if
        [hBand,wBand] = freqz(b,a,nFreqs); % Get filter response

        % Convert to filter taps to biquad coefficients
        [SOStaps,gSOStaps] = tf2sos(b,a);

        % If gain is less than unity, filter is "hot", could overflow
        if ( gSOStaps < 1 ) % So scale b's instead
            SOStaps(1,1:3) = SOStaps(1,1:3)*gSOStaps;
            gSOStaps = 1;
        end%if

        % Duplicate filter taps to match number of filters per band
        SOS = zeros(nFiltersPerBand,length(SOStaps));
        gSOS = 1;
        for iFilter = 1:nFiltersPerBand
            SOS(iFilter,:) = SOStaps;
            gSOS = gSOStaps*gSOS;
        end%for

        % Round to nFilterBits bit coefficients
        SOS = round( (2^(nFilterBits-1))*SOS ) / ( 2^(nFilterBits-1) );

        gSOS = peakAmp(iBand)*gSOS/max(abs(hBand));      % Correct gain
        expGain  = ceil( log2(gSOS) );                   % Calculate gain exponent
        fracGain = round( (2^15) * gSOS / (2^expGain) ); % Gain mantissa is always a Q15 number
        gSOS = (fracGain * (2^expGain))/(2^15);          % Actual gain for this filter

        expGain  = ceil( log2(gSOS) );                   % Calculate gain exponent
        fracGain = round( (2^15) * gSOS / (2^expGain) ); % Gain mantissa is always a Q15 number
        gSOS = (fracGain * (2^expGain))/(2^15);          % Actual gain for this filter

        SOSBands{iBand}  = SOS;   % Store filter
        GainBands(iBand) = gSOS;

        [b,a] = sos2tf(SOSBands{iBand},GainBands(iBand)); % Convert from biquad to transfer function coefficients
        [hBand,wBand] = freqz(b,a,nFreqs); % Calculate filter response
        hBandsRaw{iBand} = hBand; % Store unscaled filter response for later use

        % Accumulate total filter response over all bands
        if ( iBand == 1 )
            hTotal = SignBands(1)*hBand;
        else
            hTotal = hTotal + SignBands(iBand)*hBand;
        end%if ( iBand == 1 )

        hBand = hBandsRaw{iBand};
        iLineSpec = mod(iBand-1,length(plotLineSpecs)) + 1;

        LegendText{iBand+1} = ['Filter ',num2str(iBand)];
        figure(hFigdB);
        switch ( lower(PlotCommand) )
            case {'plot'}
                plot((Fs/2)*wBand/pi, MyDB(abs(hBand)), plotLineSpecs{iLineSpec} );
            case {'semilogx'}
                semilogx((Fs/2)*wBand/pi, MyDB(abs(hBand)), plotLineSpecs{iLineSpec} );
            otherwise
                error(['Unknown plot command: ',PlotCommand]);
        end%switch

    end%for iBand = 1:nBands

    % Calculate average and peak filter response
    iBad = find(isnan(hTotal));
    hTotal(iBad) = 0;
    AvgResponse  = mean(abs(hTotal));
    PeakResponse = max(abs(hTotal));

    % Finish with plot of overall filter response
    switch ( lower(PlotCommand) )
        case {'plot'}
            plot((Fs/2)*wBand/pi,MyDB(abs(hTotal)),'--k','LineWidth',2);
        case {'semilogx'}
            semilogx((Fs/2)*wBand/pi,MyDB(abs(hTotal)),'--k','LineWidth',2)
        otherwise
            error(['Unknown plot command: ',PlotCommand]);
    end%switch
    LegendText{end} = 'Total';
    if ( GotMatlab ) %then
        hLegendMag = legend(get(gcf,'CurrentAxes'),LegendText,'Location','Best');
    else
        hLegendMag = legend(get(gcf,'CurrentAxes'),LegendText);
    end%if
    set(hLegendMag,'FontSize',8);
    title(['Filter Design ',num2str(FilterDesignID),', Peak Filter Response: ',num2str(MyDB(PeakResponse)),' dB, over ',num2str(nBands),' bands']);
    axis([-inf,inf,-20,inf]);

    % Plot out desired gains in same color as filter responses
    for iBand = 1:nBands
        iLineSpec = mod(iBand-1,length(plotLineSpecs)) + 1;
        switch ( lower(PlotCommand) )
            case {'plot'}
                plot(FcFilters(iBand),DesiredGains(iBand),[plotLineSpecs{iLineSpec},'+'],'LineWidth',2,'MarkerSize',10);
                plot(FcFilters(iBand),DesiredGains(iBand),[plotLineSpecs{iLineSpec},'o'],'LineWidth',2,'MarkerSize',10);
            case {'semilogx'}
                semilogx(FcFilters(iBand),DesiredGains(iBand),[plotLineSpecs{iLineSpec},'+'],'LineWidth',2,'MarkerSize',10);
                semilogx(FcFilters(iBand),DesiredGains(iBand),[plotLineSpecs{iLineSpec},'o'],'LineWidth',2,'MarkerSize',10);
            otherwise
                error(['Unknown plot command: ',PlotCommand]);
        end%switch
    end%for

    % Display filter design parameters on Command Window for later comparison
    disp(' ');disp(['FILTER DESIGN ',num2str(FilterDesignID),' **************************************']);
    disp(['nFilterBits: ',num2str(nFilterBits),', Fs: ',num2str(Fs),' Hz, # Bands: ',num2str(length(FcFilters))]);
    disp(['Center Freqs [Hz]:  ', sprintf('%10g\t',FcFilters)]);
    disp(['Filter Gains [dB]:  ', sprintf('%10g\t',DesiredGains)]);
    disp(['Filter BWidth [Hz]: ',sprintf('%10g\t',BWidthFilters)]);
    QsChosen = sprintf('%10g',Qs(1));
    for iBand = 2:(length(FcFilters)-1)
        QsChosen = [QsChosen,sprintf('          NA')];
    end%for
    QsChosen = [QsChosen,sprintf('\t%10g',Qs(2))];
    disp(['LowPass/HighPass Qs:', QsChosen]);

    % Are we done?? ***********************************************************
    drawnow; %Flush plotting commands on Octave.
    DoneResponse = questdlg('Filter Design Complete?','Done?','Yes','No','No');
    switch ( DoneResponse )
        case 'Yes'
            bDone = true;
        case 'No'
            bDone = false;
        otherwise
    end%switch

end%while (~bDone)


% Scale all band filters to produce peak of unity for total response **********
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
    hBandsScaled{iBand} = SignBands(iBand)*hBand;

    % Accumulate total filter response over all bands
    if ( iBand == 1 )
        hTotal = SignBands(iBand)*hBand;
    else
        hTotal = hTotal + SignBands(iBand)*hBand;
    end%if ( iBand == 1 )

end%for iBand = 1:nBands
iBad = find(isnan(hTotal));
hTotal(iBad) = 0;
PeakScaledResponse = max(abs(hTotal));
AvgScaledResponse  = mean(abs(hTotal));


% Replot scaled filters for user inspection ***********************************
hFigScaled  = figure;
grid on; hold on;

for iBand = 1:nBands
    hBand = hBandsScaled{iBand};
    iLineSpec = mod(iBand-1,length(plotLineSpecs)) + 1;
    switch ( lower(PlotCommand) )
        case {'plot'}
            plot((Fs/2)*wBand/pi, MyDB(abs(hBand)), plotLineSpecs{iLineSpec} );
        case {'semilogx'}
            semilogx((Fs/2)*wBand/pi, MyDB(abs(hBand)), plotLineSpecs{iLineSpec} );
        otherwise
            error(['Unknown plot command: ',PlotCommand]);
    end%switch

end%for iBand = 1:nBands

switch ( lower(PlotCommand) )
    case {'plot'}
        plot((Fs/2)*wBand/pi,MyDB(abs(hTotal)),'--k','LineWidth',2);
    case {'semilogx'}
        semilogx((Fs/2)*wBand/pi,MyDB(abs(hTotal)),'--k','LineWidth',2)
    otherwise
        error(['Unknown plot command: ',PlotCommand]);
end%switch
LegendText = LegendText(2:end);
if ( GotMatlab ) %then
    hLegendMag = legend(get(gcf,'CurrentAxes'),LegendText,'Location','Best');
else
    hLegendMag = legend(get(gcf,'CurrentAxes'),LegendText);
end%if
xlabel('Frequency [Hz]');
ylabel('Filter Response [dB]');
set(hLegendMag,'FontSize',8);
title(['Filter Design ',num2str(FilterDesignID),', Peak Scaled Filter Response: ',num2str(MyDB(PeakScaledResponse)),' dB, over ',num2str(nBands),' bands']);
axis([-inf,inf,-20,inf]);
drawnow; %Flush plotting commands on Octave


% Save workspace for later use ************************************************
DefaultFileName = ['myFilters',num2str(nBands),'x',num2str(nFiltersPerBand),'_Q',num2str(nFilterBits-1),'.mat'];
[matfile, matpath] =  uiputfile({'*.mat'},'Filter Definition Workspace',DefaultFileName);
if ( matfile ~= 0 ) %then - output workspace
    pathFilterWorkspace = [matpath matfile];
    save(pathFilterWorkspace,...
        'AvgResponse','AvgScaledResponse','BWidthFilters','DesiredGains','FcFilters', ...
        'FilterDesignID','Fs','GainBands','PeakResponse','PeakScaledResponse','Qs', ...
        'SignBands','SOSBands','hBandsRaw','hBandsScaled','hTotal','nBands', ...
        'nFilterBits','nFilterOrder','nFiltersPerBand','nFreqs','wBand');
end%if ( matfile ~= 0 )


% Output header file to initialize filter array *******************************
DefaultFileName = ['myFilters',num2str(nBands),'x',num2str(nFiltersPerBand),'_Q',num2str(nFilterBits-1),'.h'];
[filterfile, filterpath] =  uiputfile({'*.h'}, 'Filter Definition Header', DefaultFileName);

if ( filterfile ~= 0 ) %then - output filter definition
    pathFilterDefn = [filterpath filterfile];
    fidFilterDefn = fopen(pathFilterDefn,'w');

    % Print out start of header
    fprintf(fidFilterDefn,'/*******************************************************************************\n');
    fprintf(fidFilterDefn,'  Audio Equalizer Library Filter Definition\n   (Created by GraphicEqualizerDesign.m)\n\n');
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
    fprintf(fidFilterDefn,'  </pre>\n\n');
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
        %    SignBands sign applied to fracGain value written out so that
        %       Yout = Yband1 + Yband2 + Yband3... + YbandN
        %    rather than
        %       Yout = Yband1 - Yband2 + Yband3... +/- YbandN
        %
        for iFilter = 1:nFiltersPerBand % for each filter...

            % Start new filter
            fprintf(fidFilterDefn,'\t\t{ //Filter %d\n',iFilter);

            % Write out {Fractional Gain, Gain Exponent}, Log2Alpha
            if ( iFilter == nFiltersPerBand ) % If last filter in the cascade, use gain
                fprintf(fidFilterDefn,'\t\t\t{%d,%d},%d,\n',SignBands(iBand)*fracGain,expGain,log2Alpha);
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


% Ask to clean up Matlab workspace or leave it for additional design iterations
ClearWorkspaceResponse = questdlg('Clear workspace at end?','Clear Workspace','Yes','No','No');
switch ( ClearWorkspaceResponse )
    case 'Yes'
        clear all
    case 'No'
    otherwise
end%switch
