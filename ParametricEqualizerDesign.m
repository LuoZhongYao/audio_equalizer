%
% Matlab script to generate parametric equalizer filters
%
% Inputs:
%   Filter bit width, 16 or 32
%   Sample Rate, in Hz
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
r = sqrt( 2 ) / 2 ; % CPB factor
nFilterOrder = 2; % All filters are single biquad IIRs


% Setup defaults for diaglog window *******************************************
if ( ~exist('nFilterBits','var') ) %then
    nFilterBits = 32; % Warning: 16 bit filters show truncation artifacts at low frequencies
end%if
if ( ~exist('Fs','var') ) %then
    Fs = 48000; % Default sample rate, in Hz
end%if
if ( ~exist('FcFilters','var') ) %then
    FcFilters = [32    64    125    250    500    1000    2000    4000    8000    16000]; % Default filter center frequencies
end%if
if ( ~exist('DesiredGains','var') ) %then
    DesiredGains = [7    7    4    1    0    4    1    5    7    4]; % Default filter gains, in dB;
end%if
if ( ~exist('BWidthFilters','var') )%then
    BWidthFilters = r * FcFilters(2:(end-1)) ; % Bandwidth for filters.  Remove low, high - since Q defines bw.
end%if
if ( ~exist('ShelvingQs','var') ) %then
    ShelvingQs = (sqrt(2)/2)*[1,1]; % Default Shelving Q factors
end%if
if ( ~exist('PlotCommand','var') ) %then
    PlotCommand = 'pitch'; % Default: logarithmic frequency axis
end%if

function draw(command, f, varargin)
    switch ( lower(command) )
        case {'plot'}
            plot(f, varargin{:} );
        case {'semilogx'}
            semilogx(f, varargin{:} );
        case {'pitch'}
            p = 69 + 12 * log2(f / 440);
            x = [20,30,40,50,60,80,100,200,300,400,500,600,800,1000,2000,3000,4000,5000,6000,8000,10000,20000]
            y = 69 + 12 * log2(x/440)
            plot(p, varargin{:} );
            set(gca, 'xtick', y)
            set(gca, 'xticklabel', x)
        otherwise
            error(['Unknown plot command: ', command]);
        end%switch
end%function draw

% Iterate filter design until satisfied ***************************************
bDone = false;
hFigdB  = figure;
while (~bDone)

    % Setup configuration dialog window
    DialogTitle     = 'Parametric Equalization Filter Design';
    DialogPrompts   = { 'Filter_Bit_Width:', ...
                        'Sample_Rate_[Hz]:', ...
                        'Center_Frequencies_[Hz]:', ...
                        'Gains_[dB]:', ...
                        'Bandwidths_[Hz]:',...
                        'Shelving_Qs_[Low,High]:', ...
                        'Plot_Command_[plot,semilogx,pitch]:'
                      };
    FcFilterDefaults   = sprintf('%10g\t',FcFilters);
    GainFilterDefaults = sprintf('%10g\t',DesiredGains);
    BandWidthDefaults  = [sprintf('         NA\t'),sprintf('%10g\t',BWidthFilters),'         NA'];
    ShelvingQsDefault   = sprintf('%10g',ShelvingQs(1));
    for iFilter = 2:(length(FcFilters)-1)
        ShelvingQsDefault = [ShelvingQsDefault,sprintf('\t        NA')];
    end%for
    ShelvingQsDefault = [ShelvingQsDefault,sprintf('\t%10g',ShelvingQs(2))];
    nColumns = max( [length(FcFilterDefaults), ...
                     length(GainFilterDefaults),  ...
                     length(BandWidthDefaults), ...
                     length(ShelvingQsDefault)] );
    DialogDefaults = { num2str(nFilterBits), ...
                       num2str(Fs), ...
                       FcFilterDefaults,...
                       GainFilterDefaults,...
                       BandWidthDefaults,...
                       ShelvingQsDefault,...
                       PlotCommand ...
                     };
    options.WindowStyle='normal';
    %             Row, Columns
    rowcolumns = [  1, length('Filter Bit Width:'); ... %nFilterbits
                    1, length('Sample_Rate_[Hz]:'); ... %Fs
                    1, nColumns; ... % Center frequencies
                    1, nColumns; ... % Gains
                    1, nColumns; ... % Bandwidths
                    1, nColumns; ... % Shelving Qs
                    1, length('Plot Command [plot,semilogx]:') % Plot command: linear or log frequency axis
                 ];

    % Get user choices
    DialogAnswers = inputdlg(DialogPrompts,DialogTitle,rowcolumns,DialogDefaults,options);
    assert( length(DialogAnswers) == 7,'Dialog terminated! Dialog Length: %d',length(DialogAnswers));
    nFilterBits   = str2num(DialogAnswers{1});
    Fs            = str2num(DialogAnswers{2});
    FcFilters     = str2num(DialogAnswers{3});
    DesiredGains   = str2num(DialogAnswers{4});
    BWidthFilters = str2num(strrep(DialogAnswers{5},'NA',''));
    ShelvingQs    = str2num(strrep(DialogAnswers{6},'NA',''));
    PlotCommand   = DialogAnswers{7};

    % Validate user input
    assert( Fs > 0,'Sampling is not positive: %g',Fs)
    nFiltersPerBand = length(FcFilters); % Single biquad for each band

    nFcFilters = length(FcFilters);
    nGainFilters = length(DesiredGains);
    nBandWidths = length(BWidthFilters);
    if ( nFcFilters ~= nGainFilters ) %then
        error(['Number of frequencies (',num2str(nFcFilters),') doesn''t match gains specified (',num2str(nGainFilters),')!']);
    end%if
    if ( (nFcFilters-2) ~= nBandWidths ) %then
        error(['Number of interior frequencies (',num2str(nFcFilters-2),') doesn''t match bandwidths specified (',num2str(nBandWidths),')!']);
    end%if

    nBands = 1;
    nFiltersPerBand = nFcFilters;

    switch ( lower(PlotCommand) )
        case 'pitch'
            %Do Nothing
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
    tempSOS = tf2sos( ones(1,nFilterOrder+1), ones(1,nFilterOrder+1) );
    tempSOS = zeros(size(tempSOS));

    SOSFilters = cell(nFcFilters,1);
    for iFilter = 1:nFiltersPerBand
        SOSFilters{iFilter} = tempSOS;
    end%for
    GainFilters = zeros(1,nFiltersPerBand);

    hFiltersRaw    = cell(1,nFiltersPerBand);


    % Start plot of filter response as dB *************************************
    clf(hFigdB, 'reset')
    draw(PlotCommand, FcFilters(1),DesiredGains(1),'k+','LineWidth',2,'MarkerSize',10);
    % switch ( lower(PlotCommand) )
    %     case {'plot'}
    %         plot(FcFilters(1),DesiredGains(1),'k+','LineWidth',2,'MarkerSize',10);
    %     case {'semilogx'}
    %         semilogx(FcFilters(1),DesiredGains(1),'k+','LineWidth',2,'MarkerSize',10);
    %     otherwise
    %         error(['Unknown plot command: ',PlotCommand]);
    % end%switch
    grid on; hold on;
    LegendText = cell(1,nFiltersPerBand+2);
    xlabel('Frequency [Hz]');
    ylabel('Filter Response [dB]');
    title(['Parametric Equalizer Filters, ',num2str(nFiltersPerBand),' Filters']);
    LegendText{1} = 'Desired';
    legend ("location", "northeastoutside");


    % Calculate coefficients for each filter ***************************
    for iFilter = 1:nFiltersPerBand

        % Calculate IIR biquad coefficients
        if     ( iFilter == 1 )      %then - Bass Shelf
            [ b, a ] = shelving( 'Bass_Shelf' , DesiredGains( 1 ) , FcFilters(1) , ShelvingQs(1) , Fs ) ;
        elseif ( iFilter == nFiltersPerBand ) %then - Treble Shelf
            [ b, a ] = shelving( 'Treble_Shelf' , DesiredGains(nFiltersPerBand) , FcFilters(nFiltersPerBand) , ShelvingQs(2) , Fs ) ;
        else                       % then - interior filter
            [ b, a ] = parmeq( DesiredGains(iFilter) , FcFilters(iFilter) , BWidthFilters(iFilter-1) , Fs ) ;
        end%if

        % Convert to filter taps to biquad coefficients
        [SOS,gSOS] = tf2sos(b,a);

        if ( gSOS < 1 ) % then scale b's instead
            SOS(1,1:3) = SOS(1,1:3)*gSOS;
            gSOS = 1;
        end%if

        % Round to 16 bit coefficients, with log2Alpha = 1
        SOS = round( (2^(nFilterBits-1))*SOS/2 ) / (2^(nFilterBits-2) );

        expGain  = ceil( log2(gSOS) );                   % Calculate gain exponent
        fracGain = round( (2^15) * gSOS / (2^expGain) ); % Gain mantissa is always a Q15 number
        gSOS = (fracGain * (2^expGain))/(2^15);          % Actual gain for this filter

        SOSFilters{iFilter} = SOS;
        GainFilters(iFilter) = gSOS;

        [b,a] = sos2tf(SOSFilters{iFilter},GainFilters(iFilter)); % Convert from biquad to transfer function coefficients
        [hFilter,wFilter] = freqz(b,a,nFreqs); % Calculate filter response
        hFiltersRaw{iFilter} = hFilter; % Store unscaled filter response for later use

        % Accumulate total filter response over all bands
        if ( iFilter == 1 )
            hTotal = hFilter;
        else
            hTotal = hTotal.*hFilter;
        end%if ( iFilter == 1 )

        hFilter = hFiltersRaw{iFilter};
        iLineSpec = mod(iFilter-1,length(plotLineSpecs)) + 1;

        LegendText{iFilter+1} = ['Filter ',num2str(iFilter)];
        figure(hFigdB);
        draw(PlotCommand, (Fs/2)*wFilter/pi, MyDB(abs(hFilter)), plotLineSpecs{iLineSpec} );
        % switch ( lower(PlotCommand) )
        %     case {'plot'}
        %         plot((Fs/2)*wFilter/pi, MyDB(abs(hFilter)), plotLineSpecs{iLineSpec} );
        %     case {'semilogx'}
        %         semilogx((Fs/2)*wFilter/pi, MyDB(abs(hFilter)), plotLineSpecs{iLineSpec} );
        %     otherwise
        %         error(['Unknown plot command: ',PlotCommand]);
        % end%switch

    end%for iFilter = 1:nFiltersPerBand

    % Calculate average and peak filter response
    AvgResponse  = mean(abs(hTotal));
    PeakResponse = max(abs(hTotal));

    % Finish with plot of overall filter response
    draw(PlotCommand, (Fs/2)*wFilter/pi, MyDB(abs(hTotal)), '--k', 'LineWidth', 2)
    % switch ( lower(PlotCommand) )
    %     case {'plot'}
    %         plot((Fs/2)*wFilter/pi,MyDB(abs(hTotal)),'--k','LineWidth',2);
    %     case {'semilogx'}
    %         semilogx((Fs/2)*wFilter/pi,MyDB(abs(hTotal)),'--k','LineWidth',2)
    %     otherwise
    %         error(['Unknown plot command: ',PlotCommand]);
    % end%switch
    LegendText{end} = 'Total';
    if ( GotMatlab ) %then
        hLegendMag = legend(get(gcf,'CurrentAxes'),LegendText,'Location','Best');
    else
        hLegendMag = legend(get(gcf,'CurrentAxes'),LegendText);
    end%if
    set(hLegendMag,'FontSize',8);
    title(['Peak Filter Response: ',num2str(MyDB(PeakResponse)),' dB, over ',num2str(nFiltersPerBand),' Filters']);

    % Plot out desired gains in same color as filter responses
    for iFilter = 1:nFiltersPerBand
        iLineSpec = mod(iFilter-1,length(plotLineSpecs)) + 1;
        draw(PlotCommand, FcFilters(iFilter),DesiredGains(iFilter),[plotLineSpecs{iLineSpec},'+'],'LineWidth',2,'MarkerSize',10);
        draw(PlotCommand, FcFilters(iFilter),DesiredGains(iFilter),[plotLineSpecs{iLineSpec},'o'],'LineWidth',2,'MarkerSize',10);
        % switch ( lower(PlotCommand) )
        %     case {'plot'}
        %         plot(FcFilters(iFilter),DesiredGains(iFilter),[plotLineSpecs{iLineSpec},'+'],'LineWidth',2,'MarkerSize',10);
        %         plot(FcFilters(iFilter),DesiredGains(iFilter),[plotLineSpecs{iLineSpec},'o'],'LineWidth',2,'MarkerSize',10);
        %     case {'semilogx'}
        %         semilogx(FcFilters(iFilter),DesiredGains(iFilter),[plotLineSpecs{iLineSpec},'+'],'LineWidth',2,'MarkerSize',10);
        %         semilogx(FcFilters(iFilter),DesiredGains(iFilter),[plotLineSpecs{iLineSpec},'o'],'LineWidth',2,'MarkerSize',10);
        %     otherwise
        %         error(['Unknown plot command: ',PlotCommand]);
        % end%switch
    end%for

    % Display filter design parameters on Command Window for later comparison
    disp(' ');disp('FILTER DESIGN **************************************');
    disp(['nFilterBits: ',num2str(nFilterBits),', Fs: ',num2str(Fs),' Hz, # Bands: ',num2str(length(FcFilters))]);
    disp(['Center Freqs [Hz]:  ', sprintf('%10g\t',FcFilters)]);
    disp(['Filter Gains [dB]:  ', sprintf('%10g\t',DesiredGains)]);
    disp(['Filter BWidth [Hz]: ', sprintf('        NA\t'),sprintf('%10g\t',BWidthFilters),'        NA']);
    ShelvingQsChosen = sprintf('%10g',ShelvingQs(1));
    for iFilter = 2:(length(FcFilters)-1)
        ShelvingQsChosen = [ShelvingQsChosen,sprintf('          NA')];
    end%for
    ShelvingQsChosen = [ShelvingQsChosen,sprintf('\t%10g',ShelvingQs(2))];
    disp(['Shelving Qs:        ', ShelvingQsChosen]);

    % Are we done?? ***********************************************************
    drawnow; %Flush plotting commands on Octave.
    DoneResponse = questdlg('Filter Design Complete?','Done?','Yes','No','No');
    switch ( DoneResponse )
        case 'Yes'
            bDone = true;
        case 'No'
        otherwise
    end%switch

end%while (~bDone)


% Save workspace for later use ************************************************
DefaultFileName = ['ParametricFilters',num2str(nBands),'x',num2str(nFiltersPerBand),'_Q',num2str(nFilterBits-1),'.mat'];
[matfile, matpath] =  uiputfile({'*.mat'},'Filter Definition Workspace',DefaultFileName);
if ( matfile ~= 0 ) %then - output workspace
    pathFilterWorkspace = [matpath filesep() matfile];
    save(pathFilterWorkspace,...
        'AvgResponse','BWidthFilters','DesiredGains','FcFilters','Fs', ...
        'GainFilters','PeakResponse','ShelvingQs','SOSFilters','hFiltersRaw', ...
        'hTotal','nBands','nFilterBits','nFilterOrder','nFiltersPerBand','nFreqs','wFilter');
end%if ( matfile ~= 0 )


% Output header file to initialize filter array *******************************
DefaultFileName = ['ParametricFilters',num2str(nBands),'x',num2str(nFiltersPerBand),'_Q',num2str(nFilterBits-1),'.h'];
[filterfile, filterpath] =  uiputfile({'*.h'}, 'Filter Definition Header', DefaultFileName);

if ( filterfile ~= 0 ) %then - output filter definition
    pathFilterDefn = [filterpath filesep() filterfile];
    fidFilterDefn = fopen(pathFilterDefn,'w');

    % Print out start of header
    fprintf(fidFilterDefn,'/*******************************************************************************\n');
    fprintf(fidFilterDefn,'  Audio Equalizer Library Filter Definition\n   (Created by ParametricEqualizerDesign.m)\n\n');
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
    fprintf(fidFilterDefn,'</pre>\n\n');
    fprintf(fidFilterDefn,'*******************************************************************************/\n');

    % Write out EQUALIZER_FILTER or EQUALIZER_FILTER_32 structure initialization to define filters
    if ( nFilterBits == 16 )
        fprintf(fidFilterDefn,'EQUALIZER_FILTER ParametricFilters[] =\n');
    else
        fprintf(fidFilterDefn,'EQUALIZER_FILTER_32 ParametricFilters[] =\n');
    end%if ( nFilterBits == 16 )
    fprintf(fidFilterDefn,'\t{\n');

    % Write out filters
    iBand = 1; % Only one band

    % Start definitions for the filters
    fprintf(fidFilterDefn,'\t  //Band %d\n',iBand);

    log2Alpha = 1; % Always one

    for iFilter = 1:nFiltersPerBand % for each filter...
        % Start new filter
        fprintf(fidFilterDefn,'\t\t{ //Filter %d\n',iFilter);

        % Derive gain variables
        gFilter = GainFilters(iFilter);
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

        % Write out {Fractional Gain, Gain Exponent}, Log2Alpha
        fprintf(fidFilterDefn,'\t\t\t{%d,%d},%d,\n',fracGain,expGain,log2Alpha);

        SOS = SOSFilters{iFilter}/(2^log2Alpha); % Scale filter coefficients
        SOS = round(SOS*2^(nFilterBits-1));  % Convert coefficients to fixed point integers

        % Write out filter taps: {b0,b1,b2}, then {a1,a2} (a1 assumed to always be 1)
        fprintf(fidFilterDefn,'\t\t\t{%d,%d,%d},\n',SOS(1,1),SOS(1,2),SOS(1,3));
        fprintf(fidFilterDefn,'\t\t\t{%d,%d},\n',SOS(1,5),SOS(1,6));

        % Write out initial values for Z's
        fprintf(fidFilterDefn,'\t\t\t{0L,0L}\n');

        % End filter definition
        fprintf(fidFilterDefn,'\t\t}, //end Filter %d\n',iFilter);

    end%for iFilter = 1:nFiltersPerBand

    % End definitions for band filters
    fprintf(fidFilterDefn,'\t  //end Band %d\n\n',iBand);

    fprintf(fidFilterDefn,'\t};\n'); % Last line of initialization.

    fclose(fidFilterDefn); % Close header file, we're done.

end%if ( filename ~= 0 )

ClearWorkspaceResponse = questdlg('Clear workspace at end?','Clear Workspace','Yes','No','No');
switch ( ClearWorkspaceResponse )
    case 'Yes'
        clear all
    case 'No'
    otherwise
end%switch

% Based on peq.m
%   Designed by Mark Ledwich, 91.211, Spring 2011
%   Modified by Stu Smith, 12 May 2011
%
% peg.m is from SRT Matlab Audio Toolkit, which is part of the course
%   91.211 Computer Science for SRT Applications
% taught by Stuart Smith (See http://www.cs.uml.edu/~stu/cs211/syllabus.html)
%
