if isempty(ver('lte')) % Check for LST install
    error('zynqRadioLTEMIMOTransmitReceive:NoLST', ...
        'Please install LTE Toolbox to run this example.');
elseif ~license('test', 'LTE_Toolbox') % Check that a valid license is present
    error('zynqRadioLTEMIMOTransmitReceive:NoLST', ...
        'A valid license for LTE Toolbox is required to run this example.');
end

%%
% The example configures all scopes and figures displayed throughout the example. 

% Setup handle for image plot
if ~exist('imFig', 'var') || ~ishandle(imFig)
    imFig = figure;
    imFig.NumberTitle = 'off';
    imFig.Name = 'Image Plot';
    imFig.Visible = 'off';
else   
    clf(imFig); % Clear figure
    imFig.Visible = 'off';
end

% Setup handle for channel estimate plots
if ~exist('hhest', 'var') || ~ishandle(hhest)
    hhest = figure('Visible','Off');
    hhest.NumberTitle = 'off';
    hhest.Name = 'Channel Estimate';
else
    clf(hhest); % Clear figure
    hhest.Visible = 'off';
end

% Setup Spectrum viewer
spectrumScope = dsp.SpectrumAnalyzer( ...
    'SpectrumType',    'Power density', ...
    'SpectralAverages', 10, ...
    'YLimits',         [-130 -40], ...
    'Title',           'Received Baseband LTE Signal Spectrum', ...
    'YLabel',          'Power spectral density');

% Setup the constellation diagram viewer for equalized PDSCH symbols
constellation = comm.ConstellationDiagram('Title','Equalized PDSCH Symbols',...
                                'ShowReferenceConstellation',false);
                            
%%

%  Initialize SDR device
txsim = struct; % Create empty structure for transmitter
txsim.SDRDeviceName = 'AD936x'; % Set SDR Device
% To update the example for FMCOMMS5, set |txsim.SDRDeviceName| to
% |'FMCOMMS5'|.
radio = sdrdev(txsim.SDRDeviceName); % Create SDR device object

%% Run Example
% You can run this example by typing the script name of this example in the
% MATLAB command window:
%
% 	|zynqRadioLTEMIMOTransmitReceiveAD9361AD9364ML|
%
% The following sections explain the design and architecture of this example, 
% and what you can expect to see as the code is executed.

%% 
% *Set Up SDR Transmitter*
%
% The |txsim| structure controls the properties of the SDR transmitter
% System object. By default, the example uses an AD936x device with two
% channels. If you are using an FMCOMMS4 RF card with your device, you have
% only one channel. Set |txsim.NTxAnts| to 1.
%
% If you have set |txsim.SDRDeviceName| to |'FMCOMMS5'| in the previous
% section, the example automatically chooses four channels.

txsim.RC = 'R.3';       % Base RMC configuration, 10 MHz bandwidth
txsim.NCellID = 88;     % Cell identity
txsim.NFrame = 700;     % Initial frame number
txsim.TotFrames = 1;    % Number of frames to generate
txsim.DesiredCenterFrequency = 2.45e9; % Center frequency in Hz
switch txsim.SDRDeviceName
    case 'AD936x'
   %    txsim.NTxAnts = 2;      % Number of transmit antennas
        txsim.NTxAnts = 1;
    case 'FMCOMMS5'
        txsim.NTxAnts = 4;
    otherwise
        error('Unknown SDR device: %s',txsim.SDRDeviceName);
end

% txsim.NTxAnts = 1;

%%
% To visualize the benefit of using multi-channel transmission
% TX gain parameter: 
txsim.Gain = -50;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MsgLength=1e6;
% Create RMC
rmc = lteRMCDL(txsim.RC);
%rmc.PDSCH.Modulation='QPSK';

% Calculate the required number of LTE frames based on the size of the
% image data
trBlkSize = rmc.PDSCH.TrBlkSizes;
txsim.TotFrames = ceil(MsgLength/sum(trBlkSize(:)));

% Customize RMC parameters
rmc.NCellID = txsim.NCellID;
rmc.NFrame = txsim.NFrame;
rmc.TotSubframes = txsim.TotFrames*10; % 10 subframes per frame
rmc.CellRefP = txsim.NTxAnts; % Configure number of cell reference ports
rmc.PDSCH.RVSeq = 0;

% Fill subframe 5 with dummy data
rmc.OCNGPDSCHEnable = 'On';
rmc.OCNGPDCCHEnable = 'On';

% If transmitting over two channels enable transmit diversity
if rmc.CellRefP >= 2
    rmc.PDSCH.TxScheme = 'TxDiversity';
    rmc.PDSCH.NLayers = txsim.NTxAnts;
    rmc.OCNGPDSCH.TxScheme = 'TxDiversity';
end

fprintf('\nGenerating LTE transmit waveform:\n')
fprintf('  Packing image data into %d frame(s).\n\n', txsim.TotFrames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SDR Transmitter

sdrTransmitter = sdrtx(txsim.SDRDeviceName);
sdrTransmitter.BasebandSampleRate = 15360000; % 15.36 Msps for default RMC (R.7) 
sdrTransmitter.CenterFrequency = txsim.DesiredCenterFrequency;
sdrTransmitter.Gain = txsim.Gain;
sdrTransmitter.ShowAdvancedProperties = true;
sdrTransmitter.BypassUserLogic = true;
sdrTransmitter.IPAddress='137.194.172.32';
sdrTransmitter.ChannelMapping =1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Set-up Bit Rx

rxsim = struct;
rxsim.RadioFrontEndSampleRate = sdrTransmitter.BasebandSampleRate; % Configure for same sample rate
                                                       % as transmitter
rxsim.RadioCenterFrequency = txsim.DesiredCenterFrequency;
rxsim.NRxAnts = txsim.NTxAnts;
rxsim.FramesPerCapture = txsim.TotFrames+1; % Number of LTE frames to capture.
                                          % Capture 1 more LTE frame than transmitted to  
                                          % allow for timing offset wraparound...
rxsim.numCaptures = 1; % Number of captures

% Derived parameters
samplesPerFrame = 10e-3*rxsim.RadioFrontEndSampleRate; % LTE frames period is 10 ms
captureTime = rxsim.FramesPerCapture * 10e-3; % LTE frames period is 10 ms

%%
% Create an SDR receiver System object with the specied properties for the
% device used for the image transmission.

rxsim.SDRDeviceName = txsim.SDRDeviceName;
sdrReceiver = sdrrx(rxsim.SDRDeviceName);
sdrReceiver.BasebandSampleRate = rxsim.RadioFrontEndSampleRate;
sdrReceiver.CenterFrequency = rxsim.RadioCenterFrequency;
sdrReceiver.OutputDataType = 'double';
sdrReceiver.ShowAdvancedProperties = true;
sdrReceiver.GainSource='Manual';
sdrReceiver.Gain=0;
sdrReceiver.IPAddress='137.194.172.32';
sdrReceiver.BypassUserLogic = true;
% Configure RX channel map
sdrReceiver.ChannelMapping = 1;

spectrumScope.SampleRate = rxsim.RadioFrontEndSampleRate;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%SNR Loop
% Calculate the SNR values, for those obtained experimentally.
SNRdB=[5.2804 5.3667 5.697 6.2024 7.1965 7.9949 8.6001 10.8099 12.1878 13.8013];
berEst = zeros(size(SNRdB));
modSize= 16; % Modulation order 
k = log2(modSize);
fo=15.36e6;
B=10e6;

for l = 1:length(SNRdB)
EbNo=SNRdB(l) -10*log10(fo)+10*log10(B);
% Reset the error and bit counters
numErrs = 0;
numBits = 0;
txsim.Gain = txsim.Gain +2;
sdrTransmitter.Gain=txsim.Gain;

while numErrs < 1000 && numBits < 20e6
    
    trData=randi([0 1],MsgLength,1);
    [eNodeBOutput,txGrid,rmc] = lteRMCDLTool(rmc,trData);
    % Scale the signal for better power output.
    powerScaleFactor = 0.8;
    if txsim.NTxAnts >= 2
         for i =1:txsim.NTxAnts
             eNodeBOutput(:,i) = eNodeBOutput(:,i).*(1/max(abs(eNodeBOutput(:,i)))*powerScaleFactor);
         end
    else
         eNodeBOutput = eNodeBOutput.*(1/max(abs(eNodeBOutput))*powerScaleFactor);
    end
    % Cast the transmit signal to int16 --- 
    % this is the native format for the SDR hardware. 
    eNodeBOutput = int16(eNodeBOutput*2^15);
    transmitRepeat(sdrTransmitter,eNodeBOutput);

    enb.PDSCH = rmc.PDSCH;
    enb.DuplexMode = 'FDD';
    enb.CyclicPrefix = 'Normal';
    enb.CellRefP = 4; 

    %%
    % The sampling rate of the signal controls the captured bandwidth. The
    % number of RBs captured is obtained from a lookup table using
    % the chosen sampling rate.

    % Bandwidth: {1.4 MHz, 3 MHz, 5 MHz, 10 MHz, 20 MHz}
    SampleRateLUT = [1.92 3.84 7.68 15.36 30.72]*1e6;
    NDLRBLUT = [6 15 25 50 100];
    enb.NDLRB = NDLRBLUT(SampleRateLUT==rxsim.RadioFrontEndSampleRate);
    if isempty(enb.NDLRB)
        error('Sampling rate not supported. Supported rates are %s.',...
               '1.92 MHz, 3.84 MHz, 7.68 MHz, 15.36 MHz, 30.72 MHz');
    end
    fprintf('\nSDR hardware sampling rate configured to capture %d LTE RBs.\n',enb.NDLRB);

    %%
    % Channel estimation is configured to be performed using cell-specific
    % reference signals. A 9-by-9 averaging window is used to minimize the
    % effect of noise.
    % Channel estimation configuration structure
    cec.PilotAverage = 'UserDefined';  % Type of pilot symbol averaging
    cec.FreqWindow = 9;                % Frequency window size in REs
    cec.TimeWindow = 9;                % Time window size in REs
    cec.InterpType = 'Cubic';          % 2D interpolation type
    cec.InterpWindow = 'Centered';     % Interpolation window type
    cec.InterpWinSize = 3;             % Interpolation window size

   %%
   % *Capture and Process Signal*

   enbDefault = enb;
   rxsim.numCaptures = rxsim.numCaptures+1;
   while rxsim.numCaptures
       % Set default LTE parameters
        enb = enbDefault;
        rxWaveform = capture(sdrReceiver, captureTime, 'Seconds');
    
        if rxsim.NRxAnts > 1
            spectrumScope.ShowLegend = true; % Turn on legend for spectrum analyzer
            spectrumScope.ChannelNames = cellfun(@(x) ['SDR Channel ' num2str(x)], num2cell(1:4), 'UniformOutput', false);
        end
    
        % Show power spectral density of captured burst
        %spectrumScope.SampleRate = rxsim.RadioFrontEndSampleRate;
        spectrumScope(rxWaveform);
    
        % Perform frequency offset correction for known cell ID
        frequencyOffset = lteFrequencyOffset(enb,rxWaveform);
        rxWaveform = lteFrequencyCorrect(enb,rxWaveform,frequencyOffset);
        fprintf('\nCorrected a frequency offset of %i Hz.\n',frequencyOffset)
    
        % Perform the blind cell search to obtain cell identity and timing offset
        %   Use 'PostFFT' SSS detection method to improve speed
        cellSearch.SSSDetection = 'PostFFT'; cellSearch.MaxCellCount = 1;
        [NCellID,frameOffset] = lteCellSearch(enb,rxWaveform,cellSearch);
        fprintf('Detected a cell identity of %i.\n', NCellID);
        enb.NCellID = NCellID; % From lteCellSearch
    
        % Sync the captured samples to the start of an LTE frame, and trim off
        % any samples that are part of an incomplete frame.
        rxWaveform = rxWaveform(frameOffset+1:end,:);
        tailSamples = mod(length(rxWaveform),samplesPerFrame);
        rxWaveform = rxWaveform(1:end-tailSamples,:);
        enb.NSubframe = 0;
        fprintf('Corrected a timing offset of %i samples.\n',frameOffset)
    
        % OFDM demodulation
        rxGrid = lteOFDMDemodulate(enb,rxWaveform);
    
        % Perform channel estimation for 4 CellRefP as currently we do not
        % know the CellRefP for the eNodeB.
        [hest,nest] = lteDLChannelEstimate(enb,cec,rxGrid);
    
        sfDims = lteResourceGridSize(enb);
        Lsf = sfDims(2); % OFDM symbols per subframe
        LFrame = 10*Lsf; % OFDM symbols per frame
        numFullFrames = length(rxWaveform)/samplesPerFrame;
    
        rxDataFrame = zeros(sum(enb.PDSCH.TrBlkSizes(:)),numFullFrames);
        recFrames = zeros(numFullFrames,1);
        rxSymbols = []; txSymbols = [];
    
        % For each frame decode the MIB, PDSCH and DL-SCH
        for frame = 0:(numFullFrames-1)
            fprintf('\nPerforming DL-SCH Decode for frame %i of %i in burst:\n', ...
                frame+1,numFullFrames)
        
            % Extract subframe #0 from each frame of the received resource grid
            % and channel estimate.
            enb.NSubframe = 0;
            rxsf = rxGrid(:,frame*LFrame+(1:Lsf),:);
            hestsf = hest(:,frame*LFrame+(1:Lsf),:,:);
               
            % PBCH demodulation. Extract resource elements (REs)
            % corresponding to the PBCH from the received grid and channel
            % estimate grid for demodulation.
            enb.CellRefP = 4;
            pbchIndices = ltePBCHIndices(enb); 
            [pbchRx,pbchHest] = lteExtractResources(pbchIndices,rxsf,hestsf);
            [~,~,nfmod4,mib,CellRefP] = ltePBCHDecode(enb,pbchRx,pbchHest,nest);
        
            % If PBCH decoding successful CellRefP~=0 then update info
            if ~CellRefP
                fprintf('  No PBCH detected for frame.\n');
                continue;
            end
            enb.CellRefP = CellRefP; % From ltePBCHDecode
            % Decode the MIB to get current frame number
            enb = lteMIB(mib,enb);
            % Incorporate the nfmod4 value output from the function
            % ltePBCHDecode, as the NFrame value established from the MIB
            % is the system frame number modulo 4.
            enb.NFrame = enb.NFrame+nfmod4;        
            % The eNodeB transmission bandwidth may be greater than the
            % captured bandwidth, so limit the bandwidth for processing
            enb.NDLRB = min(enbDefault.NDLRB,enb.NDLRB);
            % Store received frame number
            recFrames(frame+1) = enb.NFrame;  
            % Process subframes within frame (ignoring subframe 5)
            for sf = 0:9
                if sf~=5 % Ignore subframe 5
                    % Extract subframe
                    enb.NSubframe = sf;
                    rxsf = rxGrid(:,frame*LFrame+sf*Lsf+(1:Lsf),:);
                    % Perform channel estimation with the correct number of CellRefP
                    [hestsf,nestsf] = lteDLChannelEstimate(enb,cec,rxsf);
                    % PCFICH demodulation. Extract REs corresponding to the PCFICH
                    % from the received grid and channel estimate for demodulation.
                    pcfichIndices = ltePCFICHIndices(enb);
                    [pcfichRx,pcfichHest] = lteExtractResources(pcfichIndices,rxsf,hestsf);
                    [cfiBits,recsym] = ltePCFICHDecode(enb,pcfichRx,pcfichHest,nestsf);
                    % CFI decoding
                    enb.CFI = lteCFIDecode(cfiBits);
                    % Get PDSCH indices
                    [pdschIndices,pdschIndicesInfo] = ltePDSCHIndices(enb, enb.PDSCH, enb.PDSCH.PRBSet); 
                    [pdschRx, pdschHest] = lteExtractResources(pdschIndices, rxsf, hestsf);
                    % Perform deprecoding, layer demapping, demodulation and
                    % descrambling on the received data using the estimate of
                    % the channel
                    [rxEncodedBits, rxEncodedSymb] = ltePDSCHDecode(enb,enb.PDSCH,pdschRx,...
                                                   pdschHest,nestsf);
                    % Append decoded symbol to stream
                    rxSymbols = [rxSymbols; rxEncodedSymb{:}]; %#ok<AGROW>
                    % Transport block sizes
                    outLen = enb.PDSCH.TrBlkSizes(enb.NSubframe+1);
                    % Decode DownLink Shared Channel (DL-SCH)
                    [decbits{sf+1}, blkcrc(sf+1)] = lteDLSCHDecode(enb,enb.PDSCH,...
                                                   outLen, rxEncodedBits);  %#ok<SAGROW>
                    % Recode transmitted PDSCH symbols for EVM calculation                            
                    %   Encode transmitted DLSCH 
                    txRecode = lteDLSCH(enb,enb.PDSCH,pdschIndicesInfo.G,decbits{sf+1});
                    %   Modulate transmitted PDSCH
                    txRemod = ltePDSCH(enb, enb.PDSCH, txRecode);
                    %   Decode transmitted PDSCH
                    [~,refSymbols] = ltePDSCHDecode(enb, enb.PDSCH, txRemod);
                    %   Add encoded symbol to stream
                    txSymbols = [txSymbols; refSymbols{:}]; %#ok<AGROW>
                    %release(constellation); % Release previous constellation plot
                    %constellation(rxEncodedSymb{:}); % Plot current constellation
                    pause(0); % Allow constellation to repaint
                end
            end
        
            % Reassemble decoded bits
            rxdata = [];
            for i = 1:length(decbits)
                if i~=6 % Ignore subframe 5
                    rxdata = [rxdata; decbits{i}{:}]; %#ok<AGROW>
                end
            end
        
            % Store data from receive frame
            rxDataFrame(:,frame+1) = rxdata;

            % Plot channel estimate between CellRefP 0 and the receive antennae
            focalFrameIdx = frame*LFrame+(1:LFrame);
        end
        rxsim.numCaptures = rxsim.numCaptures-1;
    end

    [~,frameIdx] = min(recFrames);
    decodedRxDataStream = zeros(length(rxDataFrame(:)),1);
    frameLen = size(rxDataFrame,1);
    for n=1:numFullFrames
        currFrame = mod(frameIdx-1,numFullFrames)+1; % Get current frame index 
        decodedRxDataStream((n-1)*frameLen+1:n*frameLen) = rxDataFrame(:,currFrame);
        frameIdx = frameIdx+1; % Increment frame index
    end
    [nErrors,BitErr] = biterr(decodedRxDataStream(1:length(trData)), trData);

    n64QAMSymbols=round(length(trData)/6);
    numErrs = numErrs + nErrors;
    numBits = numBits + round(length(trData));

end
berEst(l) = numErrs/numBits;
end
% Release both the transmitter and receiver objects once reception is complete
release(sdrTransmitter);
release(sdrReceiver);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the curves obtained 

berTheory = berawgn(EbNoVec,'qam',16);
figure()
semilogy(EbNoVec,berEst)
%hold on
%semilogy(EbNoVec,berTheory)
grid
legend('Estimated BER','Theoretical BER')
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
title('BER vs SNR 16QAM MIMO')


% Perform bit error rate (BER) calculation
%bitErrorRate = comm.ErrorRate;
%err = bitErrorRate(decodedRxDataStream(1:length(trData)), trData);
%fprintf('  Bit Error Rate (BER) = %0.5f.\n', err(1));
%fprintf('  Number of bit errors = %d.\n', err(2));
%fprintf('  Number of transmitted bits = %d.\n',length(trData));


%% Things to Try
% By default, the example uses two antennas for transmission and
% reception of the LTE waveform. Depending on your hardware, you can modify
% the SDR transmitter and receiver to use a single or four antennas. To 
% observe the difference in the EVM and BER after signal reception and 
% processing, you can decrease the transmitter gain. Check for any errors 
% in the displayed received image.

%% Troubleshooting the Example
%
% For more information on troubleshooting SDR hardware and the Communications
% Toolbox Support Package for Xilinx Zynq-Based Radio, see
% <docid:xilinxzynqbasedradio_ug#bu1kohy-57 Common Problems and Fixes>.

%% Selected Bibliography
% # 3GPP TS 36.191. "User Equipment (UE) radio transmission and reception."
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network; Evolved Universal Terrestrial Radio Access (E-UTRA).

displayEndOfDemoMessage(mfilename)