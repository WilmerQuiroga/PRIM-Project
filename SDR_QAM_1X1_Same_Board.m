%% Receive Tone Signal Using Analog Devices AD9361/AD9364
%
% This example shows how to use the Xilinx(R) Zynq-Based Radio Support
% Package and Communications Toolbox(TM) software to perform a
% simple loopback of a complex sinusoid signal at RF Using Analog
% Devices(TM) AD9361/AD9364. A Direct Digital
% Synthesizer (DDS) in the FPGA generates a complex sinusoid and transmits
% it using the RF card. The transmitted signal is then received by the RF
% card and the downsampled baseband signal is visualized in MATLAB(R). This
% simple example confirms that the SDR hardware is setup correctly and
% shows how to capture RF data from SDR hardware using MATLAB.

% Copyright 2014-2020 The MathWorks, Inc.

%% Configure SDR Hardware
% If you are running this example after performing the setup of
% hardware using Support Package Installer then you can skip this section.
% Otherwise, follow <docid:xilinxzynqbasedradio_ug#bu1kohy-14  Guided Host-Radio 
% Hardware Setup> to configure your host computer to work with the SDR hardware.
% Connect an SMA loopback cable with attenuation between TX1A and RX1A (for
% FMCOMMS2 or FMCOMMS3), between TXA and RXA (for FMCOMMS4), between
% TX1A_A and RX1A_A (for FMCOMMS5), or attach appropriate antenna suitable
% for 2.4 GHz band.

%% Running the Example
% This example can be run by executing
% <matlab:edit('zynqRadioToneReceiverAD9361AD9364ML')
% zynqRadioToneReceiverAD9361AD9364ML.m>.
if ~exist('prmToneRx','var')
    prmToneRx.SDRDeviceName = 'AD936x';
    prmToneRx.IPAddress = '137.194.172.32';
end
% To update the example for FMCOMMS5, set |prmToneRx.SDRDeviceName| to
% |'FMCOMMS5'|.

%% Transmit a Tone Signal from the FPGA
% Set the Direct Digital Synthesizer (DDS) in the FPGA fabric to transmit a
% complex sinusoid to the RF card. This is provided in the FPGA for testing
% and debugging purposes.

%%
% Create a transmitter System object(TM) to configure the RF card settings. Set
% the RF card to transmit data at a center frequency of 2.4 GHz.
Fs=30e6;
RadioBasebandRate = 30e6;
CenterFrequency = 2.45e9;
RadioFrameLength= 2000;
TxChannel=[1];
RxChannel=[2]; 

sdrTransmitter = sdrtx(prmToneRx.SDRDeviceName, ...
    'IPAddress',       prmToneRx.IPAddress, ...
    'CenterFrequency', CenterFrequency, ...
    'ChannelMapping',TxChannel, ...
    'Gain', 0, ...
    'BasebandSampleRate',RadioBasebandRate);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Filtering

txFilter = comm.RaisedCosineTransmitFilter(...
    'RolloffFactor',0.25, ...
    'FilterSpanInSymbols',25, ...
    'OutputSamplesPerSymbol',2);

rxFilter = comm.RaisedCosineReceiveFilter(...
    'InputSamplesPerSymbol',2, ...
    'DecimationFactor',2, ...
    'FilterSpanInSymbols',25);

%%
% Turn on the properties related to DDS by setting
% |ShowAdvancedProperties| to true. Set the |DataSourceSelect| property of
% sdrTransmitter System object to 'DDS'. Set the tone frequency and scale for
% DDS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Signal generation
basebandSamplingRate=Fs;
test_type='mod';
symbolRate  = 15e6;
basebandOverSampling    = round(basebandSamplingRate/symbolRate);
sps=basebandOverSampling ;
NSamples_BB=20000; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
barker = comm.BarkerCode(... % For preamble
    'Length',13,'SamplesPerFrame',13);
msgLen = 10000;
numFrames = 10;
M = 4; 
frameLen = msgLen/numFrames;

preamble = (1+barker())/2;  % Length 13, unipolar
data = zeros(msgLen, 1);
for idx = 1 : numFrames
    payload = randi([0 M-1],frameLen-barker.Length,1);
    data((idx-1)*frameLen + (1:frameLen)) = [preamble; payload];
end
%data(1:end)=0;
preambleLen=50;
data(1:preambleLen)=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modSize = 4; % Modulation order for 4QAM
k = log2(modSize);  
nQAMSymbols   = round(NSamples_BB/basebandOverSampling);

inSig= randi([0 1],nQAMSymbols,k);      % generate symbols as integer
dataSym = bi2de(inSig);
dataSym(1:preambleLen)=0;

qamSig= qammod(dataSym,modSize,'UnitAveragePower',true);

%%% Baseband (digital) shaping filter %%%
rollOff     = 0.25; % (for RRC filter, single sided output BW is (1+beta)*Rsymb/2 )
symbolSpan  = 50;   % This parameter is related to both the filter length and the attenuation of the stop band
% Instanciate filter
basebandRRC = rcosdesign(rollOff,symbolSpan,basebandOverSampling,'sqrt'); 
basebandSig   = resample(qamSig,basebandOverSampling,1,basebandRRC);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Transmit
sdrTransmitter.ShowAdvancedProperties = true;
sdrTransmitter.BasebandSampleRate = RadioBasebandRate;
fin=5e6;
fin2=0.2e6;
N=2e4;              %Number of points
t=0:1/RadioBasebandRate:(N-1)/RadioBasebandRate;
transmit_signal=zeros(N,length(TxChannel));
%transmit_signal(:,1)=exp(1j*2*pi*fin*t);
transmit_signal=basebandSig;

transmitRepeat(sdrTransmitter,transmit_signal);
pause(1)
%     'GainSource',             'AGC fast Attack', ...
sdrReceiver = sdrrx(prmToneRx.SDRDeviceName, ...
    'IPAddress',              prmToneRx.IPAddress, ...
    'CenterFrequency',        CenterFrequency, ...
    'BasebandSampleRate',     RadioBasebandRate,...
    'GainSource',             'Manual', ...
    'Gain',                   0,...
    'SamplesPerFrame',        N, ...
    'ChannelMapping',         RxChannel, ...
    'OutputDataType',         'double', ...
    'ShowAdvancedProperties', true, ...
    'BypassUserLogic',        true);

RxOut = capture(sdrReceiver,2*N);
RxOutCh1=RxOut;
%Log = dsp.SignalSink;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reception
basebandComplexDigital_fir= resample(RxOutCh1,1,basebandOverSampling,basebandRRC);
% Normalize received symbols to UnitAveragePower (see qammod(inSig)) 
% Trully effective when noise and distortions are not too large compared to the useful signal
basebandComplexDigital_fir_norm= basebandComplexDigital_fir./sqrt(var(basebandComplexDigital_fir));

%Delay compensation 
delay=finddelay(qamSig,basebandComplexDigital_fir_norm);
basebandComplexDigital_fir_norm_trun=basebandComplexDigital_fir_norm(delay+preambleLen+1:delay+(N/2));
%Phase compensation
AngleR=pi/8;
basebandComplexDigital_fir_norm_trun=basebandComplexDigital_fir_norm_trun*exp(-1j*AngleR);
phoffset=mean(angle(basebandComplexDigital_fir_norm_trun)) %Print phoffset
basebandComplex_phasecorrected=basebandComplexDigital_fir_norm_trun*exp(-1j*phoffset);
%basebandComplex_phasecorrected=basebandComplexDigital_fir_norm_trun;
mean(angle(basebandComplex_phasecorrected))
%DEmodulation
qamSigDeSymbols=qamdemod(basebandComplex_phasecorrected,modSize);
DataOut = de2bi(qamSigDeSymbols,k);

resPhzData2 = qamdemod(basebandComplex_phasecorrected,4);

[resPhzTtlErr,resPhzBER] = biterr(dataSym(preambleLen+1:end),resPhzData2);



%%%%RxOut the output of the channels%%%%

release(sdrReceiver);
release(sdrTransmitter);
%release(spectrumScope);
%release(timeScope);
%release(constellation);

%% Conclusion
% In this example, you used SDR Transmitter and Receiver System objects to
% transmit a complex sinusoidal signal from the FPGA and receive it in MATLAB. You
% visualized the received signal in time, frequency and on the complex plane. By
% performing this loopback test, you can confirm that the SDR system is
% setup correctly. You can now proceed to use it in conjunction with
% Communications Toolbox to develop your baseband algorithms and
% verify using real world RF data.
displayEndOfDemoMessage(mfilename)


window_number       = 1;
lineSpec_index      = 1;
fullband_spectrum   = true;



window_number       = window_number+1;
plot_spectrum(RxOutCh1,window_number,...
               Fs,lineSpec_index,fullband_spectrum);
title('Receiver complex recombined output')


if strcmp(test_type, 'mod')
   figure (3)
   subplot(1,2,1)
   plot(qamSig,'d')
   title('constellation at TX')
   subplot(1,2,2)
   plot(basebandComplex_phasecorrected,'d')
   title('constellation at RX')
   % WARNING : the received constellation has almost no sense until phase
   %           distortion and filters delays have been thoroughly analyzed and
   %           repaired
   %           Quick and dirty solution : use FIR filters instead of Butterworth
   %           for 'modulation' case
end
