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
if ~exist('prmToneTx','var')
    prmToneTx.SDRDeviceName = 'AD936x';
    prmToneTx.IPAddress = '137.194.172.35';
end

if ~exist('prmToneRx','var')
    prmToneRx.SDRDeviceName = 'AD936x';
    prmToneRx.IPAddress = '137.194.172.35';
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

RadioBasebandRate = 30e6;
CenterFrequency = 2.4e9;
RadioFrameLength = 4000;
TxChannel=[1,2];
RxChannel=[1,2];      

sdrTransmitter = sdrtx(prmToneTx.SDRDeviceName, ...
    'IPAddress',       prmToneTx.IPAddress, ...
    'CenterFrequency', CenterFrequency, ...
    'ChannelMapping',TxChannel, ...
    'Gain', 0);

%%
% Turn on the properties related to DDS by setting
% |ShowAdvancedProperties| to true. Set the |DataSourceSelect| property of
% sdrTransmitter System object to 'DDS'. Set the tone frequency and scale for
% DDS.

sdrTransmitter.ShowAdvancedProperties = true;
sdrTransmitter.BasebandSampleRate = RadioBasebandRate;
fin_or=1e6;
BW=10e6;
Fs=RadioBasebandRate;
fin2=3e6;
N=2^16;              %Number of points
t=0:1/RadioBasebandRate:(N-1)/RadioBasebandRate;

Bin_in=round(fin_or/Fs*N);
fin=Bin_in*Fs/N;

transmit_signal=zeros(N,length(TxChannel));
transmit_signal(:,1)=exp(1j*2*pi*fin*t);

transmit_signal(:,2)=exp(1j*2*pi*fin2*t);
transmitRepeat(sdrTransmitter,transmit_signal);



sdrReceiver = sdrrx(prmToneRx.SDRDeviceName, ...
    'IPAddress',              prmToneTx.IPAddress, ...
    'CenterFrequency',        CenterFrequency, ...
    'BasebandSampleRate',     RadioBasebandRate,...
    'GainSource',             'Manual',...
    'Gain',                   0, ...
    'SamplesPerFrame',        RadioFrameLength, ...
    'ChannelMapping',         RxChannel, ...
    'OutputDataType',         'double', ...
    'ShowAdvancedProperties', true, ...
    'BypassUserLogic',        true);
RxOut = capture(sdrReceiver,2*N);

RxOut1=RxOut(:,1);
RxOut2=RxOut(:,2);
%%%%RxOut the output of the channels%%%%

plot_spectrum(RxOut1,1,RadioBasebandRate,1);


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