
%% PART 1: Transmit Diversity vs. Receive Diversity
%
% Using diversity reception is a well-known technique to mitigate the
% effects of fading over a communications link. However, it has mostly been
% relegated to the receiver end. In [ <#13 1> ], Alamouti proposes a
% transmit diversity scheme that offers similar diversity gains, using
% multiple antennas at the transmitter. This was conceived to be more
% practical as, for example, it would only require multiple antennas at the
% base station in comparison to multiple antennas for every mobile in a
% cellular communications system.
%
% This section highlights this comparison of transmit vs. receive diversity
% by simulating coherent binary phase-shift keying (BPSK) modulation over
% flat-fading Rayleigh channels. For transmit diversity, we use two
% transmit antennas and one receive antenna (2x1 notationally), while for
% receive diversity we employ one transmit antenna and two receive antennas
% (1x2 notationally).
%
% The simulation covers an end-to-end system showing the encoded and/or
% transmitted signal, channel model, and reception and demodulation of the
% received signal. It also provides the no-diversity link (single transmit-
% receive antenna case) and theoretical performance of second-order
% diversity link for comparison. It is assumed here that the channel is
% known perfectly at the receiver for all systems. We run the simulation
% over a range of Eb/No points to generate BER results that allow us to
% compare the different systems.

%%
% We start by defining some common simulation parameters
frmLen = 100;       % frame length
numPackets = 1000;  % number of packets
EbNo = 2:1:10;      % Eb/No varying to 20 dB
N = 2;              % maximum number of Tx antennas
M = 2;              % maximum number of Rx antennas
fo=15.36e6;
B=10e6;
P = 4;
k=log2(modSize);
SNR=EbNo+10*log10(fo)-10*log10(B);
SNRdB=EbNo+10*log10(k);

%%
% and set up the simulation.

				% modulation order

% Create comm.OSTBCEncoder and comm.OSTBCCombiner System objects
ostbcEnc = comm.OSTBCEncoder;
ostbcComb = comm.OSTBCCombiner;

% Create two comm.AWGNChannel System objects for one and two receive
% antennas respectively. Set the NoiseMethod property of the channel to
% 'Signal to noise ratio (Eb/No)' to specify the noise level using the
% energy per bit to noise power spectral density ratio (Eb/No). The output
% of the BPSK modulator generates unit power signals; set the SignalPower
% property to 1 Watt.

awgn1Rx = comm.AWGNChannel(...
    'NoiseMethod', 'Signal to noise ratio (Eb/No)', ...
    'SignalPower', 1);
awgn2Rx = clone(awgn1Rx);

% Create comm.ErrorRate calculator System objects to evaluate BER.
errorCalc1 = comm.ErrorRate;
errorCalc2 = comm.ErrorRate;
errorCalc3 = comm.ErrorRate;
errorCalc4 = comm.ErrorRate;

% Since the comm.AWGNChannel System objects as well as the RANDI function
% use the default random stream, the following commands are executed so
% that the results will be repeatable, i.e., same results will be obtained
% for every run of the example. The default stream will be restored at the
% end of the example.
s = rng(55408);

% Pre-allocate variables for speed
H = zeros(frmLen, N, M);
ber_noDiver  = zeros(3,length(EbNo));
ber_Alamouti = zeros(3,length(EbNo));
ber_MaxRatio = zeros(3,length(EbNo));
ber_QAMnoDiver=zeros(3,length(EbNo));
ber_thy2     = zeros(1,length(EbNo));
ber_thy3     = zeros(1,length(EbNo));
%%

% Set up a figure for visualizing BER results
fig = figure; 
grid on;
ax = fig.CurrentAxes;
hold(ax,'on');

ax.YScale = 'log';
xlim(ax,[SNRdB(1), SNRdB(end)]);
ylim(ax,[1e-4 1]);
xlabel(ax,'SNR (dB)');
ylabel(ax,'BER'); 
fig.NumberTitle = 'off';
fig.Renderer = 'zbuffer';
fig.Name = 'Transmit vs. Receive Diversity';
title(ax,'Transmit vs. Receive Diversity 4-QAM');
set(fig, 'DefaultLegendAutoUpdate', 'off');
fig.Position = figposition([15 50 25 30]);

% Loop over several EbNo points
for idx = 1:length(EbNo)
    reset(errorCalc1);
    reset(errorCalc2);
    reset(errorCalc3);
    reset(errorCalc4);
    
    % Set the EbNo property of the AWGNChannel System objects
    awgn1Rx.EbNo = EbNo(idx); 
    awgn2Rx.EbNo = EbNo(idx); 
    % Loop over the number of packets
    for packetIdx = 1:numPackets
        % Generate data vector per frame 
        data = randi([0 P-1], frmLen, 1); 
        
        % Modulate data
        modData = qammod(data,4);  
        modData2= qammod(data,4);

        % Alamouti Space-Time Block Encoder
        encData = ostbcEnc(modData);
        % Create the Rayleigh distributed channel response matrix
        %   for two transmit and two receive antennas
        H(1:N:end, :, :) = (randn(frmLen/2, N, M) + ...
                         1i*randn(frmLen/2, N, M))/sqrt(2);
        %   assume held constant for 2 symbol periods
        H(2:N:end, :, :) = H(1:N:end, :, :);
        
        % Extract part of H to represent the 1x1, 2x1 and 1x2 channels
        H11 = H(:,1,1);
        H21 = H(:,:,1)/sqrt(2);
        H12 = squeeze(H(:,1,:));
        H22=squeeze(H);
        
        % Pass through the channels
        chanOut11 = H11 .* modData;
        chanOut11QAM = H11 .* modData2;
        chanOut21 = sum(H21.* encData, 2);
        chanOut12 = H12 .* repmat(modData, 1, 2);
        
        % Add AWGN
        rxSig11QAM = awgn1Rx(chanOut11QAM);
        rxSig11 = awgn1Rx(chanOut11);
        rxSig21 = awgn1Rx(chanOut21);
        rxSig12 = awgn2Rx(chanOut12);
        
              
        % Alamouti Space-Time Block Combiner
        decData = ostbcComb(rxSig21, H21);

        % ML Detector (minimum Euclidean distance)
        demod11QAM = qamdemod(rxSig11QAM.*conj(H11),4);
        demod11 = qamdemod(rxSig11.*conj(H11),4);
        demod21 = qamdemod(decData,4);
        demod12 = qamdemod(sum(rxSig12.*conj(H12),2),4);
        
        
        % Calculate and update BER for current EbNo value
        %   for uncoded 1x1 system
        ber_noDiver(:,idx)  = errorCalc1(data, demod11);
        %   for Alamouti coded 2x1 system
        ber_Alamouti(:,idx) = errorCalc2(data, demod21);
        %   for Maximal-ratio combined 1x2 system
        ber_MaxRatio(:,idx) = errorCalc3(data, demod12);
         %   for uncoded 1x1 system QAM
        ber_QAMnoDiver(:,idx)  = errorCalc4(data, demod11QAM);
        
    end % end of FOR loop for numPackets

    % Calculate theoretical second-order diversity BER for current EbNo
    ber_thy2(idx) = berfading(EbNo(idx), 'qam', 4, 2);
    ber_thy3(idx) = berfading(EbNo(idx), 'qam', 4, 1);

    % Plot results
    semilogy(ax,SNRdB(1:idx), ber_noDiver(1,1:idx), 'r*', ...
             SNRdB(1:idx), ber_Alamouti(1,1:idx), 'go', ...
             SNRdB(1:idx), ber_MaxRatio(1,1:idx), 'bs', ...
             SNRdB(1:idx), ber_thy2(1:idx), 'm', ...
             SNRdB(1:idx), ber_thy3(1:idx));
    legend(ax,'No Diversity (1Tx, 1Rx)', 'Alamouti (2Tx, 1Rx)',...
           'Maximal-Ratio Combining (1Tx, 2Rx)', ...
           'Theoretical 2nd-Order Diversity', ...
           'Theoretical First-Order Diversity');

    drawnow;
end  % end of for loop for EbNo

% Perform curve fitting and replot the results
fitBER11 = berfit(SNRdB, ber_noDiver(1,:));
fitBER21 = berfit(SNRdB, ber_Alamouti(1,:));
fitBER12 = berfit(SNRdB, ber_MaxRatio(1,:));
semilogy(ax,SNRdB, fitBER11, 'r', SNRdB, fitBER21, 'g', SNRdB, fitBER12, 'b');
hold(ax,'off');

% Restore default stream
rng(s);