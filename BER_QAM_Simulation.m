M = 16;                 % Modulation order
M2=4;
%k = log2(M);            % Bits per symbol
EbNoVec = (5:15)';      % Eb/No values (dB)

SNRdB=[5.2804 5.3667 5.697 6.2024 7.1965 7.9949 8.6001 10.8099 12.1878 13.8013];
berEst = zeros(size(SNRdB));
berEst2 = zeros(size(SNRdB));
EbNo = zeros(size(SNRdB));
EbNo2 = zeros(size(SNRdB));
modSize= 16; % Modulation order for 4QAM
modSize2= 4;
k = log2(modSize);
k2 = log2(modSize2);
fo=15.36e6;
B=10e6;
numSymPerFrame = 100;   % Number of QAM symbols per frame

%berEst = zeros(size(EbNoVec));

for n = 1:length(SNRdB)
    EbNo(n)=SNRdB(n) -10*log10(k);
    EbNo2(n)=SNRdB(n) -10*log10(k2);
    % Reset the error and bit counters
    numErrs = 0;
    numBits = 0;
    
    numErrs2 = 0;
    numBits2 = 0;
    
    while numErrs < 1000 && numBits < 10e6
        % Generate binary data and convert to symbols
        dataIn = randi([0 1],numSymPerFrame,k);
        dataSym = bi2de(dataIn);
        % QAM Modulation 
        dataIn2 = randi([0 1],numSymPerFrame,k2);
        dataSym2 = bi2de(dataIn2);
        
        % QAM modulate using 'Gray' symbol mapping
        txSig = qammod(dataSym,M);
        txSig2 = qammod(dataSym2,M2);
        % Pass through AWGN channel
        rxSig = awgn(txSig,SNRdB(n),'measured');
        rxSig2 = awgn(txSig2,SNRdB(n),'measured');
        
        % Demodulate the noisy signal
        rxSym = qamdemod(rxSig,M);
        rxSym2 = qamdemod(rxSig2,M2);
        % Convert received symbols to bits
        dataOut = de2bi(rxSym,k);
        dataOut2 = de2bi(rxSym2,k2);
        % Calculate the number of bit errors
        nErrors = biterr(dataIn,dataOut);
        nErrors2 = biterr(dataIn2,dataOut2);
        
        % Increment the error and bit counters
        numErrs = numErrs + nErrors;
        numBits = numBits + numSymPerFrame*k;
        
        numErrs2 = numErrs2 + nErrors2;
        numBits2 = numBits2 + numSymPerFrame*k2;
    end
    
    % Estimate the BER
    berEst(n) = numErrs/numBits;
    berEst2(n) = numErrs2/numBits2;
end


berTheory = berawgn(EbNo,'qam',M);
berTheory2 = berawgn(EbNo2,'qam',M2);

figure (1)
semilogy(SNRdB,berEst)
hold on
semilogy(SNRdB,berTheory,'r')
semilogy(SNRdB,berEst2,'c')
semilogy(SNRdB,berTheory2)
grid
legend('Estimated BER 16-QAM','Theoretical BER 16-QAM','Estimated BER QAM','Theoretical BER QAM' )
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
title('BER vs SNR')

