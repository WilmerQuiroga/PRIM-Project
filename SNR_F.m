function SNR_freq=SNR_F(x_PSD,N,fin,continuousTimeSamplingRate, BW)
  % beggining and end of the useful band (here, Nyquist band)  
sig_bin=round(fin/continuousTimeSamplingRate*N)+1;
binBW=[1, round(BW/continuousTimeSamplingRate*N)+1];
 
nintsig   = 180000;              % Number of points around the useful signal to integrate
PS_freq   = sum(x_PSD(sig_bin-nintsig:sig_bin+nintsig));
PN_freq   = sum(x_PSD(binBW(1):binBW(2)))-PS_freq;                   
SNR_freq  = 10*log10(PS_freq/PN_freq); 
end