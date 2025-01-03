%% Gaussian noise custom function
function [y,n] = add_awgn_noise(x,SNR_dB)
    %[y,n]=awgn_noise(x,SNR) adds AWGN noise vector to signal ’x’ to generate a
    %resulting signal vector y of specified SNR in dB. It also returns the
    %noise vector ’n’ that is added to the signal x
    L=length(x);
    SNR = 10^(SNR_dB/10); %SNR to linear scale
    Esym=sum(abs(x).^2)/(L); %Calculate actual symbol energy N0=Esym/SNR; %Find the noise spectral density
    N0 = Esym/SNR
    if(isreal(x)),
        noiseSigma = sqrt(N0);%Standard deviation for AWGN Noise when x is real
        n = noiseSigma*randn(L,1);%computed noise else
        noiseSigma = sqrt(N0/2);%Standard deviation for AWGN Noise when x is complex n = noiseSigma*(randn(1,L)+1i*randn(1,L));%computed noise
    end
    y = x + n; %received signal
end