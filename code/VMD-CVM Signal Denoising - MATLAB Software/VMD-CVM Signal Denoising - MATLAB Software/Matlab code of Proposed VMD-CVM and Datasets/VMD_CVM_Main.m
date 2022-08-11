clc; 
clear all
close all
%% Input Parameter

win_len=32; % Window Length
NIMF=10;     % IMFs to consider
Np=36;       % No of consecutive iterations that must detect for the detection to hold
N_mon=10;     
sigL = 12;

%% Select Pfa using a decaying function e^(-k+1). 

% opt_Pfa=[1 0.6 0.1 0.05 0.01 0.005 0.001 0.001 0.0005 0.0001];
opt_Pfa = (1/2.71).^(0:NIMF-1);
% opt_Pfa =[1 0.5 0.2 0.05 0.01 0.005 0.001 0.0005 0.0001 0.00001];

%% Input Signal
disp('Enter 1 for Blocks; 2 for Bumps; 3 for Heavy Sine; 4 for Doppler;');
disp('5 for Sofar; 6 for Tai Chi and 7 to denoise the Raw ECG signal');
sig = input('Enter you choice =: ');

if sig == 1
    a = wnoise(sig,sigL)';
elseif sig == 2
    a = wnoise(sig,sigL)';
elseif sig == 3
    a = wnoise(sig,sigL)';
elseif sig == 4
    a = wnoise(sig,sigL)';
elseif sig == 5
         load sofar.mat;
    a = imresize(x,[1024,1]);
elseif sig == 6
    load taichi.mat;
    a = imresize(taichi_final(:,2),[1024,1]);
elseif sig == 7
    load ECGsig.mat;
    a=noisy(1:2000);
    a1=orig(1:2000);
end
figure;plot(a);

%% Generating the noisy signal for given input SNR

if (sig~=7)
iSNR=input('Enter your choice of input SNR =: ');
else
    iSNR= 20;
end
if(iSNR<=0)
    SNRii=0;
else
    SNRii=iSNR;
end
%%%%%%%%%%%%%%%%%% SNR = -2db
        sigma=0.1+(0.3/(SNRii+1));
 x = -1/2:1/((win_len)-1):1/2;
g = exp( -(x.^2)/(2*sigma^2) );
g = g / sum(sum(g));

snri = iSNR;
f=awgn(a,snri,'measured'); 

%% Proposed Approach

[imf rec y]=Prop_VMD_CVM(a,f,win_len,NIMF,opt_Pfa,Np); % Proposed denoising algorithm

y1 = conv (y, g, 'same'); % Post Processing

%% Computing Performance measures and plotting the denoised signal

osnry=snr(a,y);
osnry1=snr(a,y1);
if osnry < osnry1
    osnry=osnry1;
    y = y1;
end
oSNR=osnry
figure;
plot(a);hold on
plot(y);