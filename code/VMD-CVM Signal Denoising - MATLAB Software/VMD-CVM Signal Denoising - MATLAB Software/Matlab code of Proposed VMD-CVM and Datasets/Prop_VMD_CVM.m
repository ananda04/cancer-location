
function [imf1 rec sigrec]=Prop_VMD_CVM(ref,noisy,N,NIMF,Pf,Np) 
a=ref; f=noisy;

pts=length(a); % data length

sigma=std(a-f); % Estimate total noise variance for normalising IMFs

% some sample parameters for VMD
alpha = 2000;        % moderate bandwidth constraint
tau = 0;            % noise-tolerance (no strict fidelity enforcement)
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly
tol = 1e-7;

[imf1, imf_hat, omega] = VMD(f, alpha, tau, NIMF, DC, init, tol);
imf=imf1/sigma;
if size(imf,1)<NIMF
    NIMF=size(imf,1)-1;
    warning('NIMF found to be less than the size of actual imf');
end

rec=zeros(size(imf)); % Holds recoved signal on IMF level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dist=[];
% hs = ksdensity(noisy);
[nEdf,nInd]=ecdf(noisy);

for imfcnt=1:NIMF

tempx=imf(imfcnt,:);
[Edf,Ind]=ecdf(tempx);
z=cdfcalc(sort(tempx),nEdf,nInd);    
Dist(imfcnt)=cvm(z,pts); % CVM statistic
end
D=abs(diff(Dist));
n=find(D==max(D));
while n <= 5
    D(n)=0;
    n=find(D==max(D));
end
ni=10-n;
if (ni < 3)
    ni=3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imfvec=[];
for i=NIMF-ni:NIMF               % loop for all IMFs
    ch1=imf(i,:); 
    imfvec = [imfvec ch1];
end
[PfvThvec,ind_m,disn_m]= threshvspfa(imfvec,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for imfcnt=1:NIMF-3 
    TH=PfvThvec(1,:)';
    PF=PfvThvec(2,:)';
    temp=imf(imfcnt,:); % Current IMF values
    dref=disn_m;
    indref=ind_m;
    [~,indpf]=min(abs(Pf(imfcnt)-PF)); % Matches the optimal Pfa calulated for current IMF with available PFAs
    
    if ~isempty(indpf)
        thresh = TH(indpf); %extract the threshold value for current IMF
    else
        error('Pfa value given is not valid');
    end
    
    booln=zeros(size(temp)); % Holds the flag for all points for detection. e.g 1 for the detection of signal
    tau_thr = zeros(size(temp));
    %%%%%%%%%%%% Overlapping windows 
    for jj = N+1:pts-N % Loop through all the windows of a single IMF    
        x=temp(1,jj-floor(N/2)+1:jj+ceil(N/2)); % window on an IMF
        z=cdfcalc(sort(x),dref,indref); % dref is F_{eta}(x)
        
        test=cvm(z,N); % CVM statistic
            
        if test > thresh % statistic > threshold: signal present
            booln(jj)=1;
%             tau_thr(jj) = (thresh^2)/(test^2);
        end
    end
    
    % Consider detection only if it happens for atleast for length N; removes impulse-like noise! 
    D = diff([0,booln,0]);
    bg = find(D == 1); ed = find(D == -1) - 1; % bg (ed) contains the beginning/end locations of clusters of 1's 
    for ii=1:length(bg)
        if(ed(ii)-bg(ii)<Np) % If the length of clusters of 1's is less than 'N' window length
            booln(bg(ii):ed(ii))=zeros(ed(ii)-bg(ii)+1,1); % No detection of signal there!
        end
    end
    rec1=temp.*booln;
    rec(imfcnt,:)= rec1;
end
rec(1,:)=imf(1,:);
rec = sigma*rec;
sigrec = sum(rec,1)';
end