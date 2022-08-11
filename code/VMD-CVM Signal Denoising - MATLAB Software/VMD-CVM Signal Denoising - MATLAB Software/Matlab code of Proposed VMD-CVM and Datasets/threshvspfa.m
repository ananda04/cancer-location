function [PfvThvec,ind_m,disn_m]= threshvspfa(imfvec,N)

%% Estimation of noise EDF from rejected modes

MC=length(imfvec);
for j=1:floor(MC/N)             % loop for all windows
    ch=imfvec(N*(j-1)+1:N*j);    % pick the jth window
    [temp,tind]=ecdf(ch);       % calculate ECDF
    tv(:,j)=temp(2:end,1);      % store value in tv
    ti(:,j)=tind(2:end,1);      % store index in ti
end

disn_m=mean(tv,2);   % take mean value of ECDF values
ind_m=mean(ti,2);

g=0;

N=32;
thresh_min=0.001;
inc=0.001;
thresh_max=20;

%% Threshold versus Pfa curve estimation from rejected modes

threshvec = thresh_min:inc:thresh_max;

pfavec=zeros(length(threshvec),1);            % vector for storing Pfa
% imfvec=zeros(s,2,length(threshvec));          % vector for storing Pfa vs Threshold values for all IMFs

% for noofimf=IMF_start:NIMF                % for the first NIMF
% noofimf=3;
    i=1;
    g=g+1;
    x=imfvec;           % pick an IMF
    disnref=disn_m;  % pick corresponding ECDF value
    indref=ind_m;   % pick corresponding ECDF index
for thresh= threshvec         % vary threshold
    
    count_detection=0;
   
    for litcount=1:floor(MC/N)  % loop for all windows
        z=cdfcalc(sort(x(1,N*(litcount-1)+1:N*litcount)),disnref,indref);   % calculated F_eta (x)
        
        test=cvm(z,N); % CVM statistic
       
        if test > thresh                                                    % compare with threshold
            count_detection = count_detection + 1;                          % increment detection count
        end
    end
    
    Pfa = count_detection/floor(MC/N);                                      % calculate Pfa
    pfavec(i,1)=Pfa;                                                        % store Pfa in the vector
    i=i+1;
    
    if Pfa < 0.000005
        break;
    end
end


PfvThvec=[threshvec;pfavec'];                                    % store Pfa vs Threshold values for each IMF here for later use