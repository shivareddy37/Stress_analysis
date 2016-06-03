%% the function performs all major operations namely FILTERING ,ZERO MEANING,ARTIFACT REMEOVAL,FEATURE GENRATION AND FEATURE REDUCTION
% output of thefuction helps us to create the feature matrix
%% AUTHORS : SHIVA KARTHIK REDDY KOYYA , RAMESH NAIR
%% start of function 
% some values have been hard coded have to be taken care of while
% modification


function [ weight,meanMatrix,l ] = allMajor( stress_class, non_stress_class )

%% Filtering of data
SData=stress_class;
NSData=non_stress_class;

fs=256; % sampling rate  % power line noise frequency
Wn=[0.4 64]./(fs/2); %normalised cutoff frequency for the bandpass filter
filterOrder=2;
[b,a]=butter(filterOrder,Wn);% genrates filter cofficients
filteredDataS=filtfilt(b,a,SData);
filteredDataNS=filtfilt(b,a,NSData);

wo = 60/(fs/2);
bw = wo/60;
[bn,an] = iirnotch(wo,bw); %notch filter at 60 Hz
filteredDataS=filter(bn,an,filteredDataS);
filteredDataNS=filter(bn,an,filteredDataNS);


%% zero meaning data sample and normalising the features
[m,n]=size(filteredDataS);
DataS=zeros(m,n);
DataNS=zeros(m,n);


%zero meaning
for i=1:m
    DataS(i,:)=filteredDataS(i,:)-mean(filteredDataS(i,:));
    DataNS(i,:)=filteredDataNS(i,:)-mean(filteredDataNS(i,:));
end
% normalising features

%  for i=1:n
%      DataS(:,i)=(DataS(:,i)-min(DataS(:,i)))/(max(DataS(:,i))-min(DataS(:,i)));
%      DataNS(:,i)=(DataNS(:,i)-min(DataNS(:,i)))/(max(DataNS(:,i))-min(DataNS(:,i)));
%  end


%% analysing the  filtred , zero meaning data for both classes
figure(2)
plot(DataS(1:256*5,1),'r');
hold on
plot(DataNS(1:256*5,1),'b');
grid on
title('plot of filtered and Zero meaned Data for both classes channel 1 ');
legend('stress','nonstress');
hold off




%% Artificat removal
%ICA and LWT

%performing ICA using Fast ICA function
figure(3)
[icasig1, A1, W1] = fastica (DataS');
figure(4)
[icasig2,A2,W2]=fastica(DataNS');

% performing 4 level decomposition using LWT
for i=1:4
    [CA1(i,:),CD1(i,:)]=lwt(icasig1(i,:),'db4',4);
    [CA2(i,:),CD2(i,:)]=lwt(icasig2(i,:),'db4',4);
end

% power spectrumn of the independent components
% figure(5)
% T = 1/fs;                     % Sample time
% L = length(icasig1);          % Length of signal
% t = (0:L-1)*T;                % Time vector
% NFFT = 2^nextpow2(L);         % Next power of 2 from length of y
% Y = fft(icasig1(1:100),NFFT)/L;
% f = fs/2*linspace(0,1,NFFT/2+1);
% 
% % Plot single-sided amplitude spectrum.
% plot(f,2*abs(Y(1:NFFT/2+1)))
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')

% thresholding
for i=1:4
    CD1(i,:)=wthresh(CD1(i,:),'h',100);
    CD2(i,:)=wthresh(CD2(i,:),'h',100);
end

% ILWT
for i=1:4
    icasig1(i,:)=ilwt(CA1(i,:),CD1(i,:),'db4',4);
    icasig2(i,:)=ilwt(CA2(i,:),CD2(i,:),'db4',4);
end

% reconstruction of noise free signal
StressClass=pinv(W1)*icasig1;
nonStressClass=pinv(W2)*icasig1;
%plot of reconstructed signal free from artifacts
figure(6)
plot(StressClass(1,1:256*5)','r');
hold on
plot(nonStressClass(1,1:256*5)','b');
grid on
title('plot of rawdata for both classes channel 1 ');
legend('stress','nonstress');
hold off

%% feature genration and extraction
StressClass=StressClass';
nonStressClass=nonStressClass';
l=length(StressClass);
sampleLength=l/56;
B_S=reshape(StressClass,sampleLength,5,[]);
B_NS=reshape(nonStressClass,sampleLength,5,[]);
mean_S=zeros(8,5,8);
mean_NS=zeros(8,5,8);
rms_S=zeros(8,5,8);
rms_NS=zeros(8,5,8);
std_S=zeros(8,5,8);
std_NS=zeros(8,5,8);
Feature_S=[];
Feature_NS=[];

for i=1:56
    for j=1:5
        [C_S,L_S]=wavedec(B_S(:,j,i),4,'db4');
        [C_NS,L_NS]=wavedec(B_NS(:,j,i),4,'db4');
        A1_S=appcoef(C_S,L_S,'db4',1);
        A2_S=appcoef(C_S,L_S,'db4',2);
        A3_S=appcoef(C_S,L_S,'db4',3);
        A4_S=appcoef(C_S,L_S,'db4',4);
        D1_S=detcoef(C_S,L_S,1);
        D2_S=detcoef(C_S,L_S,2);
        D3_S=detcoef(C_S,L_S,3);
        D4_S=detcoef(C_S,L_S,4);
        
        
        mean_S(:,j,i)=[mean(abs(A1_S));mean(abs(A2_S));mean(abs(A3_S));mean(abs(A4_S));mean(abs(D1_S));...
            mean(abs(D2_S));mean(abs(D3_S));mean(abs(D4_S))];
        rms_S(:,j,i)=[rms(A1_S);rms(A2_S);rms(A3_S);rms(A4_S);...
            rms(D1_S);rms(D2_S);rms(D3_S);rms(D4_S);];
        std_S(:,j,i)=[std(A1_S);std(A2_S);std(A3_S);std(A4_S);std(D1_S);std(D2_S);...
            std(D3_S);std(D4_S);];
        
        
        A1_NS=appcoef(C_NS,L_NS,'db4',1);
        A2_NS=appcoef(C_NS,L_NS,'db4',2);
        A3_NS=appcoef(C_NS,L_NS,'db4',3);
        A4_NS=appcoef(C_NS,L_NS,'db4',4);
        D1_NS=detcoef(C_NS,L_NS,1);
        D2_NS=detcoef(C_NS,L_NS,2);
        D3_NS=detcoef(C_NS,L_NS,3);
        D4_NS=detcoef(C_NS,L_NS,4);
        
        
        mean_NS(:,j,i)=[mean(abs(A1_NS));mean(abs(A2_NS));mean(abs(A3_NS));mean(abs(A4_NS));mean(abs(D1_NS));...
            mean(abs(D2_NS));mean(abs(D3_NS));mean(abs(D4_NS))];
        rms_NS(:,j,i)=[rms(A1_NS);rms(A2_NS);rms(A3_NS);rms(A4_NS);...
            rms(D1_NS);rms(D2_NS);rms(D3_NS);rms(D4_NS);];
        std_NS(:,j,i)=[std(A1_NS);std(A2_NS);std(A3_NS);std(A4_NS);std(D1_NS);std(D2_NS);...
            std(D3_NS);std(D4_NS);];
        Feature_S=[Feature_S mean_S(:,j,i)' rms_S(:,j,i)' std_S(:,j,i)'];
        Feature_NS=[Feature_NS mean_NS(:,j,i)' rms_NS(:,j,i)' std_NS(:,j,i)'];
    end
    
end

Feature_S=reshape(Feature_S,[],120);
Feature_NS=reshape(Feature_NS,[],120);
S_Label=ones(56,1);
NS_Label=zeros(56,1);
DATA=[Feature_S;Feature_NS];
Label=[S_Label;NS_Label];
l=Label;





%% Run PCA

fprintf('Performing PCA...');

% Variance to maintain
per = 0.99;

% Get mean, and create matrix
meanVals    = mean(DATA,1);
meanMat     = ones(size(DATA))*diag(meanVals);

% Zero-mean
meanMatrix = (DATA - meanMat);

% Co-variance
sigma = cov(meanMatrix);

% Perform SVD
[U, E, V] = svd(sigma);

% Get the sum of the eigenvalues
tra = trace(E);

% Find the smallest number of eigenvalues to maintain specified variance
sm = 0;
for k=1:length(E)
    sm = sm + E(k,k);
    if( sm/tra > per )
        break
    end
end

% PCA matrix
weight = V(:,1:k);

fprintf('finished!, #PCs = %d\n',k);

end

