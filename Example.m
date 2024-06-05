clear;
close all;
clc

%% signal 1
load('Sig1.mat');
fs1=100;
n1=length(Sig1);
time1=(1:n1)/fs1;
fre1=(fs1/2)/(n1/2):(fs1/2)/(n1/2):(fs1/2);  % max(fre) <= fs/2

%% signal 2
load('Sig2.mat');
n2 = length(Sig2);
fs2 = 2501;
time2 = (1:n2)/fs2;
fre2=(fs2/2)/(n2/2):(fs2/2)/(n2/2):(fs2/2);

%% different Chirp rate estimators
[tfr_f1_1,tfr_t1_1,~,~,~,Ts_STFT1] = TFET(Sig1',40,1,'TFET');
[tfr_f1_2,tfr_t1_2,~,~,~,~] = TFET(Sig1',40,2,'TFET');
[tfr_f1_3,tfr_t1_3,~,~,~,~] = TFET(Sig1',40,3,'TFET');
[tfr_f1_4,tfr_t1_4,~,~,~,~] = TFET(Sig1',40,4,'TFET');
[tfr_f1_5,tfr_t1_5,~,~,TFR1_TFET,~] = TFET(Sig1',40,5,'TFET');
[tfr_f1_6,tfr_t1_6,~,~,TFR1_TFMST,~] = TFMST_Y(Sig1',40);

M=max(max(abs(Ts_STFT1)));
Ts_STFT1(find(abs(Ts_STFT1)<0.3*M))=0;

Ts_STFT1_1 = zeros(size(Ts_STFT1));
Ts_STFT1_2 = Ts_STFT1;
Ts_STFT1_1(:,96:106) = 2;
Ts_STFT1_2(find(abs(Ts_STFT1_2)~=0))=1; 
Label1 = Ts_STFT1_2 - Ts_STFT1_1 ;

Pre1 = zeros([6,size(Ts_STFT1)]);
ConfusionMat1 = zeros(6,4);
Metric1 = zeros(6,2);  % mPA and mIoU
for i = 1:6
    F = eval(['tfr_f1_',num2str(i)]);
    T = eval(['tfr_t1_',num2str(i)]);
    F(find(abs(F)~=0))=1; 
    T(find(abs(T)~=0))=-1; 
    Pre = F+T;
    [TP, FN, FP, TN] = confusionMatrix(-1/2*Ts_STFT1_1, Label1 - (-1/2*Ts_STFT1_1), Pre);
    mPA = 1/2*(TP / (TP + FP) + TN / (TN + FN));
    mIoU = 1/2*(TP / (TP + FP + FN) + TN / (TN + FN + FP));
    
    Pre1(i,:,:) = Pre;
    ConfusionMat1(i,:) = [TP, FN, FP, TN];
    Metric1(i,:) = [mPA, mIoU];
end

%% segmentation results on signal 1
figure;
for i = 1:6
    subplot(3,2,i)
    imagesc(time1,fre1,squeeze(Pre1(i,:,:)));axis xy;
end
colormap Colorcube;

figure;
for i = 1:2:12
    subplot(3,4,i)
    imagesc(time1,fre1,abs(eval(['tfr_f1_',num2str(ceil(i/2))])));axis xy;
    subplot(3,4,i+1)
    imagesc(time1,fre1,abs(eval(['tfr_t1_',num2str(ceil(i/2))])));axis xy;
end

%% different Chirp rate estimators
[tfr_f2_1,tfr_t2_1,~,~,~,Ts_STFT2] = TFET(Sig2',225,1,'TFET');
[tfr_f2_2,tfr_t2_2,~,~,~,~] = TFET(Sig2',225,2,'TFET');
[tfr_f2_3,tfr_t2_3,~,~,~,~] = TFET(Sig2',225,3,'TFET');
[tfr_f2_4,tfr_t2_4,~,~,~,~] = TFET(Sig2',225,4,'TFET');
[tfr_f2_5,tfr_t2_5,~,~,TFR2_TFET,~] = TFET(Sig2',225,5,'TFET');
[tfr_f2_6,tfr_t2_6,~,~,TFR2_TFMST,~] = TFMST_Y(Sig2',225);

M=max(max(abs(Ts_STFT2)));
Ts_STFT2(find(abs(Ts_STFT2)<0.3*M))=0;

Ts_STFT2_1 = zeros(size(Ts_STFT2));
Ts_STFT2_2 = zeros(size(Ts_STFT2));
Ts_STFT2_2(490:810,:) = Ts_STFT2(490:810,:);
Ts_STFT2_1(find(abs(Ts_STFT2)~=0)) = 1;
Ts_STFT2_2(find(abs(Ts_STFT2_2)~=0)) = 2; 
Label2 = Ts_STFT2_1 - Ts_STFT2_2;

Pre2 = zeros([6,size(Ts_STFT2)]);
ConfusionMat2 = zeros(6,4);
Metric2 = zeros(6,2);  % mPA and mIoU
for i = 1:6
    F = eval(['tfr_f2_',num2str(i)]);
    T = eval(['tfr_t2_',num2str(i)]);
    F(find(abs(F)~=0))=1; 
    T(find(abs(T)~=0))=-1; 
    Pre = F+T;
    [TP, FN, FP, TN] = confusionMatrix(-1/2*Ts_STFT2_2, Label2 - (-1/2*Ts_STFT2_2), Pre);
    mPA = 1/2*(TP / (TP + FP) + TN / (TN + FN));
    mIoU = 1/2*(TP / (TP + FP + FN) + TN / (TN + FN + FP));
    
    Pre2(i,:,:) = Pre;
    ConfusionMat2(i,:) = [TP, FN, FP, TN];
    Metric2(i,:) = [mPA, mIoU];
end

%% segmentation results on signal 2
figure;
for i = 1:6
    subplot(3,2,i)
    imagesc(time2,fre2,squeeze(Pre2(i,:,:)));axis xy;
end
colormap Colorcube;

figure;
for i = 1:2:12
    subplot(3,4,i)
    imagesc(time2,fre2,abs(eval(['tfr_f2_',num2str(ceil(i/2))])));axis xy;
    subplot(3,4,i+1)
    imagesc(time2,fre2,abs(eval(['tfr_t2_',num2str(ceil(i/2))])));axis xy;
end

%% A comparison between TFET and TFMST
figure;
subplot(2,2,1);
imagesc(time1,fre1,abs(TFR1_TFET));axis xy;
xlabel('Time / s');
ylabel('Fre / Hz');
subplot(2,2,2);
imagesc(time1,fre1,abs(TFR1_TFMST));axis xy;
xlabel('Time / s');
ylabel('Fre / Hz');
subplot(2,2,3);
imagesc(time2,fre2,abs(TFR2_TFET));axis xy;
xlabel('Time / s');
ylabel('Fre / Hz');
subplot(2,2,4);
imagesc(time2,fre2,abs(TFR2_TFMST));axis xy;
xlabel('Time / s');
ylabel('Fre / Hz');