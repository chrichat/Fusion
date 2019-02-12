R1=6;
R2=6;
B=1; %For exactly the same subject factors it can be manually set or random analogously with next line
Nt=4; %Downsampling factor
% B=randn(1,1); %Keep the two matrixes analogous 


[tfmri,teeg]=Simulate_time(1,0); % The one is for the same HRF everywhere so the asumption of the model will be correct 
%I do not think that you need something extra for the time beeing, with 0
%then we have a different HRF per source. The second 0 is if we do not want ERPs
%(trial -mode in your case is not needed either).
[sfmri,seeg]=Simulate_space();

load('data/signal_hrf.mat')
hrf=sig_hrf(:,1)';

A1=abs(randn(6,6));
A2=A1;


snr1=[0.1 0.2 0.5 1 1.4 2];
snr2=[  2  5 10 12 14 16];
rowLabels = {'Es. Sour. 1', 'Es. Sour. 2','Es. Sour. 3','Es. Sour. 4','Es. Sour. 5','Es. Sour. 1'};
columnLabels = {' Sour. 1', ' Sour. 2',' Sour. 3',' Sour. 4',' Sour. 5',' Sour. 6'};
teeg(958:end,:)=0;

s=load('C:\Users\hatzi\Dropbox\slall.mat');
corte(1,:)=mean(max(abs(corr(s.sol(1).factors.B,teeg)),[],2));
corse(1,:)=mean(max(abs(corr(s.sol(1).factors.A,seeg')),[],2));
cortf(1,:)=mean(max(abs(corr(s.sol(1).factors.E,tfmri')),[],2));
corsf(1,:)=mean(max(abs(corr(s.sol(1).factors.D,sfmri')),[],2));
NN(1)=mean(mean(s.sol(1).factors.N1));

% s2=load('C:\Users\hatzi\Dropbox\sl2.mat');
corte(2,:)=mean(max(abs(corr(s2.sol(2).factors.B,teeg)),[],2));
corse(2,:)=mean(max(abs(corr(s.sol(2).factors.A,seeg')),[],2));
cortf(2,:)=mean(max(abs(corr(s.sol(2).factors.E,tfmri')),[],2));
corsf(2,:)=mean(max(abs(corr(s.sol(2).factors.D,sfmri')),[],2));
NN(2)=mean(mean(s.sol(2).factors.N1));

corte(3,:)=mean(max(abs(corr(s.sol(3).factors.B,teeg)),[],2));
corse(3,:)=mean(max(abs(corr(s.sol(3).factors.A,seeg')),[],2));
cortf(3,:)=mean(max(abs(corr(s.sol(3).factors.E,tfmri')),[],2));
corsf(3,:)=mean(max(abs(corr(s.sol(3).factors.D,sfmri')),[],2));
NN(3)=mean(mean(s.sol(3).factors.N1));


corte(4,:)=mean(max(abs(corr(s.sol(4).factors.B,teeg)),[],2));
corse(4,:)=mean(max(abs(corr(s.sol(4).factors.A,seeg')),[],2));
cortf(4,:)=mean(max(abs(corr(s.sol(4).factors.E,tfmri')),[],2));
corsf(4,:)=mean(max(abs(corr(s.sol(4).factors.D,sfmri')),[],2));
NN(4)=mean(mean(s.sol(4).factors.N1));


s=load('C:\Users\hatzi\Dropbox\slambdasol3.mat');
corte(3,:)=mean(max(abs(corr(s.sol.factors.B,teeg)),[],2));
corse(3,:)=mean(max(abs(corr(s.sol.factors.A,seeg')),[],2));
cortf(3,:)=mean(max(abs(corr(s.sol.factors.E,tfmri')),[],2));
corsf(3,:)=mean(max(abs(corr(s.sol.factors.D,sfmri')),[],2));
NN(3,:)=mean(mean(s.sol.factors.N1));

s=load('C:\Users\hatzi\Dropbox\slambdasol4.mat');
corte(4,:)=mean(max(abs(corr(s.sol.factors.B,teeg)),[],2));
corse(4,:)=mean(max(abs(corr(s.sol.factors.A,seeg')),[],2));
cortf(4,:)=mean(max(abs(corr(s.sol.factors.E,tfmri')),[],2));
corsf(4,:)=mean(max(abs(corr(s.sol.factors.D,sfmri')),[],2));
NN(4,:)=mean(mean(s.sol.factors.N1));


s=load('C:\Users\hatzi\Dropbox\slambdasol5.mat');
corte(5,:)=mean(max(abs(corr(s.sol.factors.B,teeg)),[],2));
corse(5,:)=mean(max(abs(corr(s.sol.factors.A,seeg')),[],2));
cortf(5,:)=mean(max(abs(corr(s.sol.factors.E,tfmri')),[],2));
corsf(5,:)=mean(max(abs(corr(s.sol.factors.D,sfmri')),[],2));
NN(5,:)=mean(mean(s.sol.factors.N1));


s=load('C:\Users\hatzi\Dropbox\slambdasol6.mat');
corte(6,:)=mean(max(abs(corr(s.sol.factors.B,teeg)),[],2));
corse(6,:)=mean(max(abs(corr(s.sol.factors.A,seeg')),[],2));
cortf(6,:)=mean(max(abs(corr(s.sol.factors.E,tfmri')),[],2));
corsf(6,:)=mean(max(abs(corr(s.sol.factors.D,sfmri')),[],2));
NN(6,:)=mean(mean(s.sol.factors.N1));

s=load('C:\Users\hatzi\Dropbox\slambdasol7.mat');
corte(7,:)=mean(max(abs(corr(s.sol.factors.B,teeg)),[],2));
corse(7,:)=mean(max(abs(corr(s.sol.factors.A,seeg')),[],2));
cortf(7,:)=mean(max(abs(corr(s.sol.factors.E,tfmri')),[],2));
corsf(7,:)=mean(max(abs(corr(s.sol.factors.D,sfmri')),[],2));
NN(7,:)=mean(mean(s.sol.factors.N1));

s=load('C:\Users\hatzi\Dropbox\slambdasol8.mat');
corte(8,:)=mean(max(abs(corr(s.sol.factors.B,teeg)),[],2));
corse(8,:)=mean(max(abs(corr(s.sol.factors.A,seeg')),[],2));
cortf(8,:)=mean(max(abs(corr(s.sol.factors.E,tfmri')),[],2));
corsf(8,:)=mean(max(abs(corr(s.sol.factors.D,sfmri')),[],2));
NN(8,:)=mean(mean(s.sol.factors.N1));


