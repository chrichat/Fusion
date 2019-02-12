function [tfmri,teeg,sig_erp]=Simulate_time(hrf,erp,giv_hrf)
% HRF=1 means the same hrf for all sources
% Hrf=0 means different hrf per source specified from file signal_hrf.mat

load('data/signal_erp.mat')
load('data/signal_hrf.mat')

%% Assing given hrf
sig_hrf(:,1)=giv_hrf;

%%
sig_erp(1,:)=circshift(sig_erp(1,:),[0,4]);
sig_erp(1,1:15)=-0.5*(sig_erp(2,1:15));
sig_erp(1,:)=circshift(sig_erp(1,:),[0,4]);
sig=zeros(size(sig_erp(1,:)));
sig(1,21:34)=sig_erp(2,1:14);
sig(1,21:27)=0.5*sig_erp(4,1:7);
sig_erp(5,:)=sig;
sig_erp(4,:)=0.7*circshift(sig_erp(1,:),[0,7]);
sig_erp(4,1:14)=0.3*sig_erp(2,1:14);

%%%TOEPLITZ FOR CONV%%%
% 
% x=sigt_eeg(1,:);
% hrf=ht';
% a1=conv(hrf,x);
% 
% t_hrf=toeplitz([hrf(1) zeros(1,length(x)-1)],[hrf zeros(1,length(x)-1)]);
% a2=x*t_hrf;
% 
% plotme(sig_erp)
% plotme(sig_hrf')



Nt = 60;  % number of trials
Ne = 40;  % EEG resampling points
Nf = 240; % length of fMRI time points


repeat=1;
sige_all = [];
sigf_all = [];
sigt_all.fmri = [];
sigt_all.eeg = [];
behave_all = [];
neuro_act_all = [];
N_repe = 3;     % number of repeated times
x1 = rand(1,10000);  % event neuro active
for i3 = 1:3 
    x2(i3,:) = randn(1,10000); % rand neuro active 
end;

sigt=zeros(6,Nt);
%--------------change here make new simulation-----------
actSequ=randperm(Nt);
% save('data\signal_neuro.mat','actSequ');
% --------------------------------------------------------
N_stim = 40; % number of stim
neuro_act = zeros(6,N_stim); % stim
% ---------------------------------------------------
% below very good 
neuro_act(2,:) = 0.5*cos(2 *pi* x1(1:N_stim))+1;   % DMN 2
neuro_act(3,:) = 0.2*cos(2 *pi* x1(1:N_stim))+1;   % AN  3
neuro_act(4,:) = 0.4*cos(2 *pi* x1(1:N_stim))+1;   % ACC 4
neuro_act(5,:) = 0.6*cos(2 *pi* x1(1:N_stim))+1;   % LFC 5
neuro_act(6,:) = 0.7*cos(2 *pi* x1(1:N_stim))+1;   % RFC 6

% ------------ eeg sigt ---------
for i = 1:6
    if     i==1
        rand_neuro = x2(1,((repeat-1)*Nt+1):(repeat*Nt));
        sigt_eeg(i,:) = rand_neuro./max(abs(rand_neuro));
    elseif i==3
        rand_neuro = x2(2,((repeat-1)*Nt+1):(repeat*Nt));
        sigt_eeg(i,:) = rand_neuro./max(abs(rand_neuro)); % Temporal blind 3
    else
        sigt_eeg(i,actSequ(1:N_stim)) = 1;
        sigt_eeg(i,sigt_eeg(i,:)==1) = neuro_act(i,:);
    end;
end;


t_eeg(6,:)=zeros(size(sigt_eeg(1,:)));
t_eeg(1,[2:11,22:31,42:51])=[ones(1,30)];
t_eeg(2,[24:33])=[ones(1,10)];
t_eeg(3,:)=sigt_eeg(4,:);
t_eeg(4,:)=zeros(size(sigt_eeg(1,:)));
t_eeg(4,[9:13])=0.44*[ones(1,5)];
t_eeg(4,[29:33])=0.6*[ones(1,5)];
t_eeg(4,[49:53])=0.4*[ones(1,5)];

t_eeg(5,:)=[60:-1:1]/60;
t_eeg(6,:)=abs(sigt_eeg(1,:));
t_eeg(6,randi(60,[1,30]))=0;

t_eeg([1,5,2,6,3,4],:)=t_eeg;
sig_erp([1,2,6,5,4,3],:)=sig_erp;

figure
for i=1:6
    subplot(5,6,i)
        teeg1=sig_erp(i,:)'*t_eeg(i,:);
        imagesc(teeg1');
        teeg(i,:)=teeg1(:);
        xlabel('Trial time in number of trials')
    subplot(5,6,6+i)
        plot([1:10:400],sig_erp(i,:))
        xlabel('ERP time in msec')
    subplot(5,6,2*6+i)
        stem(t_eeg(i,:),'.')
        ylim([-0.5,1.5])
        xlabel('Trials')
    subplot(5,6,3*6+i)
        if hrf==1
            plot(sig_hrf(:,1))
        elseif hrf==0
            plot(sig_hrf(:,i))
        else 
            error ('HRF must be either 0 or 1 see help');
        end
        xlabel('Number of samples (Tr=2s)')
        axis tight;
    subplot(5,6,4*6+i)
        onset=[t_eeg(i,:); zeros(3,Nf/4)];
        onset=onset(:);
        if i==1
            onset([4*1:4*10+3,4*21:4*30+3,4*41:4*50+3],1)=[ones(1,4*30)];
        elseif i==5
            onset([4*23:4*32+3],1)=[ones(1,4*10)];  
        elseif i==6
            onset([4*8:4*12+3],1)=0.44*[ones(1,4*5)];  
            onset([4*28:4*32+3],1)=0.6*[ones(1,4*5)]; 
            onset([4*48:4*52+3],1)=0.4*[ones(1,4*5)]; 
        end  
        if erp==0
            temp(:,i)=onset;
        end
        if hrf==1
            tfmri(i,:)=((conv(onset,sig_hrf(:,1))));
        elseif hrf==0
            tfmri(i,:)=((conv(onset,sig_hrf(:,i))));
        else 
            error ('HRF must be either 0 or 1 see help');
        end
        plot(tfmri(i,:))
        xlabel('Number of samples (Tr=2s)')
        axis tight;
        ylim([-0.5,1.5])
end
        if erp==0
            clear teeg
            for i=1:6
                x1=[1:size(temp,1)]';
                x2=[1:0.25:size(temp,1)+0.75]';
                teeg(:,i)=interp1q(x1,temp(:,i),x2);
            end
        end
end
