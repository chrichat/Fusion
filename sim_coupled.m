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
teeg(958:end,:)=0;

filepath = cd;
columnLabels = {' Sour. 1', ' Sour. 2',' Sour. 3',' Sour. 4',' Sour. 5',' Sour. 6'};
autocor_eegtemp=abs(corr(teeg,teeg));
filename='autote.tex';
file     = fullfile(filepath, filename);
matrix2latex(autocor_eegtemp, file, 'rowLabels', columnLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 
autocor_eegspace=abs(corr(seeg',seeg'));
filename='autose.tex';
file     = fullfile(filepath, filename);
matrix2latex(autocor_eegspace, file, 'rowLabels', columnLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 
autocor_fmritemp=abs(corr(tfmri',tfmri'));
filename='autotf.tex';
file     = fullfile(filepath, filename);
matrix2latex(autocor_fmritemp, file, 'rowLabels', columnLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 
autocor_fmrispace=abs(corr(sfmri',sfmri'));
filename='autosf.tex';
file     = fullfile(filepath, filename);
matrix2latex(autocor_fmrispace, file, 'rowLabels', columnLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 



load('data/signal_hrf.mat')
hrf=sig_hrf(:,1)';

A1=abs(randn(6,6));
A2=A1;


snr1=[0.1 0.2 0.5 1 1.4 2];
snr2=[  2  5 10 12 14 180];
rowLabels = {'Es. Sour. 1', 'Es. Sour. 2','Es. Sour. 3','Es. Sour. 4','Es. Sour. 5','Es. Sour. 1'};

for i=1:6
    close all
    teeg(958:end,:)=0;
    Xeeg=seeg([1,2,3,4,5,6],:)'*kr(A1,teeg(:,[1,2,3,4,5,6]))'; % Source 5 is not selected - blind for EEG - unshare
    Xeeg=Xeeg + mean(std(Xeeg))*sqrt(snr1(i))*randn(size(Xeeg)); % 
    % Xeeg2=awgn(Xeeg,1,'measured');
%     ww=jader(Xeeg,6);
%     topoplotme(pinv(ww)')
%     temp2=sprintf('/home/v1cchatz/Dropbox/sims_fusi/seeg%d.png',i);
%     saveas(gcf,temp2)
%     temp=abs((ww*Xeeg));
%     temp=reshape(temp,size(temp,1),[],size(A1,2));
%     plotme(abs(mean(temp,3)))
%     temp2=sprintf('/home/v1cchatz/Dropbox/sims_fusi/teeg%d.png',i);
%     saveas(gcf,temp2)
%     cor_ica_eegtemp=abs(corr(mean(temp,3)',teeg));
%     temp2=sprintf('/home/v1cchatz/Dropbox/sims_fusi/teeg%d.tex',i);
%     matrix2latex(cor_ica_eegtemp, temp2, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 
%     cor_ica_eegspace=abs(corr(pinv(ww),seeg'));
%     temp2=sprintf('/home/v1cchatz/Dropbox/sims_fusi/seeg%d.tex',i);
%     matrix2latex(cor_ica_eegspace, temp2, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 
    Teeg=reshape(Xeeg,size(Xeeg,1),size(teeg,1),[]);

    
    Xfmri=sfmri([1,2,3,4,5,6],:)'*kr(A2,tfmri([1,2,3,4,5,6],:)')'; % Source 5 is not selected - blind for EEG - unshare
    % Xfmri=awgn(Xfmri,0.1,'measured');
    Xfmri=Xfmri + mean(std(Xfmri))*sqrt(snr2(i))*randn(size(Xfmri)); % 
    Tfmri=reshape(Xfmri,size(Xfmri,1),size(tfmri,2),[]);
%     ww=jader(Xfmri',6);
%     plotme(ww*Xfmri',70,70)
%     temp2=sprintf('/home/v1cchatz/Dropbox/sims_fusi/sfmri%d.png',i);
%     saveas(gcf,temp2)
%     temp=reshape(pinv(ww)',size(ww,1),[],size(A1,2));
%     plotme(abs(mean(temp,3)))
%     temp2=sprintf('/home/v1cchatz/Dropbox/sims_fusi/tfmri%d.png',i);
%     saveas(gcf,temp2)
%     cor_ica_fmrispac=abs(corr(mean(temp,3)',tfmri'));
%     temp2=sprintf('/home/v1cchatz/Dropbox/sfmri%d.tex',i);
%     matrix2latex(cor_ica_fmrispac, temp2, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 
%     cor_ica_fmritemp=abs(corr((ww*Xfmri')',sfmri'));
%     temp2=sprintf('/home/v1cchatz/Dropbox/tfmri%d.tex',i);
%     matrix2latex(cor_ica_fmritemp, temp2, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 

% [U1,output] = cpd_gevd(Teeg,R1);
% [U2,output] = cpd_gevd(Tfmri,R2);
    
% [Vica,Sica,sv]=mlsvd(Tfmri,[size(Tfmri,2)*size(Tfmri,3),size(Tfmri,2),size(Tfmri,3)]);

% [U,S]=mlsvd_rsi(T2,[5*R2,5*R2,size(T2,3)]);

[ sol,sol2,sol3] = coupled_tensors( Teeg,Tfmri,R1,R2,Nt,'cpd',hrf );
% U{1,1}=sol.factors.D;
% U{1,2}=sol.factors.E;
% U{1,3}=sol.factors.F;
% U=cellfun(@(u,v)v*u,U,Vica,'UniformOutput',false); %%uncompressing results
% 
% U2{1,1}=sol3.factors.D;
% U2{1,2}=sol3.factors.E;
% U2{1,3}=sol3.factors.F;
% U2=cellfun(@(u,v)v*u,U2,Vica,'UniformOutput',false); %%uncompressing results



temp3=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\sol_%d.mat',i);
save(temp3,'sol3')
topoplotme(sol2.factors.A')
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\se_uncoup%d.png',i);
saveas(gcf,temp2)
plotme(abs(sol2.factors.B'))
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\te_uncoup%d.png',i);
saveas(gcf,temp2)
cor_unc_eegtemp=abs(corr(sol2.factors.B,teeg));
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\te_uncoup%d.tex',i);
matrix2latex(cor_unc_eegtemp, temp2, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 
cor_unc_eegspace=abs(corr(sol2.factors.A,seeg'));
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\se_uncoup%d.tex',i);
matrix2latex(cor_unc_eegspace, temp2, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 


temp=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\sol2_%d.mat',i);
save(temp,'sol2')
plotme(sol3.factors.D',70,70)
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\sf_uncoup%d.png',i);
saveas(gcf,temp2)
plotme(sol3.factors.E')
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\tf_uncoup%d.png',i);
saveas(gcf,temp2)
cor_uncup_fmritemp=abs(corr(sol3.factors.E,tfmri'));
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\tf_uncoup%d.tex',i);
matrix2latex(cor_uncup_fmritemp, temp2, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 
cor_uncoup_fmrispace=abs(corr(sol3.factors.D,sfmri'));
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\sf_uncoup%d.tex',i);
matrix2latex(cor_uncoup_fmrispace, temp2, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 


temp=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\sol3_%d.mat',i);
save(temp,'sol3')
topoplotme(sol.factors.A')
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\se_coup%d.png',i);
saveas(gcf,temp2)
plotme(abs(sol.factors.B'))
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\te_coup%d.png',i);
saveas(gcf,temp2)
cor_c_eegtemp=abs(corr(sol.factors.B,teeg));
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\te_coup%d.tex',i);
matrix2latex(cor_c_eegtemp, temp2, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 
cor_c_eegspace=abs(corr(sol.factors.A,seeg'));
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\se_coup%d.tex',i);
matrix2latex(cor_c_eegspace, temp2, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 
plotme(sol.factors.D',70,70)
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\sf_coup%d.png',i);
saveas(gcf,temp2)
plotme(sol.factors.E')
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\tf_coup%d.png',i);
saveas(gcf,temp2)
cor_cup_fmritemp=abs(corr(sol.factors.E,tfmri'));
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\tf_coup%d.tex',i);
matrix2latex(cor_cup_fmritemp, temp2, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 
cor_coup_fmrispace=abs(corr(sol.factors.D,sfmri'));
temp2=sprintf('C:\\Users\\Christos\\Dropbox\\sims_fusi\\sf_coup%d.tex',i);
matrix2latex(cor_coup_fmrispace, temp2, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 

end