addpath(genpath('/home/v1cchatz/Dropbox/matlab tools'))
cd('/home/v1cchatz/Dropbox/matlab tools/BBCB_code.tar/data')

load('ss.mat')
load('sa.mat')
load('data/miscdata')

mri=cor2mni(sa.cortex75K.vc);
[~,idx]=sort(mri(:,3));
sorty_mri=mri(idx,:);
biggest_z=mode(mri(:,3)); %I am performing the simulation on the biggest slice it can change or add more slices. 
[i_ind,~]=find(sorty_mri(:,3)==biggest_z);
cutslice=sorty_mri(min(i_ind):max(i_ind),:);


[~,id_un] = unique(cutslice(:,1:2),'rows','stable'); % Find unique values
dup_id=setdiff(1:size(cutslice,1),id_un); %Find replicated voxels and then replace them

cutslice(dup_id,1)=cutslice(dup_id,1)-1; % Decrease x of same voxels.

source_amp=zeros(size(sa.cortex75K.EEG_V_fem_normal,2),size(S,1));
for i=1:size(S,1)
    AS=imresize(reshape(S(i,:),100,100),[max(cutslice(:,1))-min(cutslice(:,1))+1,max(cutslice(:,2))-min(cutslice(:,2))+1] );
    for j=1:size(i_ind)
%         source(j,i)=AS(cutslice(j,1)-min(cutslice(:,1))+1,cutslice(j,2)-min(cutslice(:,2))+1);
          x(cutslice(j,1)-min(cutslice(:,1))+1,cutslice(j,2)-min(cutslice(:,2))+1)=1;
    end
    source_amp(i_ind+10000,i)=source(:,i); 
    EEG_field_pat(:, i) = sa.cortex75K.EEG_V_fem_normal*source_amp(:, i);
end

load('tools/cm17')
dataset_string='tests1'; %where to save
for i=1:size(S,1)
 ma = max(abs(EEG_field_pat(:, i)));
 allplots_head(sa, sa.EEG_elec2head*EEG_field_pat(:, i), [-ma ma], cm17, 'A.U.', ['figures/' dataset_string '/truth/dataset/EEG_pat1'], sa.EEG_locs_3D(:, 1:3));
end


