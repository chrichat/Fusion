function [sfmri,seeg]=Simulate_space()
%=========================================================

sty={'eeg','neuro','fmri'};    % neuro fmri eeg
figure
load('ss_ho.mat')


for i=1:6
    for j=1:3
        Si=i;
        style=sty{j};
        switch lower(style)
            case 'neuro'
                %========================================================
                %                   Neuro 
                %========================================================
                ax1=subplot(3,6,2*6-6+i);
                [p,ellipse]=lei_Simulate_Space_myphantom(70,8);
                p(p==254)=0.1;
                tempimg=reshape(source2(:,i),70,70)+p;
                if i==1
                    tempimg(p==255)=1.9;
                else
                    tempimg(p==255)=1.2;
                end
                tempimg(abs(tempimg)<0.05)=0;
                imagesc(tempimg)
                xlabel('Areas of activation')
                col_scale =colormap(jet(256));
                col_scale(1,:)=[0 0 0];% 
                col_scale(256,:)=[1 1 1];% scalp
                colormap(ax1,col_scale);
                axis equal
                axis off
                
            case 'fmri'
                %========================================================
                %                   fMRI
                %========================================================
                ax2=subplot(3,6,3*6-6+i);
                [p,ellipse]=lei_Simulate_Space_myphantom(70,8);
                p(p==254)=0.1;
                tempimg=reshape(source2(:,i),70,70)+p;
                if i==1
                    tempimg(p==255)=1.9;
                else
                    tempimg(p==255)=1.2;
                end
                tempimg(abs(tempimg)<0.05)=0;
                imagesc(tempimg)
                if i==1
                    tempimg(tempimg==1.9)=0;
                else
                    tempimg(tempimg==1.2)=0;
                end
                sfmri(i,:)=tempimg(:);
                col_scale =colormap(jet(256));
                col_scale(1,:)=[0 0 0];% 
                col_scale(256,:)=[1 1 1];% scalp
                colormap(ax2,col_scale);
                xlabel('Spatial maps of activation')
                axis equal
                axis off
                
            case 'eeg'
                %========================================================
                %                   EEG
                %========================================================
                % lead field matrix
                ax31=subplot(3,6,6*1-6+i);
                load('data/leadfield.mat');
                % 2452 voxels is brain areas in one slice (70*70=4900 voxels) 
                load('data/brainMask.mat');
                % lead 6 sources index
                load('data/space.mat');
                inside=find(brainMask==1);
%                 space=space(inside,Si);
                source1=S';
                tempimg=flipdim(imresize(reshape(source1(:,i),100,100),[58,58]),1);
                source2(:,i)=reshape(padarray(tempimg,[6,6]),[],1);
                space=source2(inside,Si);
                topo=leadfield*space;
                seeg(i,:)=topo(:);
                topoplot(topo,'data/GSN128.loc',  'maplimits', ...
                    [min(min(topo)), max(max(topo))],'style','both','shrink','on');   
                xlabel('Topoplot of the electrodes')
        end
    end
end
fprintf('\nDone\n')