function trajm=ImBat_MCS_get_mic_pos(FLDnum)

% FLDnum is the session number
% micnr is the number of mics in the room
mic_type=[1,2,3,4,5,6;48,12,48,12,45,12];

stclmn=3;
pnmics=[filesep 'Users' filesep 'madeleinesnyder' filesep 'Desktop' filesep 'Audio_Calibration' filesep];
micnr=6;
trajm=[]; trajma=[];
micfilename_lst=dir(fullfile(pnmics, '*.trc'));
%for i=1:size(micfilename,2)
%    if micfilename(i).name,'Unnamed')
micnr=size(micfilename_lst,1);
for mic=1:micnr
    micfilename=micfilename_lst(mic);
    if ~isempty(micfilename)
        micfilename=[pnmics micfilename.name];
        micpos=textread(micfilename,'','headerlines',6);
        if size(micpos,2)>11
            if size(micpos,2)<10 && size(micpos,2)>6
                micmnr=2;
            elseif size(micpos,2)<7
                micmnr=1;
            else
                micmnr=3;
            end
            for micm=1:micmnr,
                if mic==6 | mic==4 | mic==2
                    markerpos=stclmn+3*(micm-1);%col(marker);
                    trajma(micm,1,mic)=micpos(end,markerpos);
                    trajma(micm,2,mic)=micpos(end,markerpos+1);
                    trajma(micm,3,mic)=micpos(end,markerpos+2);
                elseif mic==5 | mic==3 | mic==1
                    markerpos=stclmn+3*(micm-1);%col(marker);
                    trajma(micm,1,mic)=micpos(200,markerpos);
                    trajma(micm,2,mic)=micpos(200,markerpos+1);
                    trajma(micm,3,mic)=micpos(200,markerpos+2);
                else    
                    markerpos=stclmn+3*(micm-1);%col(marker);
                    trajma(micm,1,mic)=micpos(1,markerpos);
                    trajma(micm,2,mic)=micpos(1,markerpos+1);
                    trajma(micm,3,mic)=micpos(1,markerpos+2);
                end
            end
            for n=1:micmnr
                trajm(mic,n)=mean(trajma(:,2));
            end
        else
            disp('Not enough Mic data')
        end
    else
        disp('No Mic data')
    end
    trajm=[];
    for i=1:6
        trajm(i,:) = mean(trajma(:,:,i),1);
    end
end

% Plotting to figure out what micpost to use
idxs = [3,5,7];
figure(); hold on; scatter3(flightPaths34.pos(1,:,1),flightPaths34.pos(2,:,1),flightPaths34.pos(3,:,3),[],'r'); %scatter3(trajma(1,1,1)./1000,trajma(1,2,1)./1000,trajma(1,3,1)./1000,[],'b'); 
cmap=jet(12);
%for i=3:3
    idxs = [3,6,9];
    scatter3(micpos(:,idxs(1))./1000,micpos(:,idxs(2))./1000,micpos(:,idxs(3))./1000,[],cmap(i,:));
%end
end