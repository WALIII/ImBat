function trajm=get_mic_pos(FLDnum)

stclmn=3;
pnmics=['H:\Navigation\alldat\flight\' FLDnum '\Generated_C3D_files\'];
micnr=4;
trajm=[];
for mic=1:micnr
    micfilename=dir([pnmics 'Trimmed_Mic' num2str(mic) '_*.trc']);
    if ~isempty(micfilename)
        micfilename=[pnmics micfilename(1).name];
        micpos=textread(micfilename,'','headerlines',6);
        if size(micpos,2)>11
            trajma=[];
            if size(micpos,2)<10 && size(micpos,2)>6
                micmnr=2;
            elseif size(micpos,2)<7
                micmnr=1;
            else
                micmnr=3;
            end
            for micm=1:micmnr,
                markerpos=stclmn+3*(micm-1);%col(marker);
                trajma(micm,1)=micpos(1,markerpos);
                trajma(micm,2)=micpos(1,markerpos+1);
                trajma(micm,3)=micpos(1,markerpos+2);
            end
            for n=1:micmnr
                trajm(mic,n)=mean(trajma(:,n));
            end
        else
            disp('Not enough Mic data')
        end
    else
        disp('No Mic data')
    end
end
