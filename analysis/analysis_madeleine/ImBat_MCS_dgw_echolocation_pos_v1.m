clear all;
close all;
validnames={'Boo';'HalfPint';'Pi';'Polyphemus';'SeventyNine';'Stinkbaer';'Tweek'};
implantbats={'Boo';'Pi';'Tweek'};
ploth=0;%plot results
fs=192e3;%sampling rate of recs
fr=120;%frame rate of flight room vids changed MCS June 21
thresh=0.0013;%voltage threshold for echolocation clicks for earthwork mics
thresh2=0.005;%voltage threshold for echolocation clicks for knowles mics
calldist=0.015;%minimal distance between two clicks (in seconds)
errorval=0.005;%added distance in seconds to allow for slight peak detection errors
pn='H:\Navigation\alldat\neuro\BinnedClustersFlight\';
behavepn=['H:\Navigation\alldat\behave\LabNavextenddata\'];
pnbehave='H:\Navigation\alldat\behave\LabNavextend_done\';
soundpn='H:\Navigation\alldat\Labnavextendsounds\';
pnfl=['extracted' filesep 'Gen_200305_fly-1_track.mat'];
flinfo='Processed_Flight_';
soundpn = pwd;

soundflds=dir([soundpn '*20*']);

for soundfldnr=1:size(soundflds,1)%session number
    fldnum=soundflds(soundfldnr).name;
    disp(['Folder ' num2str(soundfldnr) ' of ' num2str(size(soundflds,1))])
    disp(fldnum)
    recmth=fldnum(1:2);
    recdy=fldnum(3:4);
    recy=fldnum(5:8);
    trajm=[];
    trajm=get_mic_pos(fldnum);%gets 3D position of each mic
    if size(trajm,1)==4
        trajm(5,:)=[0 2.5908e3 1.4732e3];%knowles mic on left side
        trajm(6,:)=[0 -2.5908e3 1.4732e3];%knowles mic on right side
        for bn=1:size(validnames,1)%bat number
            batname=validnames{bn};
            flch=dir([pnfl fldnum '\' flinfo batname '.mat']);
            if size(flch,1)~=0
                %load flight data
                load(pnfl)
                go=0;
                imp=0;
                if str2num(recy)==2017
                    load([pnbehave 'LabNavextend_' 'Boo' '.mat'])
                    dats=datestr(result(1,[2 end]),23);
                    dval=dats(2,:);
                    recmbat=str2num(dval(1:2));
                    dval=dats(2,:);
                    recdbat=str2num(dval(4:5));
                    if str2num(recmth)<recmbat
                        go=1;
                    elseif str2num(recmth)==recmbat
                        if str2num(recdy)<=recdbat
                            go=1;
                        end
                    end
                    if str2num(recmth)==12 && strcmp(batname,'Pi')
                        go=1;
                        imp=1;
                    end
                elseif str2num(recy)>2017
                    batnamecmp={batname};
                    l = cellfun(@(c)strcmp(c,implantbats),batnamecmp,'UniformOutput',false);
                    batv=find(l{1}==1);
                    if batv>0
                        go=1;
                        imp=1;
                    end
                end
                if go==1
                    wavflag=0;
                    %check whether recs are wav or mat files
                    soundfls=dir([soundpn fldnum '\LabNavextend_' batname '_' recy recmth recdy 'T*_tn*.mat']);
                    if size(soundfls,1)~=0
                        wavflag=1;
                    else
                        soundfls=dir([soundpn fldnum '\LabNavextend_' batname '_' recy recmth recdy 'T*_tn*.wav']);
                        if size(soundfls,1)~=0
                            wavflag=2;
                        end
                    end
                    if wavflag>0
                        disp(batname)
                        disp('...')
                        switchchan=0;
                        if str2num(recy)==2017
                            if str2num(recmth)<=7 && str2num(recdy)<=17
                                switchchan=1;
                            end
                        end
                        %load behavioral information
                        if imp==0
                            load([pnbehave 'LabNavextend_' batname '.mat'])
                        else
                            load([behavepn 'LabNavextend_' batname '.mat'])
                        end
                        trindx=strmatch([recmth recdy recy],datestr(result(1,:),'mmddyyyy'));
                        expltr=find(result(4,trindx));
                        latencies=result(end,expltr);
                        triallength=round(size(smtraj,3)/fr);
                        addwait=2;
                        if strcmp(batname,'Boo')
                            if str2num(recy)>2017
                                addwait=3;
                            end
                        end
                        %ensure number behavioral trials match number of
                        %flight trials
                        if length(brightness_vals)==size(smtraj,1)
                            %get trials according to intensity of light cue
                            for inten=1:2
                                if inten==1
                                    indxtr=find(brightness_vals<80);
                                    echopos_xl=NaN(length(indxtr),size(trajm,1),1e3);
                                else
                                    indxtr=find(brightness_vals>=80);
                                    echopos_xh=NaN(length(indxtr),size(trajm,1),1e3);
                                end
                                for soundfln=1:length(indxtr)%trial number
                                    tcorrv=1;
                                    tcorrvt=1;
                                    soundflnr=indxtr(soundfln);
                                    %correct for real trial start
                                    if soundflnr>1
                                        if latencies(soundflnr)<triallength-addwait
                                            tcorrv=floor((triallength-(latencies(soundflnr)+addwait))*fr)-1;
                                            tcorrvt=floor(triallength-(latencies(soundflnr)+addwait))-1;
                                        end
                                    end
                                    if tcorrv<1
                                        tcorrv=1;
                                    end
                                    if tcorrvt<1
                                        tcorrvt=1;
                                    end
                                    startfl=tcorrv;
                                    flight=squeeze(smtraj(soundflnr,:,startfl:end));
                                    if wavflag==2
                                        [y,fs]=audioread([soundpn fldnum '\' soundfls(soundflnr).name]);
                                    else
                                        load([soundpn fldnum '\' soundfls(soundflnr).name])
                                        y=recbuf;
                                    end
                                    %initially recchan 5 and 6 were switched,
                                    %correct for this
                                    if switchchan==1
                                        yinterim=y;
                                        yinterim(:,5)=y(:,6);
                                        yinterim(:,6)=y(:,5);
                                        y=yinterim;
                                    end
                                    mbdis=[];
                                    for micnr=1:size(y,2)-1
                                        %calculate the bat's 3D distance to each
                                        %mic for each vid frame
                                        mbdis(micnr,:)=get_3d_dist(trajm(micnr,:),flight);
                                        if ploth==1
                                            figure(micnr),subplot(4,1,1);
                                            tpl=linspace(0,length(y)/fs,length(y));
                                            plot(tpl,y(:,micnr))
                                        end
                                    end
                                    startsnd=round(tcorrvt*fs);
                                    y=y(startsnd:end,:);
                                    tpl=linspace(0,length(y)/fs,length(y));
                                    allindx=NaN(size(y,2)-1,1000);
                                    allpks=NaN(size(y,2)-1,1000);
                                    allpos=NaN(size(y,2)-1,1000);
                                    for chan=1:size(y,2)-1
                                        %correct mic's frequency response
                                        %if chan<5
                                        %load(['H:\Navigation\FLR_measurements\FLR_Mic_CompIR\Comp_IR_earthworks_' num2str(chan) '.mat'])
                                        %else
                                        %load(['H:\Navigation\FLR_measurements\FLR_Mic_CompIR\Comp_IR_knowles_' num2str(chan) '.mat'])
                                        %end
                                        % convy=filter(irc,1,y(:,chan));
                                        % bandpass signal for cleaner extraction
                                        [b,a]=butter(8,2*[10e3 40e3]./fs,'bandpass');
                                        convy=filter(b,a,y(:,chan));
                                        if ploth==1
                                            figure(chan),subplot(4,1,2);
                                            plot(tpl,convy);hold on;
                                        end
                                        indxc=[];
                                        indxex=[];
                                        if switchchan==1
                                            [pks,indxc]=findpeaks(convy,'MinPeakHeight',thresh2);
                                        else
                                            if chan<5
                                                [pks,indxc]=findpeaks(convy,'MinPeakHeight',thresh);
                                            else
                                                [pks,indxc]=findpeaks(convy,'MinPeakHeight',thresh2);
                                            end
                                        end
                                        %finds peaks in audio and then removes
                                        %peaks to close to eachother
                                        if ~isempty(indxc)
                                            for findx=1:length(indxc)-1
                                                diffvec=indxc(findx+1:end)-indxc(findx);
                                                newindx=find(diffvec<calldist*fs);
                                                if ~isempty(newindx)
                                                    pkvec=pks([(findx) (newindx+findx)']);
                                                    pkind=find(pkvec~=max(pkvec));
                                                    indxc(pkind+findx-1)= NaN;
                                                end
                                            end
                                        end
                                        %check again for too close peaks
                                        indxex=find(~isnan(indxc));
                                        if ~isempty(indxex)
                                            chvec=indxc(indxex);
                                            for findx=1:length(indxex)-1
                                                if findx<length(indxex)
                                                    diffvec=indxc(indxex(findx+1:end))-indxc(indxex(findx));
                                                    newindx=find(diffvec<calldist*fs);
                                                    if ~isempty(newindx)
                                                        pkvec=pks([indxex(findx) indxex(newindx+findx)']);
                                                        pkind=find(pkvec~=max(pkvec));
                                                        indxex(pkind+findx-1)= NaN;
                                                        indxex=indxex(find(~isnan(indxex)));
                                                    end
                                                end
                                            end
                                        end
                                        if ploth==1
                                            figure(chan),subplot(4,1,2);
                                            plot(indxc(indxex)./fs,pks(indxex),'r*');hold off;
                                            %set(gca,'ylim',[-1.3 1.3])
                                        end
                                        allindx(chan,1:length(indxex))=indxc(indxex);
                                        allpks(chan,1:length(indxex))=pks(indxex);
                                    end
                                    count=0;
                                    %check for each vid frame=bat position if at
                                    %same time point (with added delay depending on
                                    %bat's distance to mic) echolocation signal
                                    for tpt=1:size(flight,2)
                                        singpks=zeros(1,size(y,2)-1);
                                        for micnr=1:size(y,2)-1
                                            if tpt==1
                                                indx=find(allindx(micnr,:)>0 & allindx(micnr,:)<=(tpt/fr+mbdis(micnr,tpt)+errorval)*fs);
                                            else
                                                indx=find(allindx(micnr,:)>((tpt-1)/fr+mbdis(micnr,tpt-1)-errorval)*fs & allindx(micnr,:)<=(tpt/fr+mbdis(micnr,tpt)+errorval)*fs);
                                            end
                                            if ~isempty(indx)
                                                allpos(micnr,indx)=flight(1,tpt);
                                            end
                                        end
                                    end
                                    if ploth==1
                                        for micnr=1:size(y,2)-1
                                            figure(micnr),subplot(4,1,3);
                                            plot(squeeze(allpos(micnr,:)),1,'*');
                                            set(gca,'xlim',[0 6e3])
                                        end
                                        pause;
                                    end
                                    if inten==1
                                        echopos_xl(soundfln,1:size(allpos,1),1:size(allpos,2))=allpos;
                                    else
                                        echopos_xh(soundfln,1:size(allpos,1),1:size(allpos,2))=allpos;
                                    end
                                end
                            end
                            if exist(['H:\Navigation\Sound_Analyzed\' fldnum])==0
                                mkdir(['H:\Navigation\Sound_Analyzed\' fldnum])
                            end
                            if imp==0
                                save(['H:\Navigation\Sound_Analyzed\' fldnum '\Echopos_' batname '_behave.mat'],'echopos_x')
                            else
                                save(['H:\Navigation\Sound_Analyzed\' fldnum '\Echopos_' batname '_implant.mat'],'echopos_x')
                            end
                            if ploth==1
                                for micnr=1:6
                                    test=squeeze(echopos_x(:,micnr,:));
                                    mn=nanmean(test);
                                    figure(micnr);subplot(4,1,4);hist(mn,0:100:6e3)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
