function  [out] =  ImBat_SegTrajectories(Location,Location_times)
% ImBat_SegTrajectories

% Nick Dotson
% Modified by WAL3

% 06/06/19



FS = 120;

 
% event_ttls = AnalogSignals(:,2); %from motion data
% [R,LT,UT,LL,UL] = risetime(event_ttls,FS); %find times of ttl pulses in SECONDS
%for now, clip event ttls after the number of trials
 
%might want to use raw c3d to keep track of each marker (this could be
%useful for when not all 3 markers are present, which is often the case
%when the animal is at a wall)
% mx = Markers(:,1,1);
% my = Markers(:,1,2);
% mz = Markers(:,1,3);
  
mx = Location(:,1);
my = Location(:,2);
mz = Location(:,3);

 
mxi = []; myi = []; mzi = [];
%find glitches using diff
% % mxi = unique([find(abs(diff(mx)) > 50); find(mx == 0)]);
% % myi = unique([find(abs(diff(my)) > 50); find(my == 0)]);
% % mzi = unique([find(abs(diff(mz)) > 50); find(mz == 0)]);
%may or may not be necessary
 
 
%set zeros to nan
mx(find(mx == 0)) = nan;
my(find(my == 0)) = nan;
mz(find(mz == 0)) = nan;
 
mx = fillmissing(mx,'nearest');
my = fillmissing(my,'nearest');
mz = fillmissing(mz,'nearest');
 
%threshold based on speed
Vx = gradient(mx, 1/FS);
Vy = gradient(my, 1/FS);
Vz = gradient(mz, 1/FS);
speed = sqrt(Vx.^2 + Vy.^2 + Vz.^2)/1000; %in m/s
 
nonflying = find(speed < 1.5);
toofast = find(speed > 30);
%set low speed points to zero
mx([nonflying; toofast]) = nan;
my([nonflying; toofast]) = nan;
mz([nonflying; toofast]) = nan;
 
trajectories_continuous(1,:) = mx;
trajectories_continuous(2,:) = my;
trajectories_continuous(3,:) = mz;
 
figure
plot3(mx,my,mz,'LineWidth',2,'Color','k')
 
%splice out individual flights
batspeed = speed;
batspeed(nonflying) = nan;
batspeed(toofast) = nan;
 
bflying=~isnan(batspeed);
 
%for each data point, sum up the next 1s of data
allsums = [];
for bf = 1 : size(bflying,1)-FS
    allsums(bf) = sum(bflying(bf:bf+90));
end
 
[R,rLT,rUT,rLL,rUL] = risetime(allsums);
[F,fLT,fUT,fLL,fUL] = falltime(allsums);
flight_starts = round(rLT);
flight_ends = round(fLT);
 
%cut out flights and save
for nf = 1 : size(R,2)
    hold on
    plot3(mx(flight_starts(nf):flight_ends(nf)),my(flight_starts(nf):flight_ends(nf)),mz(flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color','r')
    hold on
    
    fstartxyz(nf,1) = round(nanmean(mx(flight_starts(nf):flight_starts(nf)+90)));
    fstartxyz(nf,2) = round(nanmean(my(flight_starts(nf):flight_starts(nf)+90)));
    fstartxyz(nf,3) = round(nanmean(mz(flight_starts(nf):flight_starts(nf)+90)));
    
    fendxyz(nf,1) = round(nanmean(mx(flight_ends(nf):flight_ends(nf)+90)));
    fendxyz(nf,2) = round(nanmean(my(flight_ends(nf):flight_ends(nf)+90)));
    fendxyz(nf,3) = round(nanmean(mz(flight_ends(nf):flight_ends(nf)+90)));
    
    scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),1000,'g','filled')
    hold on
    scatter3(fendxyz(nf,1),fendxyz(nf,2),fendxyz(nf,3),1000,'b','filled')
%    pause
end
 
%determine the xyz start and stop positions
out.trajectories_continuous = trajectories_continuous;
out.trajectories_continuous = trajectories_continuous;
out.flight_starts_indx = flight_starts; 
out.flight_ends_indx = flight_ends; 
out.flight_starts_times = Location_times(flight_starts); 
out.flight_ends_times = Location_times(flight_ends); 
out.fstartxyz = fstartxyz; 
out.fendxyz = fendxyz;
% cd(daydir)
 
nclusters = 4;
 
jj = jet;
figure
rng(2)
kstart = kmeans(fstartxyz,nclusters);
for nf = 1 : size(R,2)
    plot3(mx(flight_starts(nf):flight_ends(nf)),my(flight_starts(nf):flight_ends(nf)),mz(flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color',jj(kstart(nf)*4,:))
    hold on
    scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),500,jj(kstart(nf)*4,:),'filled')
    hold on
end
 
jj = jet;
figure
rng(2)
kend = kmeans(fendxyz,nclusters);
for nf = 1 : size(R,2)
    plot3(mx(flight_starts(nf):flight_ends(nf)),my(flight_starts(nf):flight_ends(nf)),mz(flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color',jj(kend(nf)*4,:))
    hold on
    scatter3(fendxyz(nf,1),fendxyz(nf,2),fendxyz(nf,3),500,jj(kend(nf)*4,:),'filled')
    hold on
end
 
 
 
%find pairs of start and endpoints with a high number of flights
rng(2)
kstart = kmeans(fstartxyz,nclusters);
rng(2)
kend = kmeans(fendxyz,nclusters);
nflights = [];
allflights = [];
npair = [];
ncounter = 0;
for ks = 1 : nclusters
    for ke = 1 : nclusters
        ncounter = ncounter + 1;
        npair(ncounter,:) = [ks ke];
        nflights(ncounter) = size(intersect(find(kstart == ks),find(kend == ke)),1)
        allflights{ncounter} = intersect(find(kstart == ks),find(kend == ke))';
    end
end
 
[~, ssf] = sort(nflights,'descend');
figure
jj = jet;
for traj = 1 : 5
    subplot(1,5,traj)
    for nf = allflights{ssf(traj)}
        plot3(mx(flight_starts(nf):flight_ends(nf)),my(flight_starts(nf):flight_ends(nf)),mz(flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color',jj(traj*10,:))
        hold on
        scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),100,'r','filled')
        hold on
        scatter3(fendxyz(nf,1),fendxyz(nf,2),fendxyz(nf,3),100,'k','filled')
        hold on
    end
end
 
figure
jj = jet;
for traj = 1 : 5
    for nf = allflights{ssf(traj)}
        plot3(mx(flight_starts(nf):flight_ends(nf)),my(flight_starts(nf):flight_ends(nf)),mz(flight_starts(nf):flight_ends(nf)),'LineWidth',1,'Color',jj(traj*10,:))
        hold on
        scatter3(fstartxyz(nf,1),fstartxyz(nf,2),fstartxyz(nf,3),100,'r','filled')
        hold on
        scatter3(fendxyz(nf,1),fendxyz(nf,2),fendxyz(nf,3),100,'k','filled')
        hold on
    end
end
