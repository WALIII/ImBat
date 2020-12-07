function ImBat_NextFlight(out_markov,FlightAlignedROI,ROI2use)

% Whats the next flight, and the time interval between it?
% ANY structure in the ROI activity?
clust_number = FlightAlignedROI.clust_number;
%id2use = find(out_markov.FlightIDVector= 2);
sekLen = [2 ];
size(sekLen)
id2use = strfind(out_markov.FlightIDVector,sekLen);;
% find the next flight, sorted by time:
%nextFlight = out_markov.FlightIDVector(id2use-(1+size(sekLen,2)));
nextFlight = out_markov.FlightIDVector(id2use+(size(sekLen,2)));

% get the ID of the clusted flight, but the ROI index:
% sort2use = out_markov.out_sort;
% nextFlightID = (sort2use(id2use));


%look at the flights 1-5;
figure();
hold on;
col = jet(5);
for i = 1:5;
id = find(nextFlight ==i);
id(id ==1) = []; % remove if first flight
idX{i} = id;% this will work on clusters that are not '1'
if size(id,2)>0
for ii = 1: size(id,2)
nextFlightTime(ii) = out_markov.FlightTimeVector(id(ii))-out_markov.FlightTimeVector(id(ii)-1);
end
histogram(nextFlightTime/120,100,'FaceColor',col(i,:),'BinWidth',5);
[a2,b2] = sort(nextFlightTime);
idX_sort{i} = b2;

end
clear nextFlightTime id
end
%TD = out_markov.FlightTimeVector(id2use+1)-out_markov.FlightTimeVector(id2use);

figure(); histogram(id2use);





% Plot calcium


ROIMAT = zscore(squeeze(FlightAlignedROI.C_raw(ROI2use,:,:)));


figure();
for i = 2;
    for ii = 1: size(idX{i},2);
    %a(ii) = find(FlightAlignedROI.cluster_idX == sort2use(idX{i}(ii)));
   a(ii) = idX{i}(ii);
    ROImat(:,ii) = ROIMAT(:,a(ii));
    end

ROImat2 = ROImat(:,idX_sort{i});
aaa2 = sum(ROImat2);
aaa2rmv2 = find(aaa2 ==0);
ROImat2(:,aaa2rmv2) = [];

aaa = sum(ROImat);
aaa2rmv = find(aaa ==0);
ROImat(:,aaa2rmv) = [];
clear aaa aaa2rmv
figure();
subplot(1,2,1)
   imagesc(ROImat',[-2 4]);
   title('sorted by time');
   subplot(1,2,2);
      imagesc(ROImat2',[-2 4]);
      title('sorted by lag');





end

