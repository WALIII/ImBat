function ROI_Data = ImBat_RepairFlightData(ROI_Data);
% Repair flight data

n = 100; % minimum length of gap to repair:
hang1 = 3;%overhang

for day2use = 1:size(ROI_Data,2)
    
    disp(['Processing day ',num2str(day2use)]);
GG = ROI_Data{day2use}.Alignment.out.Location3(:,1);
Flights_Smoothed = ROI_Data{day2use}.Alignment.out.Location3(:,:);
Flights_Original = ROI_Data{day2use}.Alignment.out.Location(:,:);



GG2 = GG;
hold on;
% find stretches of zeros:
x = diff(smooth(GG,1));
transitions = diff([0; x == 0; 0]); %find where the array goes from non-zero to zero and vice versa
runstarts = find(transitions == 1);
runends = find(transitions == -1); %one past the end
runlengths = runends - runstarts;
%keep only those runs of length n or more:
runstarts(runlengths < n) = [];
runends(runlengths < n) = [];
%expand each run into a list indices:
indices = arrayfun(@(s, e) s:e-1, runstarts, runends, 'UniformOutput', false);
indices = [indices{:}];  %concatenate the list of indices into one vector
GG2(indices) = NaN; %replace the indices with 2
indices = [indices size(GG2,1)];
%plot(GG2);

% replace zeros:
Flights_new = ROI_Data{day2use}.Alignment.out.Location2(:,:);
Flights_new(indices,:) = Flights_Smoothed(indices,:);
% figure(); hold on; plot(Flights_new(:,1)); plot(GG);

%% do it again:
clear transitions runstarts runends indices GG GG2;
GG = Flights_new(:,1);
GG2 = GG;
GG3 = GG;
% find stretches of zeros:
x = diff(smooth(GG,1));
transitions = diff([0; x == 0; 0]); %find where the array goes from non-zero to zero and vice versa
runstarts = find(transitions == 1);
runends = find(transitions == -1); %one past the end
runlengths = runends - runstarts;
%keep only those runs of length n or more:
runstarts(runlengths < n) = [];
runends(runlengths < n) = [];
%expand each run into a list indices:
indices = arrayfun(@(s, e) s:e-1, runstarts, runends, 'UniformOutput', false);
indices = [indices{:}];  %concatenate the list of indices into one vector
GG2(indices) = NaN; %replace the indices with 2
indices = [indices size(GG2,1)];


Flights_fixed = Flights_new;
% now, get rid of them..

for i = 1: size(runstarts,1)
    % check which side is bigger:
    if runstarts(i) ==1;
        b = Flights_new(runends(i)+hang1,:);
        toUse = b;
        Flights_fixed(runstarts(i):runends(i)+hang1,:) = ones(size(runstarts(i):runends(i)+hang1,2),3).*toUse;
        
    elseif runends(i) > length(GG3)-hang1
        a = Flights_new(runstarts(i)-hang1,:);
        toUse = a;
        Flights_fixed(runstarts(i)-hang1:runends(i),:) = ones(size(runstarts(i)-hang1:runends(i),2),3).*toUse;
        
    else
        a = Flights_new(runstarts(i)-hang1,:);
        b = Flights_new(runends(i)+hang1,:);
        if a(:,1)>b(:,1); toUse = a;
        else
            toUse = b;
        end
        Flights_fixed(runstarts(i)-hang1:runends(i)+hang1,:) = ones(size(runstarts(i)-hang1:runends(i)+hang1,2),3).*toUse;
        
    end
end

for i = 1:3
Flights_fixed(:,i) = smooth(Flights_fixed(:,i),50);
end
figure();
hold on;
plot(Flights_Original(:,1),'k--');

plot(Flights_new(:,1),'b--');
plot(Flights_fixed(:,1),'r');

 legend('original','smoothed','repaired');
 title(['day: ',num2str(day2use)]);
 hold off
ROI_Data{day2use}.Alignment.out.Flights_Repaired = [];
 ROI_Data{day2use}.Alignment.out.Flights_Repaired(:,:) = Flights_fixed;
 
end