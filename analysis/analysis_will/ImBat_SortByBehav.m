
function  FlightAligned_Sorted = ImBat_SortByBehav(FlightAlignedROI)

% sort the calcium matrix by the top prefered peak for each cell
FlightAligned_Sorted = FlightAlignedROI;

% for all Flights
Plotting =0;


for i = 1:size(FlightAlignedROI,2)
    
    for iii = 1:size(FlightAlignedROI{i}.C,1)
        
        % Plot ROI:
        Ca1 = zscore(squeeze(FlightAlignedROI{i}.C_raw(iii,:,:)))';
        Ca1_d = zscore(squeeze(FlightAlignedROI{i}.C(iii,:,:)))';
        Ca1_s = zscore(squeeze(FlightAlignedROI{i}.S(iii,:,:)))';
        Fl1 = FlightAlignedROI{i}.ClustFlight;
        Fl1 = Fl1(50:550,:,:);
        
        
        % find index with and without zeros and remove them:
        mIdx = sum((Ca1),2);
        ind2rmv = find(mIdx==0); % remove non-zeros..
        ind2keep = find(mIdx~=0); % save this for later..
        
        % remove non-tracking from Fl and Ca
        Ca_raw_sort= Ca1; % save this for later...
        Ca_sort= Ca1_d; % save this for later...
        Ca_S_sort= Ca1_s; % save this for later...
        
        % remove non-tracking from Fl and Ca
        
        Ca1(ind2rmv,:) = []; % remove zeros
        Ca1_d(ind2rmv,:) = []; % remove zeros
        Ca1_s(ind2rmv,:) = []; % remove zeros
        
        Fl1(:,:,ind2rmv) = [];
        
        % sort based on peak:
        %figure(); plot(mean(Ca1))
%         [~,m1] = max(mean(Ca1(:,FlightAlignedROI{i}.ROI_ON-10:FlightAlignedROI{i}.ROI_ON+mean(FlightAlignedROI{i}.FlightLength)/120*30+60)));
%         m1 = m1+FlightAlignedROI{i}.ROI_ON-10;
%         [~,idx] = sort(Ca1(:,round(m1)),'descend');
 %               Gmean = mean(Fl1(:,:,idx(1:floor(size(Fl1,3)*.05))),3);

        %% sort based on cluster
try
l = linkage(Ca1(:,100:300), 'ward');
catch
    disp('NO DATA');
    return
end
c=cluster(l,'maxclust',3);
[aa,bb]=sort(c);
% bound2use = diff(aa);
%line2plot = find(bound2use==1);
idx = (bb); % the sort 

 Gmean = mean(Fl1(:,:,(find(c==1))),3);
    clear aa bb c l
    
        
        % now get the flights that correspond to this, and calculate mean euclid distance from all flights
        
        
        for ii = 1:size(Fl1,3)
            % Euclidian distance
            FL(ii) = sum(sqrt(sum((Gmean' - Fl1(:,:,ii)') .^ 2)));
        end
        % Plotting
        
        % Sort based on flights:
        try
            [~,idx_fl] = sort(FL);
            
            % index overide:
            idx_fl = idx;
            Ca_raw_sort(ind2keep,:) = Ca1(idx_fl,:);
            Ca_sort(ind2keep,:) = Ca1_d(idx_fl,:);
            Ca_S_sort(ind2keep,:) = Ca1_s(idx_fl,:);
            
            
            FlightAligned_Sorted{i}.C_raw(iii,:,:) =  squeeze(Ca_raw_sort)';
            FlightAligned_Sorted{i}.C(iii,:,:) =  squeeze(Ca_sort)';
            FlightAligned_Sorted{i}.S(iii,:,:) =  squeeze(Ca_S_sort)';
            
            if Plotting==1;
                Plot_stuff(Ca1,idx,Gmean,Fl1,Ca_raw_sort,idx_fl);
                pause();
                clf('reset');
                close all
            end
            
            
        catch
            disp(' no calcium left...');
        end
        
        
        clear Ca_sort Ca_S_sort Ca_raw_sort S_raw_sort ind2keep idx_fl FL Gmean Fl1 idx Ca1           Ca1_d Ca1_s Fl1 mIdx ind2rmv ind2keep
        
    end
end
end

function Plot_stuff(Ca1,idx,Gmean,Fl1,Ca_raw_sort,idx_fl);


%%  Plot stuff:

figure(); imagesc(Ca1(idx,:))

figure();
hold on;
plot(Gmean(:,1),Gmean(:,2),'LineWidth',3,'color','m');
for ii = 1:size(Fl1,3)
    if ii<40;
        plot(Fl1(:,1,idx(ii)),Fl1(:,2,idx(ii)),'r');
    else
        plot(Fl1(:,1,idx(ii)),Fl1(:,2,idx(ii)),'k');
    end
end


figure(); imagesc(Ca1);
title('sorted by time');

figure(); imagesc(Ca_raw_sort)

figure();
hold on;
plot(Gmean(:,1),Gmean(:,2),'LineWidth',3,'color','m');

for ii = 1:size(Fl1,3)
    if ii<40;
        plot(Fl1(:,1,idx_fl(ii)),Fl1(:,2,idx_fl(ii)),'r');
    else
        plot(Fl1(:,1,idx_fl(ii)),Fl1(:,2,idx_fl(ii)),'k');
    end
end

end


%end
