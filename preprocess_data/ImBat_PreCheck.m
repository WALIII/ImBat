function [Extraction_Metadata] = ImBat_PreCheck(ROI_Data);

% Check to see that following are within expectations:
% 1. Check if all flight and video  end times are the same length
% 2. extraction dates are recent
% 3. Frame Rates are correct



%%  1. Check if all flight and video  end times are the same length
%% 2. extraction dates are recent


for day = 1:size(ROI_Data,2)
extraction_date{day,1} = ROI_Data{1,day}.date;
extraction_date{day,2} = ROI_Data{1,day}.ROIs.results.metadata.Extraction_date; 


A = ROI_Data{1,day}.Alignment.out.Location_time(end);
B = ROI_Data{1,day}.Alignment.out.video_times(end);

C(:,day) = A-B;
extraction_date{day,3} = C(:,day)*.33;
extraction_date{day,3} = C(:,day)*.33;
extraction_date{day,4} = ROI_Data{day}.ROIs.results.Fs;

% 
end
figure();
plot(C*.33);
title('Ending offsets')
xlabel('Days');
ylabel('Offset (seconds)');


X = find(abs(C)>120);
if size(X,2) ==0;
        disp(' ALL CLEAR: No bad dates')
        Bad_Dates = [];
else
for ii = 1:size(X);
    disp(['WARNING: Failed alignment: ',ROI_Data{1,X}.date, '    offset = ', num2str(C(:,X(ii))*.33), ' seconds..']);
end

Bad_Dates = ROI_Data{1,X}.date;
end

Extraction_Metadata.Bad_Dates = Bad_Dates;
Extraction_Metadata.metadata = extraction_date;




