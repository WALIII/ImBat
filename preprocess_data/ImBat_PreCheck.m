function [Bad_Dates] = ImBat_PreCheck(ROI_Data);
% Check if all flight data is the same length

for day = 1:size(ROI_Data,2)

A = ROI_Data{1,day}.Alignment.out.Location_time(end);
B = ROI_Data{1,day}.Alignment.out.video_times(end);

C(:,day) = A-B;
% 
end
figure();
plot(C*.33);
title('Ending offsets')
xlabel('Days');
ylabel('Offset (seconds)');


X = find(abs(C)>120);
for ii = 1:size(X);
    disp(['WARNING: Failed alignment: ',ROI_Data{1,X}.date, '    offset = ', num2str(C(:,X(ii))*.33), ' seconds..']);

end

Bad_Dates = ROI_Data{1,X}.date;
% if C ==0
%     % go to the next 
% else
%     disp('error with day XXX');
% end
% 
% clear A B;
% end


