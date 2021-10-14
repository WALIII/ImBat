
function [GG_c] = ImBat_ConsistMap(x1,y1,centersF);
% map the concistancy of the bats behavior in bins in a 2D projection



bin_size = round(abs(centersF{1}(1)-centersF{1}(2)));
bin_size = bin_size*.35;
counter = 1;
for iii = 1: 50;
for ii = 1: 50

vects{counter} = find(x1>centersF{2}(ii)-bin_size & x1<centersF{2}(ii)+bin_size     & y1>centersF{1}(iii)-bin_size & y1<centersF{2}(iii)+bin_size);
counter = counter+1;

end
end


% theta matrix:     theta(i) = atan2(y1(i+60)-y1(i),x1(i+60)-x1(i));

for iii = 1:size(vects,2);
    try
        clear ttt slots2use2 vector_index theta
        ttt = abs(diff(x1(vects{iii})))+abs(diff(y1(vects{iii})));
        [slots2use, slots2use2] = findpeaks(ttt,'MinPeakProminence',10);
        
        vector_index = [1 slots2use2];
      % figure(); hold on;
        for i = 1:size(vector_index,2)-1
%              plot(x1(vects{iii}(vector_index(i)+1:vector_index(i+1)-1)),y1(vects{iii}(vector_index(i)+1:vector_index(i+1)-1)),'r')
            theta(i) = atan2(y1(vects{iii}(vector_index(i+1)-1))-y1(vects{iii}(vector_index(i)+1)),x1(vects{iii}(vector_index(i+1)-1))-x1(vects{iii}(vector_index(i)+1)));
        end
        
% = std(theta); % std around mean
       %thetaMat(iii) = 1 - sqrt((sum(sin(theta)))^2 + (sum(cos(theta)))^2)/size(theta,2);
 thetaMat(iii) = var(sin(theta)) + var(cos(theta));
    catch
        thetaMat(iii) = NaN;
    end
end
thetaMat2 = 1-mat2gray(thetaMat);
GG = reshape(thetaMat2,50,50);
GG_c = flipud(imrotate(GG,90));
GG_c = imgaussfilt(GG_c,2.5);
figure(); imagesc((GG_c))
colormap(turbo);


%% Scrap diagnostics:
% % play a movie
% 
% figure();
% for iix = 1:counter
% hold on;
% plot(x1,y1)
% 
% plot(x1(vects{iix}),y1(vects{iix}))
% 
% pause(0.01);
% clf
% end
% 
% % make a matrix:
% G(:,1) = x1(vects{1485});
% G(:,2) = y1(vects{1485});
% 
% 
% figure(); 
% hold on; 
% plot(x1,y1)
% plot(x1(vects{1485}),y1(vects{1485}))
% 
% figure(); 
% hold on;
% plot(abs(diff(x1(vects{1485}))));  plot(abs(diff(y1(vects{1485}))));
% ttt = abs(diff(x1(vects{1485})))+abs(diff(y1(vects{1485})));
% [slots2use, slots2use2] = findpeaks(ttt,'MinPeakProminence',10);
% figure(); plot(ttt); hold on; plot(slots2use2,slots2use,'*')
% 
% vector_index = [1 slots2use2];
% 
% 
% figure(); 
% hold on;
% for i = 1:size(vector_index,2)-1
% plot(x1(vects{1485}(vector_index(i)+1:vector_index(i+1)-1)),y1(vects{1485}(vector_index(i)+1:vector_index(i+1)-1)),'r')
% theta(i) = atan2(y1(vects{1485}(vector_index(i+1)-1))-y1(vects{1485}(vector_index(i)+1)),x1(vects{1485}(vector_index(i+1)-1))-x1(vects{1485}(vector_index(i)+1)));
% end

