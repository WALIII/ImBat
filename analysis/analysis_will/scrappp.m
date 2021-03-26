


counter2 = 1;
for ii = 1:40
    GG = smooth(full((ROI_Data{1}.ROIs.results.S(ii,:))),10);
    [a2, b2] = (findpeaks(GG,'MinPeakProminence',.05,'MinPeakDistance',10));
    GG = (ROI_Data{1}.ROIs.results.C_raw(ii,:));

   figure();
    hold on;
    counter =1;
    if isempty(a2)<1
    for i = 1:size(a2,1);
        %try
        d2p1 = GG(b2(i)-1000:b2(i)+1000);
        d2pt(:,counter) = d2p1;
        plot((1:length(d2p1))/30-1000/30,d2p1,'k');
        counter = counter+1;
        %catch
           % disp('peak too close to edge');
        %end
    end
    if isempty(d2pt)<1 & size(d2pt,2)>2;
    plot((1:length(d2p1))/30-1000/30,mean(d2pt'),'LineWidth',6);
    
    % calculate FWHM
    d2peM = mean(d2pt');
    d2peM_max = max(d2peM);
    [aa] = find(d2peM>d2peM_max/2);
    Ca_FWHM(counter2) = aa(end)-aa(1)/30;
    counter2 = counter2+1;
    clear aa d2peM_max d2peM d2p b2 d2pt
  
    end
    else
    end
   
end

figure();
histogram(Ca_FWHM,20);

