function ImBat_YCplot(out_markov,flight2use);



% User inputs:
PltTxt = 1;

idx2use = find(out_markov.FlightIDVector == flight2use);

% col = colormap(lines(max(out_markov.FlightIDVector)));
% col = cat(1,[0.7 0.7 0.7],col);

cmap3 = colormap(lines(7));
cmap1 = colormap(cubehelix((max(out_markov.FlightIDVector)-7),1.5,3,4,1,[0.29,0.92]));

col = cat(1,cmap3,cmap1);
col = cat(1,[0.7 0.7 0.7],col);

% % Try CUBEHELIX
% col = colormap(cubehelix((max(out_markov.FlightIDVector)),1.5,3,4,1,[0.29,0.92]));
% col(1:3,:) = [];
% col = cat(1,col,[0.7 0.7 0.7],[0.7 0.7 0.7],[0.7 0.7 0.7]);
% col(1,:) = [0.7 0.7 0.7];


figure();
hold on;
LineW = 1/size(idx2use,2)*1000;

% resort based on:
resort_IDX = 2;
if resort_IDX ==1
    %1. length of flight:
    [a2 b2] = sort(out_markov.FL_len(idx2use),'ascend');
    idx2use = idx2use(b2);
    
elseif resort_IDX ==2
idx2use(idx2use ==1) = []; %if its the first flight, remove it
idx2use(idx2use ==size(out_markov.FlightIDVector,2)) = []; %if its the last flight, remove it

%2. flight ID pre
    [a2 b2] = sort(out_markov.FlightIDVector(idx2use-1),'ascend');
    idx2use = idx2use(b2);
    % now, re-sort subclusters
    idx4subsorting = out_markov.FlightIDVector(idx2use-1); % go through every pre-flight cluster
    %idx4subsorting = idx4subsorting(b2); % go through every pre-flight cluster
    idx4subsorting2 = out_markov.FlightIDVector(idx2use+1); %ssort based on following flight
    %idx4subsorting2 = idx4subsorting2(b2);
    fL2try = unique(idx4subsorting);
    fL2try2 = unique(idx4subsorting2);

    ind2cat = [];
    for i = 1:size(fL2try,2); % for every sub transition
        idxsort = find(idx4subsorting == fL2try(i));
        [ax bx] = sort(idx4subsorting2(idxsort));
        tempind = idx2use(idxsort);
        ind2cat = cat(2,ind2cat,tempind(bx));
        
        
        for ii = 1:size(fL2try2,2);
            Py2try{fL2try(i)}(fL2try2(ii))  = length(find(idx4subsorting2(idxsort) == fL2try2(ii)));
        end
        clear tempind
    end
    idx2use = ind2cat;
    
    
    
elseif resort_IDX ==3
    %3. flight ID post
    
    [a2 b2] = sort(out_markov.FlightIDVector(idx2use+1),'ascend');
    idx2use = idx2use(b2);
    idx4subsorting = out_markov.FlightIDVector(idx2use+1); % go through every pre-flight cluster
    fL2try = unique(idx4subsorting);
    
end
ind_temp_old = 0;

for i = 1:size(idx2use,2)
    clear Lpre Lactual Lpost;
    try
        % color
        col2use_pre = out_markov.FlightIDVector(idx2use(i)-1);
        col2use = out_markov.FlightIDVector(idx2use(i));
        col2use_post = out_markov.FlightIDVector(idx2use(i)+1);
        % % mean length
        % Lpre = round(-out_markov.ID_Len(col2use_pre)*100:0*100);
        % Lactual = round(0:out_markov.ID_Len(col2use)*100);
        % Lpost = round(out_markov.ID_Len(col2use)*100:(out_markov.ID_Len(col2use)+out_markov.ID_Len(col2use_post))*100);
        
        % individual length
        Lpre = round(-out_markov.FL_len(idx2use(i)-1)*100:0*100);
        Lactual = round(0:out_markov.FL_len(idx2use(i))*100);
        Lpost = round(out_markov.FL_len(idx2use(i))*100:(out_markov.FL_len(idx2use(i))+out_markov.FL_len(idx2use(i)+1))*100);
        
        plot(Lpre,ones(1,length(Lpre))*i,'color',col(col2use_pre,:),'LineWidth',LineW);
        plot(Lactual,ones(1,length(Lactual))*i,'color',col(col2use,:),'LineWidth',LineW)
        plot(Lpost,ones(1,length(Lpost))*i,'color',col(col2use_post,:),'LineWidth',LineW);
        if PltTxt ==1;
            ind_temp = idx4subsorting(i);
            if ind_temp_old == ind_temp;
            else
                text(-75000, i, cellstr(num2str(idx4subsorting(i))), 'FontSize', 10, 'Color', 'r');
                ind_temp_old = ind_temp;
            end
            
        end
    catch
        disp('missing data')
    end
end
xlabel('time (ms)');
ylabel('Flights');

A = out_markov.FlightIDVector(idx2use+1)';
B = unique(A);
out = [B,histc(A,B)];

% get pie
% for i = 1:length(out)
%     eventualPie(out(i,1)) = length(find(out_markov.FlightIDVector(idx2use+1)==fL2try2(i)));
% end

figure(); h = pie(out(:,2));
% change colors...
patchHand = findobj(h, 'Type', 'Patch');
for i = 1:size(patchHand,1);
    patchHand(i).FaceColor = col(out(i,1),:);
end

if resort_IDX == 2;
  figure(); 
  counter = 1;
  hold on;
  for ii = 1:size(Py2try,2)
      if sum(Py2try{ii})>10; % if a transition exists w over 10 flights
      try
      subplot(3,3,counter)
  h = pie(Py2try{ii});
% change colors...
patchHand = findobj(h, 'Type', 'Patch');
for i = 1:size(patchHand,1);
    patchHand(i).FaceColor = col(fL2try2(i),:);
end
title(['PreF: ',num2str(ii),' SumF: ',num2str(sum(Py2try{ii}))]);
      catch; end; counter = counter+1; end;
  end

end
