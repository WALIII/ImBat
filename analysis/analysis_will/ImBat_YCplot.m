function ImBat_YCplot(out_markov,flight2use);


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
%2. flight ID pre
[a2 b2] = sort(out_markov.FlightIDVector(idx2use-1),'ascend');
idx2use = idx2use(b2);
% now, re-sort subclusters
idx4subsorting = out_markov.FlightIDVector(idx2use-1); % go through every pre-flight cluster
%idx4subsorting = idx4subsorting(b2); % go through every pre-flight cluster
idx4subsorting2 = out_markov.FlightIDVector(idx2use+1); %ssort based on following flight
%idx4subsorting2 = idx4subsorting2(b2);
fL2try = unique(idx4subsorting);
ind2cat = [];
for i = 1:size(fL2try,2);
    idxsort = find(idx4subsorting == fL2try(i));
    [ax bx] = sort(idx4subsorting2(idxsort));
    tempind = idx2use(idxsort);
    
    ind2cat = cat(2,ind2cat,tempind(bx));
    clear tempind
end
idx2use = ind2cat;



elseif resort_IDX ==3
%3. flight ID post
[a2 b2] = sort(out_markov.FlightIDVector(idx2use+1),'ascend');
idx2use = idx2use(b2);
end

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
    catch
        disp('missing data')
    end
end
xlabel('time (ms)');
ylabel('Flights');

