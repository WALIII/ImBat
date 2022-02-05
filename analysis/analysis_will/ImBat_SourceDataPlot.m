% Generating graphical figures from the Source Data file
% Liberti & Schmid et al, 2022

% WAL3
% 02.01.2022

% data organization:
% Source_Data.Figure.Panel

%% Figure 1: % N/A


%% Figure 2: % N/A


%% Figure 3:
    % fig3b: ROI stability vs Flight Variability (example)
    figure(); scatter(Source_Data.Fig3.b(:,1),Source_Data.Fig3.b(:,2)); title('');

    % fig3c: ROI stability vs Flight Variability (combined data)
    figure(); scatter(Source_Data.Fig3.c(:,1),Source_Data.Fig3.c(:,2)); title('');

    % fig3d: ROI stability vs Flight Variability (combined data)
    figure(); stairs(Source_Data.Fig3.d(:,1),Source_Data.Fig3.d(:,2:3)); title('');


%% Figure 4:
    % fig4f:
    figure(); stairs(Source_Data.Fig4.f(:,1),Source_Data.Fig4.f(:,2:4)); title('');
    % fig4g: 
    figure(); stairs(Source_Data.Fig4.g(:,1),Source_Data.Fig4.g(:,2:4)); title('');

%% Extended Data Figure 1: % N/A

%% Extended Data Figure 2: % N/A

%% Extended Data Figure 3:
    % ED Fig 3b
    figure(); plotSpread(Source_Data.EDFig3.b,'showMM',4);
    % ED Fig 3c
    figure(); scatter(Source_Data.EDFig3.c(:,1),Source_Data.EDFig3.c(:,2));
    % ED Fig 3d
    figure(); stairs([1 2 3 4 5],[Source_Data.EDFig3.d',0]); ylim([ 0 100]);    
%% Extended Data Figure 4:
    % ED Fig 4c
    figure(); stairs(Source_Data.EDFig4.Xaxis,Source_Data.EDFig4.c')
    % ED Fig 4d
    figure(); stairs(Source_Data.EDFig4.Xaxis,Source_Data.EDFig4.d')
    % ED Fig 4f
    figure(); stairs(Source_Data.EDFig4.Xaxis,Source_Data.EDFig4.f')
    % ED Fig 4g
    figure(); stairs(Source_Data.EDFig4.Xaxis,Source_Data.EDFig4.g')
    % ED Fig 4i
    figure(); plotSpread(Source_Data.EDFig4.i','showMM',4);%'Labels',{'day 1','day 2','day 3','day 4','day 5'})

%% Extended Data Figure 5:
    % ED Fig 5b
    figure(); stairs([Source_Data.EDFig5.Xaxis 2],[Source_Data.EDFig5.b,[0 0]']')
    % ED Fig 4c
    figure(); stairs([Source_Data.EDFig5.Xaxis 2],[Source_Data.EDFig5.c,[0 0]']')

%% Extended Data Figure 6:
    % ED Fig 6a
    figure(); scatter(Source_Data.EDFig6.a(:,1),Source_Data.EDFig6.a(:,2))
    % ED Fig 6c
    figure(); stairs([Source_Data.EDFig6.c',0])
    % ED Fig 6d
    figure(); plotSpread(Source_Data.EDFig6.d,'showMM',3);%'Labels',{'day 1','day 2','day 3','day 4','day 5'})
    % ED Fig 6e
    figure(); stairs(Source_Data.EDFig6.e(:,1), Source_Data.EDFig6.e(:,2))
    % ED Fig 6f
    figure(); plotSpread(Source_Data.EDFig6.f','showMM',3);%'Labels',{'day 1','day 2','day 3','day 4','day 5'})
    % ED Fig 6g
    figure(); stairs([0 Source_Data.EDFig6.g.bins],Source_Data.EDFig6.g.edges); set(gca,'YDir','reverse');
    % ED Fig 6h
    figure(); hold on; for i = 1:4; stairs(Source_Data.EDFig6.h.edges{i},[Source_Data.EDFig6.h.bins{i}, 0]); end

    
%% Extended Data Figure 7: % N/A


%% Extended Data Figure 8:
    % ED Fig 8a&b
   
    figure(); 
    subplot(1,2,1);
    plot(Source_Data.EDFig8.ab(:,1),Source_Data.EDFig8.ab(:,2),'o');
    subplot(1,2,2);
    plot(Source_Data.EDFig8.ab(:,1),Source_Data.EDFig8.ab(:,3),'o');


%% Extended Data Figure 9:
    % ED Fig 9a
    figure(); stairs(Source_Data.EDFig9.a(:,1),Source_Data.EDFig9.a(:,2:3))

%% Extended Data Figure 10:
    % ED Fig 10b
    figure(); hold on; for i = 1:5; stairs(Source_Data.EDFig10.d.edges{i},[Source_Data.EDFig10.d.counts{i}, 0]); end
    % ED Fig 10d
    figure(); hold on; for i = 1:5; stairs(Source_Data.EDFig10.d.edges{i},[Source_Data.EDFig10.d.counts{i}, 0]); end
    % ED Fig 10f
    figure(); plotSpread(Source_Data.EDFig10.f,'showMM',3)
    % ED Fig 10g
    figure(); plotSpread(Source_Data.EDFig10.g,'showMM',3)
    % ED Fig 10i
    figure(); plotSpread(Source_Data.EDFig10.i,'showMM',3)


