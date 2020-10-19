function ImBat_maxProfDff_wrapperUpdate
close all;
saveFlag = 1; %save data?
flyFlag = 0; %0 = fly, 1 = rest1, 2 = rest2

if saveFlag == 1
saveDir1 = '\\169.229.54.11\server_home\users\tobias\flight\data_processed\topQualityData\analysis_done\plots\';
% Check if folder exists
if exist([saveDir1 datestr(now,'yymmdd') filesep 'maxProjAll'])>0;
    disp('Youve been working today..');
else
    mkdir([saveDir1 datestr(now,'yymmdd') filesep 'maxProjAll'])
end
saveDir = [saveDir1 datestr(now,'yymmdd') filesep 'maxProjAll' '\'];
end

g = dir('Ge*');
z = dir('Zp*');
dirTop = vertcat(g,z); %find all folders in top quality directory

figAllDff = figure('units','normalized','outerposition',[0 0.1 0.4 0.8]);
sgtitle(['Max Projection Across Days ' dirTop(1).name(1:2)]);
figAllA = figure('units','normalized','outerposition',[0 0.1 0.4 0.8]);
sgtitle(['Corr Image Across Days ' dirTop(1).name(1:2)]);
%for each bat/session
for sesh_i = 1:length(dirTop)
    
    %get meta info for each bat/day
    cd([dirTop(sesh_i).name filesep 'extracted']);
    if flyFlag == 0
    dirFly = dir('*fly*extraction*');
    elseif flyFlag == 1
        dirFly = dir('*rest1*extraction');
    elseif flyFlag == 2
        dirFly = dir('*rest2*extraction');
    end
    batName = dirFly(1).name(1:3);
    dateSesh = dirFly(1).name(5:10);
    sessionType = dirFly(1).name(12:16);
    fileName = [batName '_' dateSesh '_' sessionType];
    %load cellData, flightpaths,alignment files
    cd(dirFly(1).name);
    dirAnal = dir('analysis*');
    
    dirDff = dir([dirAnal(end).folder filesep dirAnal(end).name filesep 'ROI' filesep '*maxProject.fig']); %find the max projection image
    figDff = openfig([dirDff.folder filesep dirDff.name],'reuse');    
    axDff = gca;
    figure(1)    
    s1 = subplot(ceil(length(dirTop)/3),3,sesh_i);
    fig1 = get(axDff,'children');
    copyobj(fig1,s1);
    axis image;
    hold on;
    colormap(gray);
    set(gca,'xticklabel',[],'yticklabel',[]);
    title([batName ' ' dateSesh ' ' sessionType]);
    drawnow
    close Figure 3
    
    dirA = dir([dirAnal(end).folder filesep dirAnal(end).name filesep 'ROI' filesep '*maxProjectROI.fig']);
    figA = openfig([dirA.folder filesep dirA.name],'reuse');
    axA = gca;
    figure(2)
    s2 = subplot(ceil(length(dirTop)/3),3,sesh_i);
    fig2 = get(axA,'children');
    copyobj(fig2,s2);
    axis image;
    colormap(gray);
    hold on;
    set(gca,'xticklabel',[],'yticklabel',[]);
    title([batName ' ' dateSesh ' ' sessionType]);
    drawnow
    close Figure 3
    
    cd(dirTop(sesh_i).folder);
end

if saveFlag == 1
    saveas(figAllDff,[saveDir filesep batName '_' dateSesh '_' sessionType '_maxProjAllDays_' datestr(now,'yymmdd_HHMM') '.tif']);
    savefig(figAllDff,[saveDir filesep batName '_' dateSesh '_' sessionType '_maxProjAllDays_' datestr(now,'yymmdd_HHMM') '.fig']);
    saveas(figAllA,[saveDir filesep batName '_' dateSesh '_' sessionType '_corrImAllDays_' datestr(now,'yymmdd_HHMM') '.tif']);
    savefig(figAllA,[saveDir filesep batName '_' dateSesh '_' sessionType '_corrImAllDays_' datestr(now,'yymmdd_HHMM') '.fig']);
end


