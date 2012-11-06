function runAnalysis3(handles)
%%%%% VERSION 1.2 12/15/10
%%%%% For Windows/Mac/Unix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get input parameters from GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FilesPath: The file to analyze
filesPath = get(handles.FilesPath,'String');
% outputPath: Where output information will be written
outputPath = get(handles.OutputPath,'String');
outputFile = fullfile(outputPath,get(handles.FileName,'String'));
% fishSet: The fish to be analyzed (see processing section)
% fishSet

ppoTime = str2num(get(handles.ppoTime,'String'));
% noiseThresh
noiseThresh = str2double(get(handles.noiseThresh,'String'));
%
plotGMVOT = get(handles.plotGMVOT,'Value');
plotCOVOT = get(handles.plotCOVOT,'Value');
plotNoise = get(handles.plotNoise,'Value');
plotHeatmap = get(handles.plotHeatmap,'Value');
plotIntensities = get(handles.plotIntensities,'Value');
plotHistograms = get(handles.plotHistograms,'Value');
plotPPO = get(handles.plotPPO,'Value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check Input and Output Information 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numMovies = 0;
if (isempty(filesPath) || isempty(outputPath) || strcmp(filesPath,'No files selected...') || strcmp(outputPath,'No directory selected...'))
    errordlg('You did not specify input file(s) or an output directory');
    set(handles.goBut,'Enable','on');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load mat file and further process parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(filesPath);
if(strcmp(get(handles.fishSet,'String'),'[:]'))
    fishSet{1} = [1:size(fishDistances,1)];
else 
    fishSet = str2cell(get(handles.fishSet,'String'));
end

if(strcmp(get(handles.ppoTime,'String'),'[:]'))
    ppoTime = [1:size(fishDistances,2)];
end

fprintf('\nAnalyzing data and generating figures, please wait...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for detection and display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smoothFactor = 10;
emptyWellThresh = 75; % min % of frames with NOEP before well is considered empty 
maxNoiseThresh = noiseThresh*100; % max % of frames that TMOEP or NOEP can be detected before well is thrown out


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert from pixels/frame to mm/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('frameRate')
    frameRate = 2;
end
[numWells,numFrames] = size(fishDistances);
%To add functionality for other wells, put the plate well number in the
%first set, and the corresponding well diameter in the same spot of the second 
wellDiameterConv = containers.Map({96,48,24,12,6},{6.78,10.5,15.62,22.1,34.8});  %% 7/16 are optionally used for 96/24 wells
%pix/frame * diameter(mm)/2*radius(pix) * frameRate frames/second
if ~exist('mmConv')
    mmConv = (wellDiameterConv(numWells)/(2*unscaledRadius));
end
if ~exist('convFact')
    convFact = mmConv*frameRate;
end
fishVelocities = fishDistances*convFact;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare output structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GroupOut = {'Wells','n', 'Mean Velocity (mm/s)', 'Mean Velocity SD', ...
    'Active Velocity (mm/s)', 'Active Velocity SD', '% Time Moving',...
    '% Time Moving SD', 'Active Duration (s)', 'Active Duration SD', ...
    'Rest Duration (s)', 'Rest Duration SD'};
IndividOut = {'Set','Well','Mean Velocity(mm/s)', 'Mean Velocity SD','Active Velocity(mm/s)','% Time Moving','Active Duration (s)', 'Rest Duration(s)'};
% To properly process 
firstTime = 1;
for setNum = 1:length(fishSet)
    %Get the current fish group
    fish = fishSet{setNum};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Assess well usability
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate a list of all possible usable wells
    list = 1:numWells;
    % Find wells with acceptable noObjectError percent (NOEP) and
    % tooManyObjectError percent (TMOEP)
    NOEP = noObjectError./numFrames*100;
    okNOEP = intersect(fish,list(NOEP < maxNoiseThresh));
    TMOEP = tooManyObjectError./numFrames*100;
    okTMOEP = intersect(fish,list(TMOEP < maxNoiseThresh));
    % Find empty wells (those with NOEP > empty threshold) 
    emptyWells = intersect(fish,list(NOEP > emptyWellThresh));
    % Find clean wells (those that have ok NOEP and TMOEP error rates)
    cleanWells = intersect(okNOEP,okTMOEP);
    % Find dirty wells (those that are not clean (high NOEP and TMOEP error rates))
    dirtyWells = setdiff(fish, cleanWells);
    % Find Usable Wells (those that are clean but not empty)
    usableWells = setdiff(cleanWells, emptyWells);
    fishSet{setNum} = usableWells;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Prepare data for analysis 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Some functions don't like scalar data (eg. ribbon plots)
    %  A trick to keep singular samples multidimentional
    %  is to simply double the sample to create a group.
    %  This has no effect on mean, std, or other measures, and will show
    %  an interfish variability of 0 (which is accurate).
    %  The only odd effect of this is that the ouput lists the sample twice and
    %  calls the group size 2. The alternative is that  everything below would 
    %  have to be rewritten for the degenerative matrix case (a vector). 
    if(length(usableWells) == 1)
        usableWells(end+1)=usableWells;
    end

    %If no wells are clean, skip this group
    if (isempty(usableWells))
        fprintf('\nWell(s) %s is(are) unusable (>5%% noise, >+/-2SD, empty, or n = 1\nTry increasing the noise threshold to 1\n',num2str(usableWells));
        return;%continue;
    end
    % Save old fishVelocities (for debugging purposes, and to restore when finished)
    oldFishVelocities = fishVelocities;
    % Store Distances only for usablewells
    fishVelocities = fishVelocities(usableWells,:);
    [numWells,numFrames] = size(fishVelocities);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculate Mean & Active Velocities, % Time 
    %%% Active and Group mean wrt Time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (plotGMVOT)
        % GMVOT = group mean velocity over time
        GMVOT = mean(fishVelocities);
        if (length(GMVOT) > 100)
            GMVOT = GMVOT(1:floor(.01*length(GMVOT)):length(GMVOT));
        else
            GMVOT = GMVOT(1:length(GMVOT));
        end
        GMVOT = gaussSmooth(GMVOT,floor(length(GMVOT)/smoothFactor)+1);
    end
    individMVs= mean(fishVelocities');
    individMVSTDs = std(fishVelocities');
    %%%% The sum of all velocities divided by the number of velocites > 0 is
    %%%% how active velocity is defined
    individAVs = sum(fishVelocities')./sum(fishVelocities'>0);
    %%%% The percent time movement is the number of velocites > 0 over the
    %%%% total number of velocities
    individTPs = (sum(fishVelocities'>0)./numFrames).*100;
    %%%% NOTE: You may want to change the 0's to something else for AV and TP
    %%%% depending on how you define movement and account for noise
    %Remove any NaN or Inf values from active velocity by setting them to 0
    individAVs((or(isinf(individAVs),isnan(individAVs))))=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculate burst and rest durations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    individADs(numWells) = 0;
    individRDs(numWells) = 0;
    for x = 1:numWells
        active = regionprops(im2bw(fishVelocities(x,:)),'Area');
        rest = regionprops(not(fishVelocities(x,:)),'Area');
        individADs(x) = mean([active.Area])/2;
        individRDs(x) = mean([rest.Area])/2;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Generate CoV data and smooth it for display
    %%% Note: This can be removed for improved performance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (plotCOVOT)
        if (numFrames>100)
            timePoint = floor(0.01*numFrames);
        else
            timePoint = 1; 
        end 
        CoV(length(1:timePoint:numFrames)) = 0;
        signal(length(1:timePoint:numFrames)) = 0;
        for x = timePoint:timePoint:numFrames
            i = floor(x/timePoint);
            if (numWells == 1)
                signal = mean(fishVelocities(1:x)');
                noise = std(fishVelocities(1:x)');
                CoV(i) = nanmean(noise./signal);
            else
                signal = mean(fishVelocities(:,1:x)');
                noise = std(fishVelocities(:,1:x)');
                CoV(i) = nanmean(noise./signal);
            end
                signals(i) = nanmean(signal);
        end
        CoV = gaussSmooth(CoV,floor(length(CoV)/smoothFactor));
    end
    %keyboard;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Display group information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nVideo Frame Rate: %.2f\tDistance Conversion Factor: %.2f',frameRate, mmConv);
    fprintf('\n\n_________________ Group %s _________________________________________',num2str(setNum));
    fprintf('\n================ Well Information =================================');
    fprintf('\nEmpty wells:      %s', num2str(emptyWells));
    fprintf('\nDiscarded wells: %s', num2str(setdiff(dirtyWells,emptyWells)));
    fprintf('\nAnalyzed wells:     %s', num2str(usableWells));
    fprintf('\n# of wells analyzed: %s', num2str(length(usableWells)));
    fprintf('\n================ Performance Analysis =============================');
    fprintf('\n      Type                 \tAverage     [Standard Deviation]   ');
    fprintf('\n"No object" errors:        \t%.2f %%', nanmean(NOEP(okNOEP)));
    fprintf('\n"Too many objects" errors: \t%.2f %%', nanmean(TMOEP(okTMOEP)));
    if (plotCOVOT)
        fprintf('\nMean Cv:                   \t%.2f     \t  [%.2f]', mean(CoV), std(CoV));
    end
    fprintf('\n================ Fish Activity ====================================');
    fprintf('\n      Type           \tAverage     [Standard Deviation]   ');
    fprintf('\nMean Velocity:       \t%.2f mm/s\t  [%.2f]', nanmean(individMVs), nanstd(individMVs));
    fprintf('\nActive Velocity:     \t%.2f mm/s\t  [%.2f]', nanmean(individAVs), nanstd(individAVs));
    fprintf('\nPercent Time Moving: \t%.1f %%   \t  [%.2f]\n', nanmean(individTPs), nanstd(individTPs));
    fprintf('\nActive Duration:     \t%.2f sec \t  [%.2f]',nanmean(individADs), nanstd(individADs));
    fprintf('\nRest Duration:       \t%.2f sec \t  [%.2f]',nanmean(individRDs), nanstd(individRDs));
    fprintf('\n___________________________________________________________________\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Store info for file output and graphing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GroupOut = cat(1,GroupOut,{usableWells,length(usableWells),mean(individMVs),std(individMVs),...
        mean(individAVs), std(individAVs), mean(individTPs), std(individTPs),...
        mean(individADs), std(individADs), mean(individRDs), std(individRDs)});
    if (firstTime)
        if (plotCOVOT)
            CoVplot = CoV;
        end
        if (plotGMVOT)
            GMVOTplot = GMVOT;
        end
        usedWells = usableWells;
        usedGroups = strcat(sprintf('Group %2d', setNum));
    else
        if (plotCOVOT)
            CoVplot = cat(1,CoVplot,CoV);
        end
        if (plotGMVOT)
            GMVOTplot = cat(1,GMVOTplot,GMVOT);
        end
        usedWells = cat(2,usedWells,usableWells);
        usedGroups = cat(1,usedGroups,strcat(sprintf('Group %2d', setNum)));
    end
    for i = 1:length(usableWells)
        IndividOut = cat(1,IndividOut,{setNum,usableWells(i),individMVs(i),individMVSTDs(i),individAVs(i),individTPs(i),individADs(i),individRDs(i)});
    end

    errorPlot = NOEP+TMOEP;

    fishVelocities = oldFishVelocities;
    [numWells,numFrames] = size(fishVelocities);
    firstTime = 0;

    % clean some values
    clear individMVs individAVs individTPs individADs individRDs usableWells;
end

%%%% Correct for the ribbon plot fix
if (usedWells(1)==usedWells(end))
    usedWells = usedWells(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write the group excel file mean and stds
GroupOut = cellfun(@num2str,GroupOut,'UniformOutput',false);
fid = fopen(strcat(outputFile(1:end-4),'_GROUP.xls'), 'wt');
[M,N] = size(GroupOut);
for i=1:M
    for j=1:N
        fprintf(fid, '%s\t',GroupOut{i,j});
    end
        fprintf(fid, '\n');
end
fclose(fid);
%%% Write the individual excel file mean and stds
temp = cellfun(@num2str,IndividOut,'UniformOutput',false);
fid = fopen(strcat(outputFile(1:end-4),'_INDIVID.xls'), 'wt');
[M,N] = size(temp);
for i=1:M
    for j=1:N
        fprintf(fid, '%s\t',temp{i,j});
    end
        fprintf(fid, '\n');
end
clear temp;
fclose(fid);
%save(strcat(outputFile(1:end-4),'_ANALYSIS.mat'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Graphs and thier associated functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% The main movement overview graph
overview_graph(fishVelocities,numWells,numFrames,usedWells)
%%% The error rate graph
if (plotNoise)
    errorRate_graph(errorPlot,1:numWells,maxNoiseThresh)
end
%%% The group mean velocity over time graph
if (plotGMVOT)
    ribbon_graph(GMVOTplot,GMVOT,usedGroups,fishSet, 'Group Mean Velocity');
end
%%% The coefficient of variation over time graph
if (plotCOVOT)
    ribbon_graph(CoVplot,CoV,usedGroups,fishSet, 'CoV');
end
%%% The intensity graphs for Vm Va and T%
if (exist('fishAreas') && plotIntensities)
    intensity_graph(unscaledRadius, fishAreas, cell2mat(IndividOut(2:end,3)), usedWells, 'Mean Velocity (mm/s)');
    intensity_graph(unscaledRadius, fishAreas, cell2mat(IndividOut(2:end,5)), usedWells, 'Active Velocity (mm/s)');
    intensity_graph(unscaledRadius, fishAreas, cell2mat(IndividOut(2:end,6)), usedWells, 'Percent Time Moving');
end

if (plotHistograms)
    histogram_graph(fishVelocities, usedWells, [0:.1:30])
end

if (exist('fishCoords') && plotPPO)
    ppo_graph(fishCoords,usedWells,ppoTime)
end

%%%% Hack to generate a heatmap from the coordinates for VPConv (slow)
if (~(exist('heatMap')==1) && exist('fishCoords'))
    heatMap = zeros([max(max(fishCoords(:,:,1)))+1,max(max(fishCoords(:,:,2)))+1]);
    for i = 1:size(fishCoords,2)
        for j = 1:size(fishCoords,1)
           heatMap(floor(fishCoords(j,i,1))+1,floor(fishCoords(j,i,2))+1) = heatMap(floor(fishCoords(j,i,1)+1),floor(fishCoords(j,i,2))+1)+1;
        end
    end
end

if ((exist('heatMap')==1) && plotHeatmap)
    heatmap_graph(heatMap)
end
save(strcat(outputFile(1:end-4),'_ANALYSIS.mat'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to plot overview graph (Fig 1F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function overview_graph(fishVelocities,numWells,numFrames,usedWells)
    figure;
    imagesc(fishVelocities);
    xlabel('Time');
    ylabel('Sample');
    set(gca,'XTickLabel', '');
    set(gca,'YTickLabel',[1,4:4:(numWells*4)]);
    set(gca,'YTick',1:4:numWells);
    hold on;
    for x = 1:numWells-1
        line([1 numFrames], [x+.5 x+.5], 'color', 'w');
    end
    colormap(flipud(gray(16)));
    plot(1,usedWells,'gs','MarkerFaceColor','g','MarkerSize',5)
    hold off;
    title('');
    colorbar('YTickLabel',[0:max(max(fishVelocities))]);
    %%% User now required to manually save this image if it is wanted
    %%% saveas(gca, strcat(PathName,FileName(1:end-4),'.jpg'),'jpg');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Graph error rates (Sup 2B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function errorRate_graph(errorPlot,usedWells,maxNoiseThresh)
    %%% Graph customization
    meanLineColor = [.15 .23 .37];
    meanLineStyle = '--';
    stdLineColor = [.15 .23 .37];
    stdLineStyle = ':';
    barColor = [.89 .94 .9];
    %'–' Solid line (default)
    %'--'Dashed line
    %':'Dotted line
    %'–.' Dash-dot line
    %'none' No line
    
    figure;
    bar(errorPlot, 'FaceColor', barColor);
    %Alter the x-axis title and tick properties
    xlabel('Sample Number');
    xlim([0.5 length(usedWells)+0.5]);
    set(gca,'XTick',1:length(usedWells));
    set(gca,'XTickLabel', num2str(usedWells'));
    set(gca,'FontSize',9);
    ylabel('Combined error %');
    ylim([0 (2*maxNoiseThresh)+1]);
    title('Combined Error %s for each well');
    %Add the threshold line, mean, and +1/-1STD lines
    hold on;
        line([0 length(usedWells)+1], [2*maxNoiseThresh 2*maxNoiseThresh], 'color', 'r', 'LineStyle', '-');  
        line([0 length(usedWells)+1], [maxNoiseThresh maxNoiseThresh], 'color', 'y', 'LineStyle', '-'); 
        line([0 length(usedWells)+1], [mean(errorPlot) mean(errorPlot)], 'color', meanLineColor, 'LineStyle', meanLineStyle);
        if (mean(errorPlot)-std(errorPlot)>0)
            line([0 length(usedWells)+1], [mean(errorPlot)-std(errorPlot) mean(errorPlot)-std(errorPlot)], 'color', stdLineColor, 'LineStyle', stdLineStyle);
        end
        line([0 length(usedWells)+1], [mean(errorPlot)+std(errorPlot) mean(errorPlot)+std(errorPlot)], 'color', stdLineColor, 'LineStyle', stdLineStyle);
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to graph ribbon plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ribbon_graph(plotData, data, usedGroups, fishSet, ttitle)
    figure;
    ribbon(plotData');
    legend(gca,usedGroups);
    xlim([0.5 size(usedGroups,1)+0.5]);
    set(gca,'XTick',1:length(fishSet));
    ylabel('Time (min)');
    ylim([1 length(data)]);
    set(gca,'YTick',floor(length(data)/10):floor(length(data)/10):length(data));
    %%% Tick labels should be adjusted for frame rate, etc... disabled for now
    %set(gca,'YTickLabel', floor(numFrames/2/60/10):floor(numFrames/2/60/10):floor(numFrames/2/60));
    %%zlabel('Mean Velocity (mm/s)');
    title(ttitle);
    view(68,16); %% Change the view so the perspective is slightly off angle
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to create intensity plots (Fig 4F) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intensity_graph(unscaledRadius, fishAreas, values, usedWells, graphTitle)
    figure;
    hold on;
    %Needed for conversion from matrix to image coords. 
    maxx = max(mean(fishAreas(:,3:4)'));
    %Compute the well centers using the fish area bounding box
    fishCenters = [maxx-mean(fishAreas(:,1:2)'); mean(fishAreas(:,3:4)');];
    fishCenters = fishCenters';
    plot(fishCenters(:,2),fishCenters(:,1),'k.');

    t = linspace(0,2*pi,1000);
    r = unscaledRadius;
    j = 1;
    maxIntensity = max(values);
    minIntensity = min(values);
    for i = usedWells
        h = fishCenters(i,2);
        k = fishCenters(i,1);
        x = r*cos(t)+h;
        y = r*sin(t)+k;
        %intensity = (values(j)-minIntensity)/(maxIntensity-minIntensity);
        %intensity = values(j)/maxIntensity;
        intensity = values(j);
        set(gca,'Clim',[minIntensity,maxIntensity]);
        j = j+1;
        %fill(x,y,[intensity intensity intensity]);
        fill(x,y,intensity);
    end
    hold off;
    axis off equal;
    title(graphTitle);
    colormap(flipud(gray(128)));
    %colorbar('YTickLabel',floor([minIntensity:floor(maxIntensity-minIntensity)/8:maxIntensity].*100)/100);
    colorbar('YTickLabel',[minIntensity:((maxIntensity-minIntensity)/8):maxIntensity]);
end

function histogram_graph(fishVelocities, wells, range)
    instV = reshape(fishVelocities(wells,:), 1, size(fishVelocities,2)*length(wells)); 
    histDist = histc(instV,range);
    wHistDist = histDist.*range;
    figure;
    bar(range,wHistDist);
    xlim([min(range) max(range)]);
    %%% Allow matlab to autoscale
    %ylim([0 3000]);
end

function ppo_graph(fishCoords, wells, range)
    figure;
    if (range(1)>0 && range(end) < size(fishCoords,2))
        plotPathOverlay(fishCoords(wells,range,:));
    elseif (range(1) && range(end) >= size(fishCoords,2))
        range = range(1):size(fishCoords,2);
        plotPathOverlay(fishCoords(wells,range,:));
    else
        plotPathOverlay(fishCoords(wells,:,:));
    end
end

function heatmap_graph(heatMap)
    figure;
    % Viewpoint uses coorinates (0,0)[(1,1) in matlab] and (1,0)[(2,1) in
    % matlab] for special use, 0 these to prevent oversaturation of the
    % heatmap
    heatMap(1,1) = 0;
    heatMap(2,1) = 0;
    %colormap(hot(128));
    colormap(flipud(gray(max(max(heatMap)))));
    image(heatMap);
    daspect([1 1 1]);
    axis off;
end
