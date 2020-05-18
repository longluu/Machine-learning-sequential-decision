%%%%%%%% Data analysis of conditional decision making experiment %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
% 

subjectID = 'average';
experimentNumber = 1:7;
experimentType = 'MainExperiment';
experiment = 'Original';
session = 1;
dataAll = [];
fontSize = 15;
includeIncongruent = 0; % 0: no incongruent
                        % 1: include incongruent
                        % 2: only incongruent
lineWidth = 2;

for ii = 1 : length(experimentNumber)
    if strcmp(experiment, 'PilotData')
        dataFullPath = fullfile('Data', subjectID, experimentType, experiment, ['Session' num2str(session)], ...
                                    ['ConditionDecisionMaking-' num2str(experimentNumber(ii))]);
    elseif strcmp(experimentType, 'MotorNoise')
        dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                    ['MotorNoise-' num2str(experimentNumber(ii))]);
    elseif strcmp(experimentType, 'PerceptNoise')
        dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                    ['PerceptNoise-' num2str(experimentNumber(ii))]);                                
    elseif strcmp(experimentType, 'TrainArray')
        dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                    ['TrainArray-' num2str(experimentNumber(ii))]);  
    else
        dataFullPath = fullfile('Data', subjectID, experimentType, [experiment num2str(session)], ...
                                    [experiment '-' num2str(experimentNumber(ii))]);          
    end
    load(dataFullPath);
    if ii == 1
        dataAll = dataTotal;
        if strcmp(experimentType, 'PerceptNoise') && isfield(params, 'staircaseTrack')
            staircaseTrack = params.staircaseTrack;
        end
    else
        dataAll = [dataAll; dataTotal];
        if strcmp(experimentType, 'PerceptNoise') && isfield(params, 'staircaseTrack')
            staircaseTrack = [staircaseTrack; params.staircaseTrack];
        end
    end
end
biasMotor = zeros(1,15);  
if strcmp(experimentType, 'MotorNoise')
    a = dataAll(:,7);
    dataAll(isnan(a),:) = [];    
end

%% Extract data
% Create the data array.  
%  Column 1: true angle difference (population)
%  Column 2: bar presentation time
%  Column 3: SOA
%  Column 4: bar noise level
%  Column 5: bar reference angle
%  Column 6 to 8: data (categorical decision and estimation)
%  Column 9: true angle difference (sample) 
%  Column 10: decision time
%  Column 11: estimation time
stdNoiseLevel = unique(dataAll(:,4));
angleDiff = (unique(dataAll(:,1)))';
dataZeroDiff = dataAll(dataAll(:,1) == 0,:);
angleEstimate = dataAll(:, 7);
decisionTimeAll = dataAll(:, 10);
estimateTimeAll = dataAll(:, 11);
madEstimateTime = 1.5*median(abs(estimateTimeAll - median(estimateTimeAll)));
estimateTimeAll(abs(estimateTimeAll-median(estimateTimeAll)) > 10*madEstimateTime) = NaN; 
if strcmp(experimentType, 'MotorNoise')
    stdNoiseLevel = 0;
end

%% Average across conditions
angleTrue = dataAll(:,5) - dataAll(:,1);
angleTrue(angleTrue>180 | angleTrue<0) = angleTrue(angleTrue>180 | angleTrue<0)...
                                            - 180 * sign(angleTrue(angleTrue>180 | angleTrue<0));
indexFixing = find(abs(angleTrue-angleEstimate) > 70);
angleEstimate(indexFixing) = angleEstimate(indexFixing) + sign(angleTrue(indexFixing)-angleEstimate(indexFixing))*180;    
angleDiffEst = dataAll(:, 5) - angleEstimate;
indexAdjust1 = abs(angleDiffEst)>90;
indexAdjust2 = sign(angleDiffEst) ~= sign(dataAll(:, 1));
indexAdjust = indexAdjust1 & indexAdjust2;
angleDiffEst(indexAdjust) = angleDiffEst(indexAdjust) - 180 * sign(angleDiffEst(indexAdjust));
biasAll = angleTrue - angleEstimate;
binaryFeedbackAll = dataAll(:,8);  
binaryDecisionAll = dataAll(:,6);  
if includeIncongruent == 0
    if sum(isnan(binaryFeedbackAll)) == length(binaryFeedbackAll)
        indicatorConsistent = ones(length(binaryFeedbackAll));
    else
        indicatorConsistent = sign(angleDiffEst) == binaryFeedbackAll;
    end
    angleDiffEst(~indicatorConsistent) = NaN;
elseif includeIncongruent == 2
    indicatorConsistent = sign(angleDiffEst) == binaryFeedbackAll;
    angleDiffEst(indicatorConsistent) = NaN;    
end

angleDiffEstimate = cell(length(stdNoiseLevel), length(angleDiff));
decisionTime = cell(length(stdNoiseLevel), length(angleDiff));
estimateTime = cell(length(stdNoiseLevel), length(angleDiff));
biasEstimate = cell(length(stdNoiseLevel), length(angleDiff));
binaryDecisionSelf = cell(length(stdNoiseLevel), length(angleDiff));
binaryFeedback = binaryDecisionSelf;
for ii = 1 : length(stdNoiseLevel)
    for jj = 1 : length(angleDiff)
        tempBias = biasAll(dataAll(:,4)==stdNoiseLevel(ii) ...
                           & dataAll(:,1)==angleDiff(jj), :);
        biasEstimate{ii,jj} = tempBias(:);
        tempangleDiffEstimate = angleDiffEst(dataAll(:,4)==stdNoiseLevel(ii) ...
                                            & dataAll(:,1)==angleDiff(jj), :);
        binaryDecisionSelf{ii,jj} = binaryDecisionAll(dataAll(:,4)==stdNoiseLevel(ii) ...
                                        & dataAll(:,1)==angleDiff(jj), :);  
        binaryFeedback{ii,jj} = binaryFeedbackAll(dataAll(:,4)==stdNoiseLevel(ii) ...
                                        & dataAll(:,1)==angleDiff(jj), :);  
        angleDiffEstimate{ii,jj} = tempangleDiffEstimate;
        decisionTime{ii,jj} = decisionTimeAll(dataAll(:,4)==stdNoiseLevel(ii) ...
                                            & dataAll(:,1)==angleDiff(jj), :);
        estimateTime{ii,jj} = estimateTimeAll(dataAll(:,4)==stdNoiseLevel(ii) ...
                                            & dataAll(:,1)==angleDiff(jj), :);
                                        
    end
end
binaryDecision = binaryDecisionSelf;  

%% Plot heat map of reaction time
% Correct trial
figure
hold on

% Decision RT
maxEstimate = 40;
minEstimate = -40;
yAxis = minEstimate:2:maxEstimate;
xAxis = angleDiff;
lengthXaxis = length(angleDiff(1):angleDiff(end)) + diff(angleDiff(end-1:end)) -1;
histEstimate = cell(1, length(stdNoiseLevel));
for ii = 1 : length(stdNoiseLevel)
    histEstimate_temp = zeros(length(yAxis), lengthXaxis);
    counter = 1;
    for jj = 1 :length(angleDiff)
        tempAngleDiffEst = angleDiffEstimate{ii,jj};
        indExclude = binaryDecision{ii,jj} ~= binaryFeedback{ii,jj};
        tempAngleDiffEst(indExclude) = [];
        tempDecisionTime = decisionTime{ii,jj};
        tempDecisionTime(indExclude) = [];
        
        tempHeatMap = computeHeatMap(tempAngleDiffEst, tempDecisionTime, yAxis); 
        if angleDiff(jj) < 0
            histEstimate_temp(:, counter:counter+1) = repmat(tempHeatMap', [1 2]); 
            counter = counter + 2;
        else
            histEstimate_temp(:, counter:counter+4) = repmat(tempHeatMap', [1 5]); 
            counter = counter + 5;            
        end
    end
    histEstimate{ii} = histEstimate_temp;
    subplot(2, length(stdNoiseLevel), ii)
    hold on
    imagesc(histEstimate_temp, 'AlphaData', ~isnan(histEstimate_temp))
    colormap hot
    set(gca,'color', 0.7*[1 1 1]); 
    colorbar
    axis xy
    plot([2 lengthXaxis-2], [find(yAxis == xAxis(1)) find(yAxis == xAxis(end))], ':', 'Color', 'white', 'LineWidth', 1.5)
    plot([1 lengthXaxis], [round(length(yAxis)/2) round(length(yAxis)/2)],  ':', 'Color', 'white', 'LineWidth', 1.5)
    plot([15 15], [1 length(yAxis)],  ':', 'Color', 'white', 'LineWidth', 1.5)    
    xDisplay = -10:10:30;
    xTick = round((xDisplay + 12) * lengthXaxis / (angleDiff(end)-angleDiff(1)));
    xTick(xDisplay<=0) = xTick(xDisplay<=0) + 2;
    yDisplay = round(linspace(yAxis(1) ,yAxis(end),5));
    set(gca, 'ylim', [1 length(yAxis)], 'xlim', [1 lengthXaxis], ...
        'XTick', xTick , 'XTickLabel', num2cell(xDisplay),...
        'YTick', find(ismember(yAxis, yDisplay)), 'YTickLabel', num2cell(yDisplay), ...
        'FontSize', fontSize)
end

% Estimation RT
histEstimate = cell(1, length(stdNoiseLevel));
for ii = 1 : length(stdNoiseLevel)
    histEstimate_temp = zeros(length(yAxis), lengthXaxis);
    counter = 1;
    for jj = 1 :length(angleDiff)
        tempAngleDiffEst = angleDiffEstimate{ii,jj};
        indExclude = binaryDecision{ii,jj} ~= binaryFeedback{ii,jj};
        tempAngleDiffEst(indExclude) = [];
        tempEstimateTime = estimateTime{ii,jj};
        tempEstimateTime(indExclude) = [];
        
        tempHeatMap = computeHeatMap(tempAngleDiffEst, tempEstimateTime, yAxis); 
        if angleDiff(jj) < 0
            histEstimate_temp(:, counter:counter+1) = repmat(tempHeatMap', [1 2]); 
            counter = counter + 2;
        else
            histEstimate_temp(:, counter:counter+4) = repmat(tempHeatMap', [1 5]); 
            counter = counter + 5;            
        end
    end
    histEstimate{ii} = histEstimate_temp;
    subplot(2, length(stdNoiseLevel), ii+length(stdNoiseLevel))
    hold on
    imagesc(histEstimate_temp, 'AlphaData', ~isnan(histEstimate_temp))
    colormap hot
    set(gca,'color', 0.7*[1 1 1]); 
    colorbar
    axis xy
    plot([2 lengthXaxis-2], [find(yAxis == xAxis(1)) find(yAxis == xAxis(end))], ':', 'Color', 'white', 'LineWidth', 1.5)
    plot([1 lengthXaxis], [round(length(yAxis)/2) round(length(yAxis)/2)],  ':', 'Color', 'white', 'LineWidth', 1.5)
    plot([15 15], [1 length(yAxis)],  ':', 'Color', 'white', 'LineWidth', 1.5)    
    xDisplay = -10:10:30;
    xTick = round((xDisplay + 12) * lengthXaxis / (angleDiff(end)-angleDiff(1)));
    xTick(xDisplay<=0) = xTick(xDisplay<=0) + 2;
    yDisplay = round(linspace(yAxis(1) ,yAxis(end),5));
    set(gca, 'ylim', [1 length(yAxis)], 'xlim', [1 lengthXaxis], ...
        'XTick', xTick , 'XTickLabel', num2cell(xDisplay),...
        'YTick', find(ismember(yAxis, yDisplay)), 'YTickLabel', num2cell(yDisplay), ...
        'FontSize', fontSize)
end

% Incorrect trial
figure
hold on

% Decision RT
maxEstimate = 40;
minEstimate = -40;
yAxis = minEstimate:2:maxEstimate;
xAxis = angleDiff;
lengthXaxis = length(angleDiff(1):angleDiff(end)) + diff(angleDiff(end-1:end)) -1;
histEstimate = cell(1, length(stdNoiseLevel));
for ii = 1 : length(stdNoiseLevel)
    histEstimate_temp = zeros(length(yAxis), lengthXaxis);
    counter = 1;
    for jj = 1 :length(angleDiff)
        tempAngleDiffEst = angleDiffEstimate{ii,jj};
        indExclude = binaryDecision{ii,jj} == binaryFeedback{ii,jj};
        tempAngleDiffEst(indExclude) = [];
        tempDecisionTime = decisionTime{ii,jj};
        tempDecisionTime(indExclude) = [];
        
        tempHeatMap = computeHeatMap(tempAngleDiffEst, tempDecisionTime, yAxis); 
        if angleDiff(jj) < 0
            histEstimate_temp(:, counter:counter+1) = repmat(tempHeatMap', [1 2]); 
            counter = counter + 2;
        else
            histEstimate_temp(:, counter:counter+4) = repmat(tempHeatMap', [1 5]); 
            counter = counter + 5;            
        end
    end
    histEstimate{ii} = histEstimate_temp;
    subplot(2, length(stdNoiseLevel), ii)
    hold on
    imagesc(histEstimate_temp, 'AlphaData', ~isnan(histEstimate_temp))
    colormap hot
    set(gca,'color', 0.7*[1 1 1]); 
    colorbar
    axis xy
    plot([2 lengthXaxis-2], [find(yAxis == xAxis(1)) find(yAxis == xAxis(end))], ':', 'Color', 'white', 'LineWidth', 1.5)
    plot([1 lengthXaxis], [round(length(yAxis)/2) round(length(yAxis)/2)],  ':', 'Color', 'white', 'LineWidth', 1.5)
    plot([15 15], [1 length(yAxis)],  ':', 'Color', 'white', 'LineWidth', 1.5)    
    xDisplay = -10:10:30;
    xTick = round((xDisplay + 12) * lengthXaxis / (angleDiff(end)-angleDiff(1)));
    xTick(xDisplay<=0) = xTick(xDisplay<=0) + 2;
    yDisplay = round(linspace(yAxis(1) ,yAxis(end),5));
    set(gca, 'ylim', [1 length(yAxis)], 'xlim', [1 lengthXaxis], ...
        'XTick', xTick , 'XTickLabel', num2cell(xDisplay),...
        'YTick', find(ismember(yAxis, yDisplay)), 'YTickLabel', num2cell(yDisplay), ...
        'FontSize', fontSize)
end

% Estimation RT
histEstimate = cell(1, length(stdNoiseLevel));
for ii = 1 : length(stdNoiseLevel)
    histEstimate_temp = zeros(length(yAxis), lengthXaxis);
    counter = 1;
    for jj = 1 :length(angleDiff)
        tempAngleDiffEst = angleDiffEstimate{ii,jj};
        indExclude = binaryDecision{ii,jj} == binaryFeedback{ii,jj};
        tempAngleDiffEst(indExclude) = [];
        tempEstimateTime = estimateTime{ii,jj};
        tempEstimateTime(indExclude) = [];
        
        tempHeatMap = computeHeatMap(tempAngleDiffEst, tempEstimateTime, yAxis); 
        if angleDiff(jj) < 0
            histEstimate_temp(:, counter:counter+1) = repmat(tempHeatMap', [1 2]); 
            counter = counter + 2;
        else
            histEstimate_temp(:, counter:counter+4) = repmat(tempHeatMap', [1 5]); 
            counter = counter + 5;            
        end
    end
    histEstimate{ii} = histEstimate_temp;
    subplot(2, length(stdNoiseLevel), ii+length(stdNoiseLevel))
    hold on
    imagesc(histEstimate_temp, 'AlphaData', ~isnan(histEstimate_temp))
    colormap hot
    set(gca,'color', 0.7*[1 1 1]); 
    colorbar
    axis xy
    plot([2 lengthXaxis-2], [find(yAxis == xAxis(1)) find(yAxis == xAxis(end))], ':', 'Color', 'white', 'LineWidth', 1.5)
    plot([1 lengthXaxis], [round(length(yAxis)/2) round(length(yAxis)/2)],  ':', 'Color', 'white', 'LineWidth', 1.5)
    plot([15 15], [1 length(yAxis)],  ':', 'Color', 'white', 'LineWidth', 1.5)    
    xDisplay = -10:10:30;
    xTick = round((xDisplay + 12) * lengthXaxis / (angleDiff(end)-angleDiff(1)));
    xTick(xDisplay<=0) = xTick(xDisplay<=0) + 2;
    yDisplay = round(linspace(yAxis(1) ,yAxis(end),5));
    set(gca, 'ylim', [1 length(yAxis)], 'xlim', [1 lengthXaxis], ...
        'XTick', xTick , 'XTickLabel', num2cell(xDisplay),...
        'YTick', find(ismember(yAxis, yDisplay)), 'YTickLabel', num2cell(yDisplay), ...
        'FontSize', fontSize)
end

%% Plot reaction time as a function of stimulus (collapsed across estimate)
color_default = get(gca,'colororder');
close

% Correct trial
rtDecisionCorrect = NaN(length(stdNoiseLevel), length(angleDiff));
rtEstimateCorrect = NaN(length(stdNoiseLevel), length(angleDiff));
rtEstCorrect_individualTrial = cell(1, 2);
n_trials_correct = NaN(length(stdNoiseLevel), length(angleDiff));
std_correct = NaN(length(stdNoiseLevel), length(angleDiff));
rtEstimate_sem_correct = NaN(length(stdNoiseLevel), length(angleDiff));
hLegend = NaN(1, length(stdNoiseLevel));
figure
for ii = 1 : length(stdNoiseLevel)
    tempRtDecision_sem = NaN(1, length(angleDiff));
    for jj = 1 : length(angleDiff)
        indExclude = binaryDecision{ii,jj} ~= binaryFeedback{ii,jj};
        tempRTDecision = decisionTime{ii, jj};
        tempRTDecision(indExclude) = [];
        tempRTEstimate = estimateTime{ii, jj};
        tempRTEstimate(indExclude) = [];
        tempRTDecision(isnan(tempRTDecision)) = [];
        tempRTEstimate(isnan(tempRTEstimate)) = [];
        rtDecisionCorrect(ii, jj) = mean(tempRTDecision);
        rtEstimateCorrect(ii, jj) = mean(tempRTEstimate);
        rtEstCorrect_individualTrial{ii} = [rtEstCorrect_individualTrial{ii} tempRTEstimate'];
        tempRtDecision_sem(jj) = std(tempRTDecision) / sqrt(length(tempRTDecision));
        rtEstimate_sem_correct(ii, jj) = std(tempRTEstimate) / sqrt(length(tempRTEstimate));
        n_trials_correct(ii, jj) = length(tempRTEstimate);
        std_correct(ii, jj) = std(tempRTEstimate); 
    end

    % Plot
    subplot(2, 2, 1)
    hold on
    hShade = shadedErrorBar(angleDiff, rtDecisionCorrect(ii, :), tempRtDecision_sem,... 
                         {'Color', color_default(ii, :), 'LineWidth', lineWidth});     
    hLegend(ii) = hShade.mainLine;
    xlabel('Stimulus orientation')
    ylabel('Reaction time')
    title('Decision')
    xlim([min(angleDiff) max(angleDiff)])

    subplot(2, 2, 2)
    hold on
    shadedErrorBar(angleDiff, rtEstimateCorrect(ii, :), rtEstimate_sem_correct(ii, :),... 
                         {'Color', color_default(ii, :), 'LineWidth', lineWidth});        
    xlabel('Stimulus orientation')
    ylabel('Reaction time')
    title('Estimation')
    xlim([min(angleDiff) max(angleDiff)])   
    ylim([1.5 5.5])
end
legend(hLegend, {'Low noise', 'High noise'})

% Incorrect trial
rtDecisionIncorrect = NaN(length(stdNoiseLevel), length(angleDiff));
rtEstimateIncorrect = NaN(length(stdNoiseLevel), length(angleDiff));
n_trials_incorrect = NaN(length(stdNoiseLevel), length(angleDiff));
std_incorrect = NaN(length(stdNoiseLevel), length(angleDiff));
rtEstimate_sem_incorrect = NaN(length(stdNoiseLevel), length(angleDiff));
rtEstIncorrect_individualTrial = cell(1, 2);
for ii = 1 : length(stdNoiseLevel)
    tempRtDecision_sem = NaN(1, length(angleDiff));
    tempRtEstimate_sem = NaN(1, length(angleDiff));
    for jj = 1 : length(angleDiff)
        indExclude = binaryDecision{ii,jj} == binaryFeedback{ii,jj};
        tempRTDecision = decisionTime{ii, jj};
        tempRTDecision(indExclude) = [];
        tempRTEstimate = estimateTime{ii, jj};
        tempRTEstimate(indExclude) = [];
        tempRTDecision(isnan(tempRTDecision)) = [];
        tempRTEstimate(isnan(tempRTEstimate)) = [];
        rtDecisionIncorrect(ii, jj) = mean(tempRTDecision);
        rtEstimateIncorrect(ii, jj) = mean(tempRTEstimate);
        rtEstIncorrect_individualTrial{ii} = [rtEstIncorrect_individualTrial{ii} tempRTEstimate'];
        tempRtDecision_sem(jj) = std(tempRTDecision) / sqrt(length(tempRTDecision));
        rtEstimate_sem_incorrect(ii, jj) = std(tempRTEstimate) / sqrt(length(tempRTEstimate));
        n_trials_incorrect(ii, jj) = length(tempRTEstimate);
        std_incorrect(ii, jj) = std(tempRTEstimate); 
    end

    % Plot
    subplot(2, 2, 3)
    hold on
    hShade = shadedErrorBar(angleDiff, rtDecisionIncorrect(ii, :), tempRtDecision_sem,... 
                         {'Color', color_default(ii, :), 'LineWidth', lineWidth});     
    hLegend(ii) = hShade.mainLine;
    xlabel('Stimulus orientation')
    ylabel('Reaction time')
    title('Decision')
    xlim([min(angleDiff) max(angleDiff)])

    subplot(2, 2, 4)
    hold on
    shadedErrorBar(angleDiff, rtEstimateIncorrect(ii, :), rtEstimate_sem_incorrect(ii, :),... 
                         {'Color', color_default(ii, :), 'LineWidth', lineWidth});        
    xlabel('Stimulus orientation')
    ylabel('Reaction time')
    title('Estimation')
    xlim([min(angleDiff) max(angleDiff)]) 
    ylim([1.5 5.5])
end

figure
subplot(3, 2, 1)
plot(angleDiff, n_trials_correct, 'o-')
xlabel('Stimulus orientation')
ylabel('Trials')
title('Correct trials')
xlim([min(angleDiff) max(angleDiff)])    

subplot(3, 2, 2)
plot(angleDiff, n_trials_incorrect, 'o-')
xlabel('Stimulus orientation')
ylabel('Trials')
title('Incorrect trials')
xlim([min(angleDiff) max(angleDiff)])    

subplot(3, 2, 3)
plot(angleDiff, std_correct, 'o-')
xlabel('Stimulus orientation')
ylabel('Std')
title('Correct trials')
xlim([min(angleDiff) max(angleDiff)])    

subplot(3, 2, 4)
plot(angleDiff, std_incorrect, 'o-')
xlabel('Stimulus orientation')
ylabel('Std')
title('Incorrect trials')
xlim([min(angleDiff) max(angleDiff)])   

subplot(3, 2, 5)
plot(angleDiff, rtEstimate_sem_correct, 'o-')
xlabel('Stimulus orientation')
ylabel('SEM')
title('Correct trials')
xlim([min(angleDiff) max(angleDiff)])    

subplot(3, 2, 6)
plot(angleDiff, rtEstimate_sem_incorrect, 'o-')
xlabel('Stimulus orientation')
ylabel('SEM')
title('Incorrect trials')
xlim([min(angleDiff) max(angleDiff)])   

% Plot the average for each panel in the figure above
figure
hold on
rt_correct = NaN(2, 3);
rt_incorrect = NaN(2, 3);
n_boot_sample = 10000;
for ii = 1 : length(stdNoiseLevel)
    rt_correct(ii, 1) = mean(rtEstCorrect_individualTrial{ii});
    bootstat = bootstrp(n_boot_sample, @mean, rtEstCorrect_individualTrial{ii});
    sem = std(bootstat);
    rt_correct(ii, 2) = rt_correct(ii, 1) - sem;
    rt_correct(ii, 3) = rt_correct(ii, 1) + sem;
%     rt_correct(ii, 2:3) = bootci(n_boot_sample, @mean, rtEstCorrect_individualTrial{ii});
    
    rt_incorrect(ii, 1) = mean(rtEstIncorrect_individualTrial{ii});
    bootstat = bootstrp(n_boot_sample, @mean, rtEstIncorrect_individualTrial{ii});
    sem = std(bootstat);
    rt_incorrect(ii, 2) = rt_incorrect(ii, 1) - sem;
    rt_incorrect(ii, 3) = rt_incorrect(ii, 1) + sem;    
%     rt_incorrect(ii, 2:3) = bootci(n_boot_sample, @mean, rtEstIncorrect_individualTrial{ii});    
end
colorIndex = [1 0 0; 0 1 0];
errorBarGraph([rt_correct(:, 1)'; rt_incorrect(:, 1)'], [rt_correct(:, 2)'; rt_incorrect(:, 2)'],...
                    [rt_correct(:, 3)'; rt_incorrect(:, 3)'], colorIndex)
xlabel('Trial type')
ylabel('Reaction time (s)')
ylim([2.4 3.2])

rtEstCorrect_all = [rtEstCorrect_individualTrial{1} rtEstCorrect_individualTrial{2}];
rtEstIncorrect_all = [rtEstIncorrect_individualTrial{1} rtEstIncorrect_individualTrial{2}];
rtEstCorrect_mean = mean(rtEstCorrect_all);
rtEstIncorrect_mean = mean(rtEstIncorrect_all);
fprintf('Correct: %4.3f, Incorrect: %4.3f, Diff: %4.3f \n', rtEstCorrect_mean, rtEstIncorrect_mean,...
                    rtEstIncorrect_mean - rtEstCorrect_mean)

rtEstLowNoise_all = [rtEstCorrect_individualTrial{1} rtEstIncorrect_individualTrial{1}];
rtEstHighNoise_all = [rtEstCorrect_individualTrial{2} rtEstIncorrect_individualTrial{2}];
rtEstLowNoise_mean = mean(rtEstLowNoise_all);
rtEstHighNoise_mean = mean(rtEstHighNoise_all);
fprintf('Low noise: %4.3f, High noise: %4.3f, Diff: %4.3f \n', rtEstLowNoise_mean, rtEstHighNoise_mean,...
                    rtEstHighNoise_mean - rtEstLowNoise_mean)
                
%% Plot reaction time as a function of estimate (collapsed across stimulus)
% Correct trial
binCenter = -40:5:40;
rtDecisionCorrect = NaN(length(stdNoiseLevel), length(binCenter));
rtEstimateCorrect = NaN(length(stdNoiseLevel), length(binCenter));
for ii = 1 : length(stdNoiseLevel)
    indExclude = [binaryDecision{ii,:}] ~= [binaryFeedback{ii,:}];
    indExclude = indExclude(:);
    tempRTDecision = [decisionTime{ii, :}];
    tempRTDecision = tempRTDecision(:);
    tempRTDecision(indExclude) = [];
    tempRTEstimate = [estimateTime{ii, :}];
    tempRTEstimate = tempRTEstimate(:);
    tempRTEstimate(indExclude) = [];
    estimateData = [angleDiffEstimate{ii, :}];
    estimateData = estimateData(:);
    estimateData(indExclude) = [];
    rtDecisionCorrect(ii, :) = computeHeatMap(estimateData, tempRTDecision, binCenter);
    rtEstimateCorrect(ii, :) = computeHeatMap(estimateData, tempRTEstimate, binCenter);
end

figure
hold on
subplot(2, 2, 1)
plot(binCenter, rtDecisionCorrect, 'o-')
xlabel('Estimated orientation')
ylabel('Reaction time')
title('Decision')
legend('Low noise', 'High noise')
xlim([min(binCenter) max(binCenter)])

subplot(2, 2, 2)
plot(binCenter, rtEstimateCorrect, 'o-')
xlabel('Estimated orientation')
ylabel('Reaction time')
title('Estimation')
xlim([min(binCenter) max(binCenter)])

% Incorrect trial
binCenter = -40:5:40;
rtDecisionIncorrect = NaN(length(stdNoiseLevel), length(binCenter));
rtEstimateIncorrect = NaN(length(stdNoiseLevel), length(binCenter));
for ii = 1 : length(stdNoiseLevel)
    indExclude = [binaryDecision{ii,:}] == [binaryFeedback{ii,:}];
    indExclude = indExclude(:);
    tempRTDecision = [decisionTime{ii, :}];
    tempRTDecision = tempRTDecision(:);
    tempRTDecision(indExclude) = [];
    tempRTEstimate = [estimateTime{ii, :}];
    tempRTEstimate = tempRTEstimate(:);
    tempRTEstimate(indExclude) = [];
    estimateData = [angleDiffEstimate{ii, :}];
    estimateData = estimateData(:);
    estimateData(indExclude) = [];
    rtDecisionIncorrect(ii, :) = computeHeatMap(estimateData, tempRTDecision, binCenter);
    rtEstimateIncorrect(ii, :) = computeHeatMap(estimateData, tempRTEstimate, binCenter);
end

subplot(2, 2, 3)
plot(binCenter, rtDecisionIncorrect, 'o-')
xlabel('Estimated orientation')
ylabel('Reaction time')
title('Decision')
xlim([min(binCenter) max(binCenter)])

subplot(2, 2, 4)
plot(binCenter, rtEstimateIncorrect, 'o-')
xlabel('Estimated orientation')
ylabel('Reaction time')
title('Estimation')
xlim([min(binCenter) max(binCenter)])

% Plot the average for each panel in the figure above
figure
hold on
subplot(2, 2, 1)
bar([1 2], nanmedian(rtDecisionCorrect, 2))
xlabel('Stimulus noise')
ylabel('Reaction time')
title('Decision')
set(gca, 'XTick', 1:2, 'XTickLabel', {'Low noise', 'High noise'})
ylim([0 yMax])

subplot(2, 2, 2)
bar([1 2], nanmedian(rtEstimateCorrect, 2))
xlabel('Stimulus noise')
ylabel('Reaction time')
title('Estimation')
set(gca, 'XTick', 1:2, 'XTickLabel', {'Low noise', 'High noise'})
ylim([0 yMax])


% Incorrect trial
subplot(2, 2, 3)
bar([1 2], nanmedian(rtDecisionIncorrect, 2))
xlabel('Stimulus noise')
ylabel('Reaction time')
title('Decision')
set(gca, 'XTick', 1:2, 'XTickLabel', {'Low noise', 'High noise'})
ylim([0 yMax])

subplot(2, 2, 4)
bar([1 2], nanmedian(rtEstimateIncorrect, 2))
xlabel('Stimulus noise')
ylabel('Reaction time')
title('Estimation')
set(gca, 'XTick', 1:2, 'XTickLabel', {'Low noise', 'High noise'})
ylim([0 yMax])