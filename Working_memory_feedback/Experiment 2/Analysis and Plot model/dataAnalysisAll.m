%%%%%%%% Data analysis of conditional decision making experiment %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
% 

subjectID = 'ln';
experimentNumber = 1;
experimentType = 'MainExperiment';
experiment = 'Original';
session = 1;
dataAll = [];
fontSize = 20;
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

% % Remove the horizontal reference trials
% cutoffPoint = 50;
% indRefRemove = (dataAll(:,5) < cutoffPoint) | (dataAll(:,5) > 180 - cutoffPoint);
% angleDiffEst(indRefRemove) = NaN;

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


%% Density plots of estimates
figure
hold on

% Correct trials
maxEstimate = 40;
minEstimate = -40;
yAxis = minEstimate:maxEstimate;
xAxis = angleDiff;
lengthXaxis = length(angleDiff(1):angleDiff(end)) + diff(angleDiff(end-1:end)) -1;
histEstimate = cell(1, length(stdNoiseLevel));
for ii = 1 : length(stdNoiseLevel)
    histEstimate_temp = zeros(length(yAxis), lengthXaxis);
    counter = 1;
    for jj = 1 :length(angleDiff)
        tempAngleDiffEst = angleDiffEstimate{ii,jj};
        tempBinaryDecision = binaryDecision{ii,jj};
        tempFeedBack = binaryFeedback{ii,jj};
        tempAngleDiffEst(tempBinaryDecision ~= tempFeedBack) = NaN;
        if angleDiff(jj) < 0
            histEstimate_temp(:, counter:counter+1) = repmat((hist(tempAngleDiffEst, yAxis))', [1 2]); 
            counter = counter + 2;
        else
            histEstimate_temp(:, counter:counter+4) = repmat((hist(tempAngleDiffEst, yAxis))', [1 5]); 
            counter = counter + 5;            
        end
    end
    histEstimate_temp = max(histEstimate_temp(:)) - histEstimate_temp;
    histEstimate{ii} = histEstimate_temp;
    subplot(2, length(stdNoiseLevel), ii)
    hold on
    imagesc(histEstimate_temp)
    colormap gray
    axis xy
    plot([2 lengthXaxis-2], [find(yAxis == xAxis(1)) find(yAxis == xAxis(end))], 'b:', 'LineWidth', 1.5)
    plot([1 lengthXaxis], [round(length(yAxis)/2) round(length(yAxis)/2)],  'k:', 'LineWidth', 1.5)
    plot([15 15], [1 length(yAxis)],  'k:', 'LineWidth', 1.5)    
    xDisplay = -10:10:30;
    xTick = round((xDisplay + 12) * lengthXaxis / (angleDiff(end)-angleDiff(1)));
    xTick(xDisplay<=0) = xTick(xDisplay<=0) + 2;
    yDisplay = round(linspace(yAxis(1) ,yAxis(end),5));
    set(gca, 'ylim', [1 length(yAxis)], 'xlim', [1 lengthXaxis], ...
        'XTick', xTick , 'XTickLabel', num2cell(xDisplay),...
        'YTick', find(ismember(yAxis, yDisplay)), 'YTickLabel', num2cell(yDisplay), ...
        'FontSize', 20)
end
 
% Incorrect trials
maxEstimate = 40;
minEstimate = -40;
yAxis = minEstimate:maxEstimate;
xAxis = angleDiff;
lengthXaxis = length(angleDiff(1):angleDiff(end)) + diff(angleDiff(end-1:end)) -1;
histEstimate = cell(1, length(stdNoiseLevel));
for ii = 1 : length(stdNoiseLevel)
    histEstimate_temp = zeros(length(yAxis), lengthXaxis);
    counter = 1;
    for jj = 1 :length(angleDiff)
        tempAngleDiffEst = angleDiffEstimate{ii,jj};
        tempBinaryDecision = binaryDecision{ii,jj};
        tempFeedBack = binaryFeedback{ii,jj};
        tempAngleDiffEst(tempBinaryDecision == tempFeedBack) = NaN;
        if angleDiff(jj) < 0
            histEstimate_temp(:, counter:counter+1) = repmat((hist(tempAngleDiffEst, yAxis))', [1 2]); 
            counter = counter + 2;
        else
            histEstimate_temp(:, counter:counter+4) = repmat((hist(tempAngleDiffEst, yAxis))', [1 5]); 
            counter = counter + 5;            
        end
    end
    histEstimate_temp = max(histEstimate_temp(:)) - histEstimate_temp;
    histEstimate{ii} = histEstimate_temp;
    if strcmp(experimentType, 'MotorNoise')
        subplot(2, length(stdNoiseLevel), ii+1)
    else
        subplot(2, length(stdNoiseLevel), ii+2)
    end
    hold on
    imagesc(histEstimate_temp)
    colormap gray
    axis xy
    plot([2 lengthXaxis-2], [find(yAxis == xAxis(1)) find(yAxis == xAxis(end))], 'b:', 'LineWidth', 1.5)
    plot([1 lengthXaxis], [round(length(yAxis)/2) round(length(yAxis)/2)],  'k:', 'LineWidth', 1.5)
    plot([15 15], [1 length(yAxis)],  'k:', 'LineWidth', 1.5)    
    xDisplay = -10:10:30;
    xTick = round((xDisplay + 12) * lengthXaxis / (angleDiff(end)-angleDiff(1)));
    xTick(xDisplay<=0) = xTick(xDisplay<=0) + 2;
    yDisplay = round(linspace(yAxis(1) ,yAxis(end),5));
    set(gca, 'ylim', [1 length(yAxis)], 'xlim', [1 lengthXaxis], ...
        'XTick', xTick , 'XTickLabel', num2cell(xDisplay),...
        'YTick', find(ismember(yAxis, yDisplay)), 'YTickLabel', num2cell(yDisplay), ...
        'FontSize', 20)
end
% save('densityEstimateExp3', 'histEstimate', 'xAxis', 'yAxis')


%% Plot the smoothed raw data
% hScatter = figure;
% figPos = [0.1, 0.2, 0.8, 0.6];
% set(hScatter,'Units','normalized','Position',figPos)
% hold on
% maxYplot = 40;
% imageData = zeros(length(-maxYplot:maxYplot), length(-22:22), length(stdNoiseLevel));
% for ii = 1 : length(stdNoiseLevel)
%     for jj = 1 :length(angleDiff)
%         tempAngleDiffEst = angleDiffEstimate{ii,jj};
%         tempAngleDiffEst(abs(tempAngleDiffEst)>90) = tempAngleDiffEst(abs(tempAngleDiffEst)>90) ...
%                                 - 180 * sign(tempAngleDiffEst(abs(tempAngleDiffEst)>90));
%         tempBinaryDecision = binaryDecision{ii,jj};
%         tempFeedBack = binaryFeedback{ii,jj};
%         indexCorrect = tempBinaryDecision == tempFeedBack;
%         tempAngleDiffEst(abs(tempAngleDiffEst)>40 | indexCorrect) = [];                    
%         for kk = 1 : length(tempAngleDiffEst)
%             if ~isnan(tempAngleDiffEst(kk)) && (round(tempAngleDiffEst(kk))+1<maxYplot)
%             imageData(maxYplot-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) = ...
%                 imageData(maxYplot-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) + 1;
%             end
%         end
%     end
%     for jj = 1 :length(angleDiff)
%         imageData(:, 3*jj-2, ii) = imageData(:, 3*jj-1, ii);
%         imageData(:, 3*jj, ii) = imageData(:, 3*jj-1, ii);
%     end
%     myfilter = fspecial('gaussian', [2 1], 1);
%     smoothImage = imfilter(imageData(:,:,ii), myfilter, 'replicate');
%     imageData(:,:,ii) = round(smoothImage*255/max(smoothImage(:)));
%     
%     subplot(1,length(stdNoiseLevel),ii)
%     hold on
%     tempImage = uint8(imageData(:,:,ii));
%     [height, width] = size(tempImage);
%     widthPlot = round(1*width);
%     tempImage = imresize(tempImage, [height widthPlot], 'bicubic');
%     imshow(tempImage)
%     axis on
%     hold on
%     set(gca,'FontSize',fontSize)
%     set(gca, 'ylim', [1 height], 'xlim', [1 widthPlot], ...
%         'XTick', round(linspace(1,widthPlot,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
%         'YTick', round(linspace(1,height,5)), 'YTickLabel', num2cell(round(linspace(maxYplot,-maxYplot,5))))
%     plot([1 widthPlot], [round(length(-maxYplot:maxYplot)/2)+22 round(length(-maxYplot:maxYplot)/2)-22], '--b', 'LineWidth', 1.1)
%     if strcmp(experimentType, 'MainExperiment')
%         title(['Stimulus noise = ' num2str(stdNoiseLevel(ii))])
%     else
%         title('Motor noise')
%     end
%     xlabel('True angle (degree)')
%     ylabel('Estimated angle (degree)')
% end

%% Plot the average
legendName = cell(1,length(stdNoiseLevel));
colorName = {'g', 'r', 'b', 'cyan', 'magenta', 'y'};
maxYplot = 40;
maxXplot = max(params.barAngleDiff);
minYplot = -40;
minXplot = min(params.barAngleDiff);

indexColor = 1;
hLegend = NaN(1, length(stdNoiseLevel));
boxcarLength = 1;
if strcmp(experimentType, 'MainExperiment')
    % Plot true angle
    hAverageCorrect = figure;
    figPos = [400, 100, 700, 700];
    set(hAverageCorrect,'Units','pixels','Position',figPos)
    hold on
    set(gca,'FontSize',fontSize)
    hAverageIncorrect = figure;
    figPos = [400, 100, 700, 700];
    set(hAverageIncorrect,'Units','pixels','Position',figPos)
    hold on
    set(gca,'FontSize',fontSize)    
    rangeCollapse = round(length(angleDiff)/2);
    estimateCorrectData = cell(length(stdNoiseLevel), rangeCollapse);
    estimateIncorrectData = cell(length(stdNoiseLevel), rangeCollapse);
    biasMeanCorrect = NaN(length(stdNoiseLevel), rangeCollapse);
    biasSEMCorrect = NaN(length(stdNoiseLevel), rangeCollapse);
    biasMeanIncorrect = NaN(length(stdNoiseLevel), rangeCollapse);
    biasSEMIncorrect = NaN(length(stdNoiseLevel), rangeCollapse);
    
    for ii = 1 : length(stdNoiseLevel)
        hold on
        set(gca,'FontSize',fontSize)
        tempAngleDiffEstimate = [angleDiffEstimate{ii,:}];
        tempIndicator1 = sign(tempAngleDiffEstimate) ~= repmat(sign(angleDiff), size(tempAngleDiffEstimate,1), 1);
        tempIndicator2 = abs(tempAngleDiffEstimate)>90;
        tempIndicator1(:,rangeCollapse) = 0;
        tempIndicator2(:,rangeCollapse) = 0;
        tempAngleDiffEstimate(tempIndicator1 & tempIndicator2) = tempAngleDiffEstimate(tempIndicator1 & tempIndicator2) - ...
                                            180 * sign(tempAngleDiffEstimate(tempIndicator1 & tempIndicator2));
        tempAngleDiffEstimate(tempIndicator1 & ~tempIndicator2) = NaN;
        tempBinaryDecision = [binaryDecision{ii,:}];
        tempFeedBack = [binaryFeedback{ii,:}];
        indexCorrect = tempBinaryDecision == tempFeedBack;
        
        % Calculate bias for correct trials
        tempAngleDiffEstimateCW_Correct = NaN(size(tempAngleDiffEstimate));
        tempAngleDiffEstimateCW_Correct(tempBinaryDecision == 1 & indexCorrect) = tempAngleDiffEstimate(tempBinaryDecision == 1 & indexCorrect);
        tempAngleDiffEstimateCWAve = squeeze(nanmean(tempAngleDiffEstimateCW_Correct,1));
        tempAngleDiffEstimateCWAve(1:rangeCollapse-1) = NaN;
        tempAngleDiffEstimateCWSEM = squeeze(nanstd(tempAngleDiffEstimateCW_Correct,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCW_Correct),1));
%         tempBias1 = tempAngleDiffEstimateCW_Correct(:,rangeCollapse:end) - repmat(angleDiff(rangeCollapse:end),size(tempAngleDiffEstimateCW_Correct,1),1);
%         tempBias1 = tempBias1 - repmat(biasMotor(rangeCollapse:end), size(tempBias1,1),1);
        
        tempAngleDiffEstimateCCW_Correct = NaN(size(tempAngleDiffEstimate));
        tempAngleDiffEstimateCCW_Correct(tempBinaryDecision == -1 & indexCorrect) = tempAngleDiffEstimate(tempBinaryDecision == -1 & indexCorrect);
        tempAngleDiffEstimateCCWAve = squeeze(nanmean(tempAngleDiffEstimateCCW_Correct,1));
        tempAngleDiffEstimateCCWAve(rangeCollapse+1:end) = NaN;
        tempAngleDiffEstimateCCWSEM = squeeze(nanstd(tempAngleDiffEstimateCCW_Correct,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCCW_Correct),1));
%         tempBias2 = tempAngleDiffEstimateCCW_Correct(:,1:rangeCollapse) - repmat(angleDiff(1:rangeCollapse),size(tempAngleDiffEstimateCCW_Correct,1),1);
%         tempBias2 = tempBias2 - repmat(biasMotor(1:rangeCollapse), size(tempBias2,1),1);
%         tempBias = [tempBias1;-tempBias2(:,sort(1:size(tempBias2,2),'descend'))];
%         tempBiasMean = nanmean(tempBias,1);
%         bias = [ones(1,boxcarLength-1)*tempBiasMean(1) tempBiasMean ones(1,boxcarLength-1)*tempBiasMean(end)];
%         smoothBias = conv(bias, ones(1,boxcarLength)/boxcarLength, 'full');
%         biasMeanCorrect(ii,:) = smoothBias(boxcarLength:boxcarLength+rangeCollapse-1);
%         biasSEMCorrect(ii,:) = squeeze(nanstd(tempBias,1))./sqrt(nansum(~isnan(tempBias),1)); 
        
        
        % Plot the estimates for correct trials
        figure(hAverageCorrect)
        hold on
        hShade = shadedErrorBar(angleDiff, tempAngleDiffEstimateCWAve, tempAngleDiffEstimateCWSEM,... 
                     {'Color', colorName{ii}, 'LineWidth', lineWidth});        
        hLegend(ii) = hShade.mainLine;
        shadedErrorBar(angleDiff, tempAngleDiffEstimateCCWAve, tempAngleDiffEstimateCCWSEM,... 
                             {'Color', colorName{ii}, 'LineWidth', lineWidth});        
        legendName{ii} = ['noise = ' num2str(stdNoiseLevel(ii))];
        xlabel('True angle (degree)')
        ylabel('Angle estimate (degree)')  
        title(['Subject ' upper(subjectID) ]);
        plot([angleDiff(1) angleDiff(end)],  [angleDiff(1) angleDiff(end)], '--k')
        plot([angleDiff(1) angleDiff(end)],  [0 0], '--k')
        xlim([minXplot maxXplot])
        ylim([minYplot maxYplot])
        box off
        text(6, -5, '+/-1SEM', 'FontSize', 15)
        
        % Calculate bias for incorrect trials
        tempAngleDiffEstimateCW_Incorrect = NaN(size(tempAngleDiffEstimate));
        tempAngleDiffEstimateCW_Incorrect(tempBinaryDecision == -1 & ~indexCorrect) = tempAngleDiffEstimate(tempBinaryDecision == -1 & ~indexCorrect);
        tempAngleDiffEstimateCWAve = squeeze(nanmean(tempAngleDiffEstimateCW_Incorrect,1));
        tempAngleDiffEstimateCWAve(1:rangeCollapse-1) = NaN;
        tempAngleDiffEstimateCWSEM = squeeze(nanstd(tempAngleDiffEstimateCW_Incorrect,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCW_Incorrect),1));
%         tempBias1 = tempAngleDiffEstimateCW_Incorrect(:,rangeCollapse:end) - repmat(angleDiff(rangeCollapse:end),size(tempAngleDiffEstimateCW_Incorrect,1),1);
%         tempBias1 = tempBias1 - repmat(biasMotor(rangeCollapse:end), size(tempBias1,1),1);
        
        tempAngleDiffEstimateCCW_Incorrect = NaN(size(tempAngleDiffEstimate));
        tempAngleDiffEstimateCCW_Incorrect(tempBinaryDecision == 1 & ~indexCorrect) = tempAngleDiffEstimate(tempBinaryDecision == 1 & ~indexCorrect);
        tempAngleDiffEstimateCCWAve = squeeze(nanmean(tempAngleDiffEstimateCCW_Incorrect,1));
        tempAngleDiffEstimateCCWAve(rangeCollapse+1:end) = NaN;
        tempAngleDiffEstimateCCWSEM = squeeze(nanstd(tempAngleDiffEstimateCCW_Incorrect,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCCW_Incorrect),1));
%         tempBias2 = tempAngleDiffEstimateCCW_Incorrect(:,1:rangeCollapse) - repmat(angleDiff(1:rangeCollapse),size(tempAngleDiffEstimateCCW_Incorrect,1),1);
%         tempBias2 = tempBias2 - repmat(biasMotor(1:rangeCollapse), size(tempBias2,1),1);
%         tempBias = [tempBias1;-tempBias2(:,sort(1:size(tempBias2,2),'descend'))];
%         tempBiasMean = nanmean(tempBias,1);
%         bias = [ones(1,boxcarLength-1)*tempBiasMean(1) tempBiasMean ones(1,boxcarLength-1)*tempBiasMean(end)];
%         smoothBias = conv(bias, ones(1,boxcarLength)/boxcarLength, 'full');
%         biasMeanIncorrect(ii,:) = smoothBias(boxcarLength:boxcarLength+rangeCollapse-1);
%         biasSEMIncorrect(ii,:) = squeeze(nanstd(tempBias,1))./sqrt(nansum(~isnan(tempBias),1));
        
        
        % Plot the estimates for incorrect trials
        figure(hAverageIncorrect)
        hold on
        hShade = shadedErrorBar(angleDiff, tempAngleDiffEstimateCWAve, tempAngleDiffEstimateCWSEM,... 
                     {'Color', colorName{ii}, 'LineWidth', lineWidth});        
        hLegend(ii) = hShade.mainLine;
        shadedErrorBar(angleDiff, tempAngleDiffEstimateCCWAve, tempAngleDiffEstimateCCWSEM,... 
                             {'Color', colorName{ii}, 'LineWidth', lineWidth});        
        legendName{ii} = ['noise = ' num2str(stdNoiseLevel(ii))];
        xlabel('True angle (degree)')
        ylabel('Angle estimate (degree)')  
        title(['Subject ' upper(subjectID) ]);
        plot([angleDiff(1) angleDiff(end)],  [angleDiff(1) angleDiff(end)], '--k')
        plot([angleDiff(1) angleDiff(end)],  [0 0], '--k')
        xlim([minXplot maxXplot])
        ylim([minYplot maxYplot])
        box off
        text(6, -5, '+/-1SEM', 'FontSize', 15)        
    end
    legend(hLegend, legendName, 'Location', 'NorthWest')
    
elseif strcmp(experimentType, 'TrainArray')
    figure
    stdMotor = NaN(1, length(angleDiff));
    meanMotor = NaN(1, length(angleDiff));
    set(gca,'FontSize',15)
    hold on
    for jj = 1 : length(angleDiff)
        tempEst = [angleDiffEstimate{:,jj}];
        tempEst = tempEst(:);
        tempEst(abs(tempEst)>90) = tempEst(abs(tempEst)>90) - 180 * sign(tempEst(abs(tempEst)>90));
        meanMotor(jj) = nanmean(tempEst);
        stdMotor(jj) = nanstd(tempEst);
        errorbar(angleDiff(jj), ...
                meanMotor(jj), 2*stdMotor(jj)/sqrt(sum(~isnan(tempEst))),...
                'bo', 'MarkerFaceColor','b', 'MarkerSize', 9)
    end
    biasMotor = meanMotor - angleDiff
    xlim([-23 23])
    ylim([-23 23])
    xlabel('\theta_{true} (degree)')
    ylabel('\theta_{estimate} (degree)')  
    title('Motor noise');
    plot([angleDiff(1) angleDiff(end)],  [angleDiff(1) angleDiff(end)], '--r', 'LineWidth', 2)
elseif strcmp(experimentType, 'MotorNoise')  
    figure
    if strcmp(experimentType, 'MotorNoise')
        stdNoiseLevel = 0;
    end
    stdMotor = NaN(1, length(angleDiff));
    meanMotor = NaN(1, length(angleDiff));
    set(gca,'FontSize',15)
    hold on
    for jj = 1 : length(angleDiff)
        tempEst = [angleDiffEstimate{:,jj}];
        tempEst = tempEst(:);
        tempEst(abs(tempEst)>90) = tempEst(abs(tempEst)>90) - 180 * sign(tempEst(abs(tempEst)>90));
        meanMotor(jj) = nanmean(tempEst);
        stdMotor(jj) = nanstd(tempEst);
        errorbar(angleDiff(jj), ...
                meanMotor(jj), 2*stdMotor(jj)/sqrt(sum(~isnan(tempEst))),...
                'bo', 'MarkerFaceColor','b', 'MarkerSize', 9)
    end
    biasMotor = meanMotor - angleDiff
    xlim([min(angleDiff) max(angleDiff)])
    ylim([min(angleDiff) max(angleDiff)])
    xlabel('\theta_{true} (degree)')
    ylabel('\theta_{estimate} (degree)')  
    title('Motor noise');
    plot([angleDiff(1) angleDiff(end)],  [angleDiff(1) angleDiff(end)], '--r', 'LineWidth', 2)
end

%% Calculate percent clockwise
if ~strcmp(experimentType, 'MotorNoise')  
    percentCW = NaN(length(stdNoiseLevel),length(params.barAngleDiff));
    nTrialsPerCondition = NaN(length(stdNoiseLevel),length(params.barAngleDiff));
    binaryDecisionAll = dataAll(:, 6);
    for ii = 1 : length(stdNoiseLevel)
        for jj = 1 : length(percentCW)
            indexAll = dataAll(:,1) == params.barAngleDiff(jj) ...
                        & dataAll(:,4) == stdNoiseLevel(ii);
            indexCW = dataAll(:,1) == params.barAngleDiff(jj) ...
                        & dataAll(:,4) == stdNoiseLevel(ii) ...
                        & binaryDecisionAll == 1;
            nTrialsPerCondition(ii,jj) = sum(~isnan(binaryDecisionAll(indexAll)));
            percentCW(ii,jj) = 100*nansum(indexCW)/nTrialsPerCondition(ii,jj);
        end
    end
    
    % Fit using cumulative Gaussian and plot the psychometric curve
    lineWidth = 3;
    hPsychCurve = figure;
    set(gca, 'FontSize', fontSize)
    hold on
    angleDiff = params.barAngleDiff;
    legendName = cell(1,length(stdNoiseLevel));
    colorName = {'g','r', 'b', 'cyan', 'magenta', 'y'};
    fixedPSE = 0;
    zeroLapse = 1;
    hLegend = NaN(1, length(stdNoiseLevel));
    deltaThreshold = NaN(1, length(stdNoiseLevel));
    disp('******* Cumulative Gaussian fit *******')
    [fitParameterCumGauss, negLogLHCumGauss, hLegend] = CumGaussFit(percentCW, nTrialsPerCondition, angleDiff, hPsychCurve, stdNoiseLevel, colorName, hLegend, fixedPSE);    
end

%% Plot the performance of subject
% if ~strcmp(experimentType, 'MainExperiment')
%     if strcmp(experimentType, 'PerceptNoise')
%         maxScorePerTrial = 1;
%     elseif strcmp(experimentType, 'MotorNoise') || strcmp(experimentType, 'TrainArray')
%         a = params.costFunctionParam(1);
%         b = params.costFunctionParam(2);
%         maxScorePerTrial = a/b;
%     end
%     scoreTemporal = 100 * params.score(2:end) / (maxScorePerTrial*params.nTrialDisplayScore);
%     scoreTotal = mean(scoreTemporal(scoreTemporal~=0));
%     figure
%     hold on
%     set(gca, 'Fontsize', fontSize)
%     plot(smooth(scoreTemporal,5), 'o-')
%     xlabel('Time')
%     ylabel('Normalized score')
%     title(['Average score = ' num2str(scoreTotal)])
% else 
%     % Calculate the median absolute deviation (median(|Xi - median(X)|) of estimate
%     madAll = [];
%     lengthRunningWindow = 100; % number of trials
%     lengthMadCumulative = round(size(dataAll, 1) / 100);
%     madCumulative = NaN(1, lengthMadCumulative);
%     for kk = 1 : lengthMadCumulative
%         madAll = [];
%         indEnd = kk*lengthRunningWindow;
%         if indEnd > size(dataAll, 1)
%             indEnd = size(dataAll, 1);
%         end
%         indSelect = params.trialOrder(1:indEnd);
%         dataSelect = dataAll(indSelect, :);
%         tempAngleReference1 = dataSelect(:,5);
%         angleEstimate = dataSelect(:, 7);
%         for ii = 1 : length(stdNoiseLevel)
%             for jj = 1 : length(angleDiff)
%                 tempAngleReference1 = tempAngleReference1(dataSelect(:,4)==stdNoiseLevel(ii) ...
%                                                     & dataSelect(:,1)==angleDiff(jj), :);
%                 tempAngleEstimate = angleEstimate(dataSelect(:,4)==stdNoiseLevel(ii) ...
%                                                     & dataSelect(:,1)==angleDiff(jj), :); 
%                 tempEst = tempAngleReference1-tempAngleEstimate;
%                 tempIndicator1 = sign(tempEst) ~= sign(angleDiff(jj));
%                 tempIndicator2 = abs(tempEst)>90;
%                 tempIndicator1(:,rangeCollapse) = 0;
%                 tempIndicator2(:,rangeCollapse) = 0;
%                 tempEst(tempIndicator1 & tempIndicator2) = tempEst(tempIndicator1 & tempIndicator2) ...
%                                     - 180 * sign(tempEst(tempIndicator1 & tempIndicator2));
%                 tempEst(tempIndicator1 & ~tempIndicator2) = NaN;
%                 tempEst(isnan(tempEst)) = [];
%                 if length(tempEst) > 1
%                     absDeviation = abs(tempEst - nanmedian(tempEst));
%                     madAll = [madAll absDeviation'];
%                 end
%             end
%         end
%         madCumulative(kk) = median(madAll);
%     end
    % Decision task
%     maxScorePerTrial = 1;
%     scoreTemporal = 100 * params.score(1,2:end) / (maxScorePerTrial*params.nTrialDisplayScore);
%     scoreTotal = mean(scoreTemporal(scoreTemporal~=0));
%     figure
%     hold on
%     set(gca, 'Fontsize', fontSize)
%     plot(smooth(scoreTemporal,1), 'o-')
%     xlabel('Time')
%     ylabel('Normalized score')
%     title(['Average score = ' num2str(scoreTotal)])
    
%     % Estimation task
%     figure;
%     plot(params.score(2,:), 'o-')
%     xlabel('Time')
%     ylabel('Score')
% end


%% Plot the response-divided version of smoothed raw data (red for CW and green for CCW)
hScatterColor = figure;
figPos = [0.1, 0.2, 0.8, 0.6];
set(hScatterColor,'Units','normalized','Position',figPos)
hold on
maxYplotColor = 60;
imageDataCW = zeros(length(-maxYplotColor:maxYplotColor), length(-22:22), length(stdNoiseLevel));
imageDataCCW = zeros(length(-maxYplotColor:maxYplotColor), length(-22:22), length(stdNoiseLevel));
imageDataBlue = zeros(length(-maxYplotColor:maxYplotColor), length(-22:22));
for ii = 1 : length(stdNoiseLevel)
    for jj = 1 :length(angleDiff)
        tempAngleDiffEst = angleDiffEstimate{ii,jj};
        tempAngleDiffEst(abs(tempAngleDiffEst)>90) = tempAngleDiffEst(abs(tempAngleDiffEst)>90) ...
                                - 180 * sign(tempAngleDiffEst(abs(tempAngleDiffEst)>90));
        tempBinaryDecision = binaryDecision{ii,jj};
        for kk = 1 : length(tempAngleDiffEst)
            if ~isnan(tempAngleDiffEst(kk)) && (round(tempAngleDiffEst(kk))+1<maxYplot)
            if tempBinaryDecision(kk) == -1;
                imageDataCCW(maxYplotColor-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) = ...
                    imageDataCCW(maxYplotColor-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) + 1;
            elseif tempBinaryDecision(kk) == 1;
                imageDataCW(maxYplotColor-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) = ...
                    imageDataCW(maxYplotColor-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) + 1;
            end
            end
        end
    end
    imageData = cat(3, imageDataCW(:,:,ii), imageDataCCW(:,:,ii), imageDataBlue);
    myfilter = fspecial('gaussian', [3 2], 1);
    smoothImage = imfilter(imageData, myfilter, 'replicate');
    smoothImage(:,:,1) = smoothImage(:,:,1)*255/max(max(smoothImage(:,:,1)));
    smoothImage(:,:,2) = smoothImage(:,:,2)*255/max(max(smoothImage(:,:,2)));    
    imageData = round(smoothImage);
    
    subplot_tight(1,length(stdNoiseLevel),ii,[0.03 0.03])
    hold on
    tempImage = uint8(imageData);
    [height, width, channel] = size(tempImage);
    widthPlot = 2*width;
    tempImage = imresize(tempImage, [height widthPlot], 'bicubic');
    imshow(tempImage)
    [height,~,~] = size(tempImage);
    axis on
    hold on
    set(gca,'FontSize',fontSize)
    set(gca, 'ylim', [1 height], 'xlim', [1 widthPlot], ...
        'XTick', round(linspace(3,widthPlot-2,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
        'YTick', round(linspace(1,height,7)), 'YTickLabel', num2cell(round(linspace(maxYplotColor,-maxYplotColor,7))))
    if ii==2 || ii==3
        set(gca, 'YTickLabel',{''})
    end
    plot([1 widthPlot], [round(length(-maxYplotColor:maxYplotColor)/2)+22 round(length(-maxYplotColor:maxYplotColor)/2)-22], '--b', 'LineWidth', lineWidth)
%     plot([width+1 width+1], [1 2*maxYplot], '--b', 'LineWidth', 1.1)  
%     xlabel('True angle (degree)')
    if ii == 1
        ylabel('Estimated angle (degree)', 'FontSize', 23)
    end
    title(['Stimulus noise = ' num2str(stdNoiseLevel(ii))], 'FontSize', 23)
end
tightfig

%% Plot the histogram
% The histogram of bias 
biasHistFig = figure;
figPos = [0.1, 0.1, 0.8, 0.8];
set(biasHistFig,'Units','normalized','Position',figPos)
hold on
rangeAngle = 21;
for ii = 1 : length(stdNoiseLevel)
    tempBias = [];
    for jj = 1 %: length(angleDiff)
        if (angleDiff(jj)>=-rangeAngle) && (angleDiff(jj)<=rangeAngle)
            tempBias = [tempBias; biasEstimate{ii,jj}];
        end
    end 
    subplot(1, length(stdNoiseLevel), ii)
    hold on
    
    % Fit by MLE
    pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
                         p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);
    pStart = .5;
    muStart = [mean(tempBias(tempBias<0)) mean(tempBias(tempBias>0))];
    sigmaStart = [std(tempBias(tempBias<0)) std(tempBias(tempBias>0))];
    start = [pStart muStart sigmaStart];
    lb = [0 -Inf -Inf 0 0];
    ub = [1 Inf Inf Inf Inf];
    options = statset('MaxIter',1000, 'MaxFunEvals',2000);
    paramEsts = mle(tempBias, 'pdf',pdf_normmixture, 'start',start, ...
                              'lower',lb, 'upper',ub, 'options',options);
                          
    % Plot the histogram and fit curve
    [n, x] = hist(tempBias,20); 
    n = n / ((x(2)-x(1)) * length(tempBias));
    bar(x,n)
    xFit = linspace(min(x), max(x),500);
    pdfFit = pdf_normmixture(xFit,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
    [~, peakInd] = findpeaks(pdfFit, 'npeaks', 2);
    plot(xFit,pdfFit,'r-', 'LineWidth', 1.5)
    plot([0 0], [0 max(n)+0.005], 'k--', 'LineWidth', 1.1)
    for kk = 1 : length(peakInd)
        plot([xFit(peakInd(kk)) xFit(peakInd(kk))], [0 pdfFit(peakInd(kk))], 'g--', 'LineWidth', 1.1)
    end
    ylim([0 max(n)+0.005])
    xlabel('Bias (degree)')  
    title(['Stimulus noise = ' num2str(stdNoiseLevel(ii))])  
end

% The histogram of estimate 
estHistFig = figure;
figPos = [0.1, 0.1, 0.8, 0.8];
set(estHistFig,'Units','normalized','Position',figPos)
hold on
rangeAngle = 21;
for ii = 1 : length(stdNoiseLevel)
    tempBias = [];
    for jj = 1 %: length(angleDiff)
        if (angleDiff(jj)>=-rangeAngle) && (angleDiff(jj)<=rangeAngle)
            tempBias = [tempBias; biasEstimate{ii,jj}];
        end
    end 
    subplot(1, length(stdNoiseLevel), ii)
    hold on
    
    % Fit by MLE
    pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
                         p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);
    pStart = .5;
    muStart = [mean(tempBias(tempBias<0)) mean(tempBias(tempBias>0))];
    sigmaStart = [std(tempBias(tempBias<0)) std(tempBias(tempBias>0))];
    start = [pStart muStart sigmaStart];
    lb = [0 -Inf -Inf 0 0];
    ub = [1 Inf Inf Inf Inf];
    options = statset('MaxIter',1000, 'MaxFunEvals',2000);
    paramEsts = mle(tempBias, 'pdf',pdf_normmixture, 'start',start, ...
                              'lower',lb, 'upper',ub, 'options',options);
                          
    % Plot the histogram and fit curve
    [n, x] = hist(tempBias,20); 
    n = n / ((x(2)-x(1)) * length(tempBias));
    bar(x,n,'style', 'hist')
    xFit = linspace(min(x), max(x),500);
    pdfFit = pdf_normmixture(xFit,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
    [~, peakInd] = findpeaks(pdfFit, 'npeaks', 2);
    plot(xFit,pdfFit,'r-', 'LineWidth', 1.5)
    plot([0 0], [0 max(n)+0.005], 'k--', 'LineWidth', 1.1)
    for kk = 1 : length(peakInd)
        plot([xFit(peakInd(kk)) xFit(peakInd(kk))], [0 pdfFit(peakInd(kk))], 'g--', 'LineWidth', 1.1)
    end
    ylim([0 max(n)+0.005])
    xlabel('Bias (degree)')  
    title(['Stimulus noise = ' num2str(stdNoiseLevel(ii))])
  
end

%% Compute the average params value
subjID = {'AverageLapse' };
paramsMean = NaN(length(subjID), 8);
figure;
hold on
for kk = 1 : length(subjID);
    fileName = ['FitResult-' subjID{kk} '.txt'];
    fileID = fopen(fileName);
    myFile = textscan(fileID,'%s','delimiter','\n');
    myFile = myFile{1};
    saveNextLine = 0;
    counter = 1;
    paramsAll = NaN(30, 8);
    for ii = 1 : length(myFile)
        tempLine = myFile{ii};
        if saveNextLine
            paramsAll(counter, :) = str2num(tempLine);
            counter = counter + 1;
            saveNextLine = 0;
        else
            indMatch = strfind(tempLine,'Iteration');
            if ~isempty(indMatch)
                saveNextLine = 1;
            end
        end
    end
    paramsMean(kk, :) = mean(paramsAll{kk}, 1);    
    
    plot(kk*ones(length(paramsAll(:, 6)), 1), paramsAll(:, 6), 'o');    
    plot(kk, mean(paramsAll(:, 6)), '*', 'MarkerSize', 10);  
    fprintf('%9.2f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f\r\n', ...
        paramsMean(1), paramsMean(2), paramsMean(3), paramsMean(4), paramsMean(5), paramsMean(6), paramsMean(7), paramsMean(8));        
end
xlim([0 length(subjID)+2])

%% Display the fit parameters (9 params)
subjIDAll = {'LN'};
for kk = 1 : length(subjIDAll)
    subjID = subjIDAll{kk};
    fileName = ['FitResult-' subjID '.txt'];
    fileID = fopen(fileName);
    myFile = textscan(fileID,'%s','delimiter','\n');
    myFile = myFile{1};
    saveNextLine = 0;
    counter = 1;
    paramsAll = NaN(30, 9);
    for ii = 1 : length(myFile)
        tempLine = myFile{ii};
        if saveNextLine
            paramsAll(counter, :) = str2num(tempLine);
            counter = counter + 1;
            saveNextLine = 0;
        else
            indMatch = strfind(tempLine,'Iteration');
            if ~isempty(indMatch)
                saveNextLine = 1;
            end
        end
    end

    [~, indexSort] = sort(paramsAll(:, 1));
    paramsSort = paramsAll(indexSort, :);
    paramsSort = [paramsSort; mean(paramsSort(1:20, :), 1)];
    currentFolder = pwd;
    fileNameRoot = 'FitResult';
    fileNumber = GetNextDataFileNumber(currentFolder, '.txt');
    fileID = fopen([fileNameRoot '-' num2str(fileNumber) '.txt'],'w');
    for ii = 1 : size(paramsSort, 1)-1
        fprintf(fileID, '%9.2f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f %9.4f\r\n', ...
            paramsSort(ii, 1), paramsSort(ii, 2), paramsSort(ii, 3), paramsSort(ii, 4),...
            paramsSort(ii, 5), paramsSort(ii, 6), paramsSort(ii, 7), paramsSort(ii, 8), paramsSort(ii, 9));     
    end
    ii = ii + 1;
    fprintf(fileID, '\n');
    fprintf(fileID, '%9.2f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f %9.4f\r\n', ...
        paramsSort(ii, 1), paramsSort(ii, 2), paramsSort(ii, 3), paramsSort(ii, 4),...
        paramsSort(ii, 5), paramsSort(ii, 6), paramsSort(ii, 7), paramsSort(ii, 8), paramsSort(ii, 9));     

    fclose(fileID);
end

%% Display the fit parameters (10 params)
subjIDAll = {'1'};
for kk = 1 : length(subjIDAll)
    subjID = subjIDAll{kk};
    fileName = ['FitResult-' subjID '.txt'];
    fileID = fopen(fileName);
    myFile = textscan(fileID,'%s','delimiter','\n');    
    myFile = myFile{1};
    saveNextLine = 0;
    counter = 1;
    paramsAll = NaN(30, 10);
    for ii = 1 : length(myFile)
        tempLine = myFile{ii};
        if saveNextLine
            paramsAll(counter, :) = str2num(tempLine);
            counter = counter + 1;
            saveNextLine = 0;
        else
            indMatch = strfind(tempLine,'Iteration');
            if ~isempty(indMatch)
                saveNextLine = 1;
            end
        end
    end

    [~, indexSort] = sort(paramsAll(:, 1));
    paramsSort = paramsAll(indexSort, :);
    paramsSort = [paramsSort; mean(paramsSort(1:20, :), 1)];
    currentFolder = pwd;
    fileNameRoot = 'FitResult';
    fileNumber = GetNextDataFileNumber(currentFolder, '.txt');
    fileID = fopen([fileNameRoot '-' num2str(fileNumber) '.txt'],'w');
    for ii = 1 : size(paramsSort, 1)-1
        fprintf(fileID, '%9.2f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f %9.4f %9.4f\r\n', ...
            paramsSort(ii, 1), paramsSort(ii, 2), paramsSort(ii, 3), paramsSort(ii, 4),...
            paramsSort(ii, 5), paramsSort(ii, 6), paramsSort(ii, 7), paramsSort(ii, 8), paramsSort(ii, 9), paramsSort(ii, 10));     
    end
    ii = ii + 1;
    fprintf(fileID, '\n');
    fprintf(fileID, '%9.2f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f %9.4f %9.4f\r\n', ...
        paramsSort(ii, 1), paramsSort(ii, 2), paramsSort(ii, 3), paramsSort(ii, 4),...
        paramsSort(ii, 5), paramsSort(ii, 6), paramsSort(ii, 7), paramsSort(ii, 8), paramsSort(ii, 9), paramsSort(ii, 10));     

    fclose(fileID);
end

%% Plot the fit parameters
% parameters: stdNoiseLevel(1:3)  LapseRate    Prior     MemNoise   Smoothness 

paramModel =       [4.3425    6.2248           0.0000     19.1377    -9.6045   3.5283    2.0902    0.9927    0.6813;
                    8.7205    8.8878           0.0000     33.0312   -21.2586   1.2076    1.8928    0.9215    0.4348;
                    8.3099    8.7268           0.0000     13.9033   -12.6251   0.6127    2.7094    0.9468    0.5099;
                    6.3379    8.4823           0.0000     19.3091   -12.7690   0.9699    2.6041    0.5558    0.4201;
                    6.4379    9.9076           0.0000     32.5355   -17.6389   3.0964    1.5830    0.2067    0.4086;
                    9.7284   15.1847           0.0000     56.3034   -42.2707   0.4378    4.0136    0.9881    0.5523;
                    7.8510    9.8641           0.0000     22.8922   -18.0641   0.8543    3.9069    0.9646    0.4949;
                    6.0243    8.0650           0.0000     32.4649   -19.9032   5.3989    2.4021    0.9986    0.5303];


noiseSensoryExp1 = paramModel(:, 1:2);
noiseMemoryExp1 = paramModel(:, 6);
noiseAll = [noiseSensoryExp1 noiseMemoryExp1];
priorRange = abs(paramModel(:, 4:5));
pCW = paramModel(:, 9);

% Sensory + memory noise
fontSize = 23;
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
hFig = figure;
hAx1 = subplot(1, 3, 1);
[~, hPanel] = errorbar_groups(noiseAll', zeros(size(noiseAll')), zeros(size(noiseAll')), ...
                'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx1);
set(hPanel, 'EdgeColor', 'none')
set(gca, 'FontSize', fontSize)

% Prior range
colorName = {'Black', 'Gray'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
hAx2 = subplot(1, 3, 2);
hold on
[~, hPanel] = errorbar_groups(priorRange', zeros(size(priorRange')), zeros(size(priorRange')), ...
                'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx2);
subplot(1, 3, 2);
hold on            
plot([0 15], [30 30], '--')
plot([0 15], [12 12], '--')            
set(hPanel, 'EdgeColor', 'none')
set(gca, 'FontSize', fontSize)

% pCW
colorName = {'SlateGray'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
hAx3 = subplot(1, 3, 3);
hold on
[~, hPanel] = errorbar_groups(pCW', zeros(size(pCW')), zeros(size(pCW')), ...
                'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx3);
subplot(1, 3, 3);
hold on            
plot([0 15], [0.5 0.5], '--')
ylim([0 1])
set(hPanel, 'EdgeColor', 'none')
set(gca, 'FontSize', fontSize)