%%%%%%%% Data analysis of conditional decision making experiment %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
% 

subjectID = 'average';
experimentNumber = 1:5;
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
%  Column 1: true angle difference (population)
%  Column 2: bar presentation time
%  Column 3: SOA
%  Column 4: bar noise level
%  Column 5: bar reference angle
%  Column 6: subject's categorical decision
%  Column 7: subject's estimate
%  Column 8: true category
%  Column 9: true angle difference (sample) 
stdNoiseLevel = unique(dataAll(:,4));
angleDiff = (unique(dataAll(:,1)))';
dataZeroDiff = dataAll(dataAll(:,1) == 0,:);
angleEstimate = dataAll(:, 7);
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
    indicatorConsistent = sign(angleDiffEst) == binaryFeedbackAll;
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
histEstimate = cell(1, length(stdNoiseLevel));
for ii = 1 : length(stdNoiseLevel)
    histEstimate_temp = zeros(length(yAxis), length(xAxis));
    for jj = 1 :length(angleDiff)
        tempAngleDiffEst = angleDiffEstimate{ii,jj};
        tempBinaryDecision = binaryDecision{ii,jj};
        tempFeedBack = binaryFeedback{ii,jj};
        tempAngleDiffEst(tempBinaryDecision ~= tempFeedBack) = NaN;
        histEstimate_temp(:, jj) = hist(tempAngleDiffEst, yAxis);        
    end
    histEstimate_temp = max(histEstimate_temp(:)) - histEstimate_temp;
    histEstimate{ii} = histEstimate_temp;
    subplot(2, length(stdNoiseLevel), ii)
    hold on
    imagesc(histEstimate_temp)
    colormap gray
    axis xy
    plot([1 length(xAxis)], [find(yAxis == xAxis(1)) find(yAxis == xAxis(end))], 'b:', 'LineWidth', 1.5)
    plot([1 length(xAxis)], [round(length(yAxis)/2) round(length(yAxis)/2)],  'k:', 'LineWidth', 1.5)
    plot([round(length(xAxis)/2) round(length(xAxis)/2)], [1 length(yAxis)],  'k:', 'LineWidth', 1.5)    
    xDisplay = -22:11:22;
    yDisplay = round(linspace(yAxis(1) ,yAxis(end),5));
    set(gca, 'ylim', [1 length(yAxis)], 'xlim', [0 length(xAxis)+1], ...
        'XTick', linspace(0, length(xAxis)+1, 5) , 'XTickLabel', num2cell(xDisplay),...
        'YTick', find(ismember(yAxis, yDisplay)), 'YTickLabel', num2cell(yDisplay), ...
        'FontSize', 20)
end

% Incorrect trials
histEstimate = cell(1, length(stdNoiseLevel));
for ii = 1 : length(stdNoiseLevel)
    histEstimate_temp = zeros(length(yAxis), length(xAxis));
    for jj = 1 :length(angleDiff)
        tempAngleDiffEst = angleDiffEstimate{ii,jj};
        tempBinaryDecision = binaryDecision{ii,jj};
        tempFeedBack = binaryFeedback{ii,jj};
        tempAngleDiffEst(tempBinaryDecision == tempFeedBack) = NaN;
        histEstimate_temp(:, jj) = hist(tempAngleDiffEst, yAxis);        
    end
    histEstimate_temp = max(histEstimate_temp(:)) - histEstimate_temp;
    histEstimate{ii} = histEstimate_temp;
    subplot(2, length(stdNoiseLevel), length(stdNoiseLevel)+ii)
    hold on
    imagesc(histEstimate_temp)
    colormap gray
    axis xy
    plot([1 length(xAxis)], [find(yAxis == xAxis(1)) find(yAxis == xAxis(end))], 'b:', 'LineWidth', 1.5)
    plot([1 length(xAxis)], [round(length(yAxis)/2) round(length(yAxis)/2)],  'k:', 'LineWidth', 1.5)
    plot([round(length(xAxis)/2) round(length(xAxis)/2)], [1 length(yAxis)],  'k:', 'LineWidth', 1.5)    
    xDisplay = -22:11:22;
    yDisplay = round(linspace(yAxis(1) ,yAxis(end),5));
    set(gca, 'ylim', [1 length(yAxis)], 'xlim', [0 length(xAxis)+1], ...
        'XTick', linspace(0, length(xAxis)+1, 5) , 'XTickLabel', num2cell(xDisplay),...
        'YTick', find(ismember(yAxis, yDisplay)), 'YTickLabel', num2cell(yDisplay), ...
        'FontSize', 12)
    xlabel('True orientation (degree)')
    ylabel('Estimated orientation (degree)')    
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
maxYplot = max(params.barAngleDiff);
maxXplot = max(params.barAngleDiff);
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
        tempBias1 = tempAngleDiffEstimateCW_Correct(:,rangeCollapse:end) - repmat(angleDiff(rangeCollapse:end),size(tempAngleDiffEstimateCW_Correct,1),1);
        tempBias1 = tempBias1 - repmat(biasMotor(rangeCollapse:end), size(tempBias1,1),1);
        
        tempAngleDiffEstimateCCW_Correct = NaN(size(tempAngleDiffEstimate));
        tempAngleDiffEstimateCCW_Correct(tempBinaryDecision == -1 & indexCorrect) = tempAngleDiffEstimate(tempBinaryDecision == -1 & indexCorrect);
        tempAngleDiffEstimateCCWAve = squeeze(nanmean(tempAngleDiffEstimateCCW_Correct,1));
        tempAngleDiffEstimateCCWAve(rangeCollapse+1:end) = NaN;
        tempAngleDiffEstimateCCWSEM = squeeze(nanstd(tempAngleDiffEstimateCCW_Correct,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCCW_Correct),1));
        tempBias2 = tempAngleDiffEstimateCCW_Correct(:,1:rangeCollapse) - repmat(angleDiff(1:rangeCollapse),size(tempAngleDiffEstimateCCW_Correct,1),1);
        tempBias2 = tempBias2 - repmat(biasMotor(1:rangeCollapse), size(tempBias2,1),1);
        tempBias = [tempBias1;-tempBias2(:,sort(1:size(tempBias2,2),'descend'))];
        tempBiasMean = nanmean(tempBias,1);
        bias = [ones(1,boxcarLength-1)*tempBiasMean(1) tempBiasMean ones(1,boxcarLength-1)*tempBiasMean(end)];
        smoothBias = conv(bias, ones(1,boxcarLength)/boxcarLength, 'full');
        biasMeanCorrect(ii,:) = smoothBias(boxcarLength:boxcarLength+rangeCollapse-1);
        biasSEMCorrect(ii,:) = squeeze(nanstd(tempBias,1))./sqrt(nansum(~isnan(tempBias),1)); 
        
        tempEstimateCollapse = [-fliplr(tempAngleDiffEstimateCCW_Correct); tempAngleDiffEstimateCW_Correct];
        for jj = 1 : rangeCollapse
            estimateCorrectData{ii,jj} = tempEstimateCollapse(:,jj+rangeCollapse-1);
        end
        
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
        xlim([-maxXplot maxXplot])
        ylim([-maxYplot maxYplot])
        box off
        text(6, -5, '+/-1SEM', 'FontSize', 15)
        
        % Calculate bias for incorrect trials
        tempAngleDiffEstimateCW_Incorrect = NaN(size(tempAngleDiffEstimate));
        tempAngleDiffEstimateCW_Incorrect(tempBinaryDecision == -1 & ~indexCorrect) = tempAngleDiffEstimate(tempBinaryDecision == -1 & ~indexCorrect);
        tempAngleDiffEstimateCWAve = squeeze(nanmean(tempAngleDiffEstimateCW_Incorrect,1));
        tempAngleDiffEstimateCWAve(1:rangeCollapse-1) = NaN;
        tempAngleDiffEstimateCWSEM = squeeze(nanstd(tempAngleDiffEstimateCW_Incorrect,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCW_Incorrect),1));
        tempBias1 = tempAngleDiffEstimateCW_Incorrect(:,rangeCollapse:end) - repmat(angleDiff(rangeCollapse:end),size(tempAngleDiffEstimateCW_Incorrect,1),1);
        tempBias1 = tempBias1 - repmat(biasMotor(rangeCollapse:end), size(tempBias1,1),1);
        
        tempAngleDiffEstimateCCW_Incorrect = NaN(size(tempAngleDiffEstimate));
        tempAngleDiffEstimateCCW_Incorrect(tempBinaryDecision == 1 & ~indexCorrect) = tempAngleDiffEstimate(tempBinaryDecision == 1 & ~indexCorrect);
        tempAngleDiffEstimateCCWAve = squeeze(nanmean(tempAngleDiffEstimateCCW_Incorrect,1));
        tempAngleDiffEstimateCCWAve(rangeCollapse+1:end) = NaN;
        tempAngleDiffEstimateCCWSEM = squeeze(nanstd(tempAngleDiffEstimateCCW_Incorrect,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCCW_Incorrect),1));
        tempBias2 = tempAngleDiffEstimateCCW_Incorrect(:,1:rangeCollapse) - repmat(angleDiff(1:rangeCollapse),size(tempAngleDiffEstimateCCW_Incorrect,1),1);
        tempBias2 = tempBias2 - repmat(biasMotor(1:rangeCollapse), size(tempBias2,1),1);
        tempBias = [tempBias1;-tempBias2(:,sort(1:size(tempBias2,2),'descend'))];
        tempBiasMean = nanmean(tempBias,1);
        bias = [ones(1,boxcarLength-1)*tempBiasMean(1) tempBiasMean ones(1,boxcarLength-1)*tempBiasMean(end)];
        smoothBias = conv(bias, ones(1,boxcarLength)/boxcarLength, 'full');
        biasMeanIncorrect(ii,:) = smoothBias(boxcarLength:boxcarLength+rangeCollapse-1);
        biasSEMIncorrect(ii,:) = squeeze(nanstd(tempBias,1))./sqrt(nansum(~isnan(tempBias),1));
        
        tempEstimateCollapse = [-fliplr(tempAngleDiffEstimateCCW_Incorrect); tempAngleDiffEstimateCW_Incorrect];
        for jj = 1 : rangeCollapse
            estimateIncorrectData{ii,jj} = tempEstimateCollapse(:,jj+rangeCollapse-1);
        end   
        
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
        xlim([-maxXplot maxXplot])
        ylim([-maxYplot maxYplot])
        box off
        text(6, -5, '+/-1SEM', 'FontSize', 15)        
    end
    legend(hLegend, legendName, 'Location', 'NorthWest')
    
    % Plot the bias 
    hBiasCorrect = figure;
    figPos = [400, -40, 800, 800];
    set(hBiasCorrect,'Units','pixels','Position',figPos)
    hold on
    for ii = 1 : length(stdNoiseLevel)
        hold on
        set(gca,'FontSize',fontSize)
        hShade = shadedErrorBar(angleDiff(rangeCollapse:end), biasMeanCorrect(ii,:), biasSEMCorrect(ii,:),... 
                             {'Color', colorName{ii}, 'LineWidth', lineWidth});        
        hLegend(ii) = hShade.mainLine;
        legendName{ii} = ['noise = ' num2str(stdNoiseLevel(ii))];
        xlabel('Absolute true angle (degree)')
        ylabel('Bias (degree)') 
%         title(['Subject ' upper(subjectID) ]);
        xlim([0 maxXplot])
        if strcmp(experiment, 'ControlReplication')
            ylim([-4 20])
            set(gca, 'XTick', [-4:4:20])
        else
            ylim([-10 15])
            set(gca, 'YTick', [-10:5:15])            
        end
        plot([0 maxXplot], [0 0], 'k--', 'LineWidth', lineWidth)
    end
    legend(hLegend, legendName, 'Location', 'NorthEast')
    
    hBiasError = figure;
    figPos = [400, -40, 800, 800];
    set(hBiasError,'Units','pixels','Position',figPos)
    hold on
    for ii = 1 : length(stdNoiseLevel)
        hold on
        set(gca,'FontSize',fontSize)
        hShade = shadedErrorBar(angleDiff(rangeCollapse:end), biasMeanIncorrect(ii,:), biasSEMIncorrect(ii,:),... 
                             {'Color', colorName{ii}, 'LineWidth', lineWidth});        
        hLegend(ii) = hShade.mainLine;
        legendName{ii} = ['noise = ' num2str(stdNoiseLevel(ii))];
        xlabel('Absolute true angle (degree)')
        ylabel('Bias (degree)') 
%         title(['Subject ' upper(subjectID) ]);
        xlim([0 maxXplot])
        if strcmp(experiment, 'ControlReplication')
            ylim([-4 20])
            set(gca, 'XTick', [-4:4:20])
        else
            ylim([-10 15])
            set(gca, 'YTick', [-10:5:15])            
        end
        plot([0 maxXplot], [0 0], 'k--', 'LineWidth', lineWidth)
    end
    legend(hLegend, legendName, 'Location', 'NorthEast')
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
    xlim([-23 23])
    ylim([-23 23])
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


%% Compute the average params value
subjID = {'AverageLapse' };
paramsAll = cell(1, length(subjID));
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
    paramsAll{kk} = paramsAll;
    paramsMean(kk, :) = mean(paramsAll{kk}, 1);    
    
    plot(kk*ones(length(paramsAll(:, 6)), 1), paramsAll(:, 6), 'o');    
    plot(kk, mean(paramsAll(:, 6)), '*', 'MarkerSize', 10);  
    fprintf('%9.2f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f\r\n', ...
        paramsMean(1), paramsMean(2), paramsMean(3), paramsMean(4), paramsMean(5), paramsMean(6), paramsMean(7), paramsMean(8));        
end
xlim([0 length(subjID)+2])

%% Display the fit parameters (9 params)
subjIDAll = {'5'};
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

%% Plot the model's parameters
% parameters: stdNoiseLevel(1:3)  LapseRate    Prior     MemNoise   Smoothness 

paramModel =       [2.6500    6.0895           0.0000     22.2852     1.6506   0.9414    2.0976;
                    3.0023    9.7384           0.0000     34.4053     0.0615   0.9480    3.1069;
                    4.6136   10.4165           0.0000     29.8375     0.1325   0.9940    3.8106;
                    7.7094   11.9114           0.0000     55.7419     0.0083   0.2850    3.8551;
                    5.1033   10.3703           0.0000     46.6421     4.7921   0.8187    3.3313;
                    5.2200	  9.5650           0          48.1703     0.055    0.9768    3.3234];

noiseSensoryExp1 = paramModel(:, 1:2);
noiseMemoryExp1 = paramModel(:, 5);
noiseAll = [noiseSensoryExp1 noiseMemoryExp1];
priorRange = paramModel(:, 4);

% Plot the parameters with bars
fontSize = 23;
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
hFig = figure;
hAx1 = subplot(1, 2, 1);
[~, hPanel] = errorbar_groups(noiseAll', zeros(size(noiseAll')), zeros(size(noiseAll')), ...
                'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx1);
set(hPanel, 'EdgeColor', 'none')
set(gca, 'FontSize', fontSize)

hAx2 = subplot(1, 2, 2);
[~, hPanel] = errorbar_groups(priorRange', zeros(size(priorRange')), zeros(size(priorRange')), ...
                'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx2);
set(hPanel, 'EdgeColor', 'none')
set(gca, 'FontSize', fontSize)
