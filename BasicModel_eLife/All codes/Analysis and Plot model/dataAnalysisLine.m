%%%%%%%% Data analysis of conditional decision making experiment %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
% 

subjectID = 'average';
experimentNumber = 1:5;
experimentType = 'MainExperiment';
experiment = 'ControlReplication';
session = 1;
extractModelParam = 1;
dataAll = [];
fontSize = 20;
includeIncongruentTrials = 0; % 0: no incongruent trials
                              % 1: include incongruent trials
                              % 2: incongruent trials only
collapseScatterPlot = 0;                              
trialOrder = [];
indexTrialOrder = [];
maxTrialOrder = 0;
includeIncorrectTrial = 1;
earlyTrial = 0;
percentInclude = 0.33;

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
    trialOrder = [trialOrder params.trialOrder+maxTrialOrder];
    maxTrialOrder = max(trialOrder);
    indexTrialOrder = [indexTrialOrder 1:length(params.trialOrder)];
    tempTrialOrder = params.trialOrder;
%     if earlyTrial == 1
%         dataTotal(tempTrialOrder(length(dataTotal) - round(percentInclude * length(dataTotal)) : end), 7) = NaN;
%     else
%         dataTotal(tempTrialOrder(1 : round(percentInclude * length(dataTotal))), 7) = NaN;
%     end
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
if includeIncongruentTrials == 1
    dataName = 'Incongruent';
else
    dataName = '';
end
if strcmp(subjectID, 'average')
    subjectName = '';
else
    subjectName = upper(subjectID);
end
if strcmp(experiment, 'ControlReplication')
    estimateDataName = ['estimateDataExp' num2str(1) subjectName dataName];
    discriminationDataName = ['discriminationDataExp' num2str(1) subjectName dataName];    
elseif strcmp(experiment, 'Original')
    estimateDataName = ['estimateDataExp' num2str(2) subjectName dataName];
    discriminationDataName = ['discriminationDataExp' num2str(2) subjectName dataName];
elseif strcmp(experiment, 'DecisionGiven')
    estimateDataName = ['estimateDataExp' num2str(3) subjectName dataName];
end
paramsFit = [3.37      4.87      9.27           0.0053      22.33       1.43     0.99	  2.56];
stdNoiseLevelAll = paramsFit(:, 1:3);
priorRangeAll = paramsFit(:, 5);
smoothFactorAll = paramsFit(:, 7);
lapseRate = paramsFit(:, 4);
windowPrior = 1;

%% Extract data
%  Column 1: angle differences (theta1 - theta2)
%  Column 2: bar presentation time
%  Column 3: SOA
%  Column 4: bar Stimulus noise
%  Column 5: bar reference angle
%  Column 6: CW/CCW OR Red/Green
%  Column 7: estimated angle
%  Column 8: given CW/CCW
stdNoiseLevel = unique(dataAll(:,4));
angleDiff = (unique(dataAll(:,1)))';
angleEstimate = dataAll(:, 7);

%% Average across conditions
angleTrue = dataAll(:,5) - dataAll(:,1);
angleTrue(angleTrue>180 | angleTrue<0) = angleTrue(angleTrue>180 | angleTrue<0)...
                                            - 180 * sign(angleTrue(angleTrue>180 | angleTrue<0));
indexCorrect = find(abs(angleTrue-angleEstimate) > 70);
angleEstimate(indexCorrect) = angleEstimate(indexCorrect) + sign(angleTrue(indexCorrect)-angleEstimate(indexCorrect))*180;    
angleDiffEst = dataAll(:, 5) - angleEstimate;
indexAdjust = abs(angleDiffEst)>90;
angleDiffEst(indexAdjust) = angleDiffEst(indexAdjust) - 180 * sign(angleDiffEst(indexAdjust));
biasAll = angleDiffEst - dataAll(:, 1);
biasAll(abs(biasAll)>50) = NaN;

if strcmp(experiment, 'DecisionGiven')
    binaryDecisionAll = dataAll(:,8);
    indicatorConsistent = sign(angleDiffEst)==dataAll(:,8);
elseif strcmp(experimentType, 'MotorNoise') || strcmp(experiment, 'TrainArray')
    indicatorConsistent = ones(length(biasAll), 1);
    binaryDecisionAll = dataAll(:,6);    
else
    binaryDecisionAll = dataAll(:,6);
    indicatorConsistent = sign(angleDiffEst)==dataAll(:,6);    
end
sum(~indicatorConsistent(abs(angleDiffEst)>7))
100*nansum(~indicatorConsistent) / sum(~isnan(indicatorConsistent))
if includeIncongruentTrials == 1
    indicatorConsistentKeep = indicatorConsistent;
    indicatorConsistent = ones(length(biasAll), 1);
elseif includeIncongruentTrials == 2
    indicatorConsistentKeep = indicatorConsistent;
    indicatorConsistent = ~indicatorConsistent;    
end
angleDiffEstimate = cell(length(stdNoiseLevel), length(angleDiff));
biasEstimate = cell(length(stdNoiseLevel), length(angleDiff));
binaryDecisionSelf = cell(length(stdNoiseLevel), length(angleDiff));
binaryDecisionGiven = binaryDecisionSelf;
colorError = NaN(length(stdNoiseLevel), length(angleDiff));
for ii = 1 : length(stdNoiseLevel)
    for jj = 1 : length(angleDiff)
        indexSelect = dataAll(:,4)==stdNoiseLevel(ii) & dataAll(:,1)==angleDiff(jj);
        tempBias = biasAll(indexSelect, :);
        biasEstimate{ii,jj} = tempBias(:);                       
        indexConsistent = indicatorConsistent(indexSelect, :);                                
        tempEst = angleDiffEst(indexSelect);
        tempEst(~indexConsistent) = NaN;                                        
        angleDiffEstimate{ii,jj}  = tempEst;
        if strcmp(experiment, 'DecisionGiven')
            binaryDecisionSelf{ii,jj} = sign(angleDiffEstimate{ii,jj});
            binaryDecisionGiven{ii,jj} = binaryDecisionAll(indexSelect, :);     
            binaryDecision = binaryDecisionGiven; 
            tempColorCorrect = dataAll(indexSelect, 6);
            colorError(ii, jj) = 100 * sum(tempColorCorrect == 0) / sum(~isnan(tempColorCorrect));
        else
            tempDecision = binaryDecisionAll(indexSelect, :);
            tempDecision(~indexConsistent) = NaN;   
            binaryDecisionSelf{ii,jj} = tempDecision;
            binaryDecision = binaryDecisionSelf;                            
        end
    end
end

%% Plot the color error
% if strcmp(experiment, 'DecisionGiven')
%     figure
%     for ii = 1 : length(stdNoiseLevel)
%         subplot(1, length(stdNoiseLevel), ii)
%         plot(colorError(ii, :), 'o-', 'MarkerSize', 15)
%         ylim([0 5.5])
%         box off
%         set(gca, 'XTick', 1:length(angleDiff), 'XTickLabel', num2cell(angleDiff))
%         xlabel('Stimulus orientation')
%         ylabel('Color error (%)')
%         title(['Noise: ' num2str(stdNoiseLevel(ii)) ', Average error: ' num2str(mean(colorError(ii, :)))])
%     end
% end

%% Plot the histogram
% figure; 
% binCenter = linspace(-40, 40, 25);
% subplot(1, 3, 1)
% set(gca,'FontSize',fontSize)
% hist(angleDiffEst(dataAll(:,4)==stdNoiseLevel(1)), binCenter)
% xlim([-40 40])
% subplot(1, 3, 2)
% set(gca,'FontSize',fontSize)
% hist(angleDiffEst(dataAll(:,4)==stdNoiseLevel(2)), binCenter)
% xlim([-40 40])
% subplot(1, 3, 3)
% set(gca,'FontSize',fontSize)
% hist(angleDiffEst(dataAll(:,4)==stdNoiseLevel(3)), binCenter)
% xlim([-40 40])
% 
hScatter = figure;
figPos = [0.1, 0.2, 0.8, 0.6];
set(hScatter,'Units','normalized','Position',figPos)
hold on
maxYplot = 40;
histEstimate = zeros(length(-maxYplot:maxYplot), length(-22:22), length(stdNoiseLevel));
for ii = 1 : length(stdNoiseLevel)
    for jj = 1 :length(angleDiff)
        tempAngleDiffEst = angleDiffEstimate{ii,jj};
        tempAngleDiffEst(abs(tempAngleDiffEst)>35) = [];
        for kk = 1 : length(tempAngleDiffEst)
            if ~isnan(tempAngleDiffEst(kk)) && (round(tempAngleDiffEst(kk))+1<maxYplot)
            histEstimate(maxYplot-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) = ...
                histEstimate(maxYplot-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) + 1;
            end
        end
    end

    for jj = 1 :length(angleDiff)
        histEstimate(:, 3*jj-2, ii) = histEstimate(:, 3*jj-1, ii);
        histEstimate(:, 3*jj, ii) = histEstimate(:, 3*jj-1, ii);
    end
    myfilter = fspecial('gaussian', [1 1], 1);
    smoothImage = imfilter(histEstimate(:,:,ii), myfilter, 'replicate');
    histEstimate(:,:,ii) = round(smoothImage*255/max(smoothImage(:)));

    figure(hScatter)
    subplot(1,length(stdNoiseLevel),ii)
    hold on
    tempImage = uint8(histEstimate(:,:,ii));
    [height, width] = size(tempImage);
    widthPlot = round(1*width);
    tempImage = imresize(tempImage, [height widthPlot], 'bicubic');
    tempImage = max(tempImage(:)) - tempImage;
    imshow(tempImage)
%     sumImage = sumImage + tempImage;
    axis on
    hold on
    set(gca,'FontSize',fontSize)
    set(gca, 'ylim', [1 height], 'xlim', [1 widthPlot], ...
        'XTick', round(linspace(1,widthPlot,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
        'YTick', round(linspace(1,height,11)), 'YTickLabel', num2cell(round(linspace(maxYplot,-maxYplot,11))))
    plot([1 widthPlot], [round(length(-maxYplot:maxYplot)/2)+22 round(length(-maxYplot:maxYplot)/2)-22], '--b', 'LineWidth', 1.1)
%     plot([width+1 width+1], [1 2*maxYplot], '--b', 'LineWidth', 1.1)    
    if strcmp(experimentType, 'MainExperiment')
        title(['Stimulus noise = ' num2str(stdNoiseLevel(ii))])
    else
        title('Motor noise')
    end
    xlabel('True angle (degree)')
    if ii == 1
        ylabel('Estimated angle (degree)')
    end
end


%% Density plots of estimates
maxEstimate = 40;
minEstimate = -40;
yAxis = minEstimate:maxEstimate;
xAxis = angleDiff;

% Check the histogram plot
figure
for ii = 1 : length(stdNoiseLevel)
    tempAngleDiffEst = angleDiffEstimate{ii, angleDiff == 3};
    tempHist = hist(tempAngleDiffEst, yAxis);  
    subplot(1, 3, ii)
    bar(yAxis, tempHist)
    xlim([-40 40])
end
if collapseScatterPlot
    histEstimate = zeros(length(yAxis), length(xAxis));
else
    histEstimate = cell(1, length(stdNoiseLevel));
end
figure
hold on
for ii = 1 : length(stdNoiseLevel)
    histEstimate_temp = zeros(length(yAxis), length(xAxis));
    for jj = 1 :length(angleDiff)
        tempAngleDiffEst = angleDiffEstimate{ii,jj};
        histEstimate_temp(:, jj) = hist(tempAngleDiffEst, yAxis);        
    end
    histEstimate_temp = max(histEstimate_temp(:)) - histEstimate_temp;
    if collapseScatterPlot
        histEstimate = histEstimate + histEstimate_temp;
    else
        histEstimate{ii} = histEstimate_temp;
        subplot(1, 3, ii)
        imagesc(histEstimate_temp)
        colormap gray
        axis xy
    end
end
if collapseScatterPlot
    imagesc(histEstimate)
    colormap gray
    axis xy
    axis off
end
save('densityEstimateExp1_Incongruent', 'histEstimate', 'xAxis', 'yAxis')

%% Plot the average
hAverage = figure;
lineWidth = 2;
figPos = [400, 100, 700, 700];
set(hAverage,'Units','pixels','Position',figPos)
hold on
legendName = cell(1,length(stdNoiseLevel));
colorName = {'g', 'r', 'b', 'cyan', 'magenta', 'y'};
maxYplot = max(params.barAngleDiff);
maxXplot = max(params.barAngleDiff);
indexColor = 1;
hLegend = NaN(1, length(stdNoiseLevel));
boxcarLength = 1;
if strcmp(experimentType, 'MainExperiment')
    % Plot true angle
    rangeCollapse = round(length(angleDiff)/2);
    estimateCorrectData = cell(length(stdNoiseLevel), rangeCollapse);
    estimateIncorrectData = cell(length(stdNoiseLevel), rangeCollapse);
    biasMeanCorrect = NaN(length(stdNoiseLevel), rangeCollapse);
    biasSEMCorrect = NaN(length(stdNoiseLevel), rangeCollapse);
    biasMeanIncorrect = NaN(length(stdNoiseLevel), rangeCollapse);
    biasSEMIncorrect = NaN(length(stdNoiseLevel), rangeCollapse);
    biasMeanAll = NaN(length(stdNoiseLevel), rangeCollapse);
    biasSEMAll = NaN(length(stdNoiseLevel), rangeCollapse);
    
    stdMeanCorrect = NaN(length(stdNoiseLevel), rangeCollapse);
    stdSEMCorrect = NaN(length(stdNoiseLevel), rangeCollapse);
    stdMeanIncorrect = NaN(length(stdNoiseLevel), rangeCollapse);
    stdSEMIncorrect = NaN(length(stdNoiseLevel), rangeCollapse);
    stdMeanAll = NaN(length(stdNoiseLevel), rangeCollapse);
    stdSEMAll = NaN(length(stdNoiseLevel), rangeCollapse);        
    signOrientation = repmat([-ones(1,rangeCollapse) ones(1,rangeCollapse-1)],length([angleDiffEstimate{1,:}]),1);
    
    for ii = 1 : length(stdNoiseLevel)
        hold on
        set(gca,'FontSize',fontSize)
        tempAngleDiffEstimate = [angleDiffEstimate{ii,:}];
        tempBinaryDecision = [binaryDecision{ii,:}];
        indexCorrect = tempBinaryDecision == signOrientation;
        indexCorrect(:,rangeCollapse) = 1;
        
        % Calculate bias for correct trials
        if ~includeIncorrectTrial
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
        else
            tempAngleDiffEstimateCW_Correct = NaN(size(tempAngleDiffEstimate));
            tempAngleDiffEstimateCW_Correct(tempBinaryDecision == 1) = tempAngleDiffEstimate(tempBinaryDecision == 1);
            tempAngleDiffEstimateCWAve = squeeze(nanmean(tempAngleDiffEstimateCW_Correct,1));
            tempAngleDiffEstimateCWSEM = squeeze(nanstd(tempAngleDiffEstimateCW_Correct,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCW_Correct),1));
            tempBias1 = tempAngleDiffEstimateCW_Correct(:,rangeCollapse:end) - repmat(angleDiff(rangeCollapse:end),size(tempAngleDiffEstimateCW_Correct,1),1);
            tempBias1 = tempBias1 - repmat(biasMotor(rangeCollapse:end), size(tempBias1,1),1);

            tempAngleDiffEstimateCCW_Correct = NaN(size(tempAngleDiffEstimate));
            tempAngleDiffEstimateCCW_Correct(tempBinaryDecision == -1 ) = tempAngleDiffEstimate(tempBinaryDecision == -1);
            tempAngleDiffEstimateCCWAve = squeeze(nanmean(tempAngleDiffEstimateCCW_Correct,1));
            tempAngleDiffEstimateCCWSEM = squeeze(nanstd(tempAngleDiffEstimateCCW_Correct,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCCW_Correct),1));            
        end
            
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
            tempEst = tempEstimateCollapse(:,jj+rangeCollapse-1);
            tempEst(isnan(tempEst)) = [];
            stdMeanCorrect(ii, jj) = std(tempEst);
            bootstat = bootstrp(1000, @std, tempEst);
            stdSEMCorrect(ii, jj) = std(bootstat);
        end
        
        % Plot the estimates
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
        box on
        text(6, -5, '+/-1SEM', 'FontSize', 15)
        
        % Calculate bias for incorrect trials
        tempAngleDiffEstimateCW_Incorrect = NaN(size(tempAngleDiffEstimate));
        tempAngleDiffEstimateCW_Incorrect(tempBinaryDecision == 1 & ~indexCorrect) = tempAngleDiffEstimate(tempBinaryDecision == 1 & ~indexCorrect);
        tempAngleDiffEstimateCWAve = squeeze(nanmean(tempAngleDiffEstimateCW_Incorrect,1));
        tempAngleDiffEstimateCWAve(rangeCollapse:end) = NaN;
        tempAngleDiffEstimateCWSEM = squeeze(nanstd(tempAngleDiffEstimateCW_Incorrect,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCW_Incorrect),1));
        tempBias1 = tempAngleDiffEstimateCW_Incorrect(:,1:rangeCollapse) - repmat(angleDiff(rangeCollapse:end),size(tempAngleDiffEstimateCW_Incorrect,1),1);
        tempBias1 = tempBias1 - repmat(biasMotor(rangeCollapse:end), size(tempBias1,1),1);
        
        tempAngleDiffEstimateCCW_Incorrect = NaN(size(tempAngleDiffEstimate));
        tempAngleDiffEstimateCCW_Incorrect(tempBinaryDecision == -1 & ~indexCorrect) = tempAngleDiffEstimate(tempBinaryDecision == -1 & ~indexCorrect);
        tempAngleDiffEstimateCCWAve = squeeze(nanmean(tempAngleDiffEstimateCCW_Incorrect,1));
        tempAngleDiffEstimateCCWAve(1:rangeCollapse-1) = NaN;
        tempAngleDiffEstimateCCWSEM = squeeze(nanstd(tempAngleDiffEstimateCCW_Incorrect,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCCW_Incorrect),1));
        tempBias2 = tempAngleDiffEstimateCCW_Incorrect(:,rangeCollapse:end) - repmat(angleDiff(1:rangeCollapse),size(tempAngleDiffEstimateCCW_Incorrect,1),1);
        tempBias2 = tempBias2 - repmat(biasMotor(1:rangeCollapse), size(tempBias2,1),1);
        tempBias = [tempBias1;-tempBias2(:,sort(1:size(tempBias2,2),'descend'))];
        tempBiasMean = nanmean(tempBias,1);
        bias = [ones(1,boxcarLength-1)*tempBiasMean(1) tempBiasMean ones(1,boxcarLength-1)*tempBiasMean(end)];
        smoothBias = conv(bias, ones(1,boxcarLength)/boxcarLength, 'full');
        biasMeanIncorrect(ii,:) = smoothBias(boxcarLength:boxcarLength+rangeCollapse-1);
        biasSEMIncorrect(ii,:) = squeeze(nanstd(tempBias,1))./sqrt(nansum(~isnan(tempBias),1));
        
        tempEstimateCollapse = [-tempAngleDiffEstimateCW_Incorrect; fliplr(tempAngleDiffEstimateCCW_Incorrect)];
        for jj = 1 : rangeCollapse
            estimateIncorrectData{ii,jj} = tempEstimateCollapse(:,jj);
            tempEst = tempEstimateCollapse(:,jj);
            tempEst(isnan(tempEst)) = [];
            if length(tempEst) > 4
                stdMeanIncorrect(ii, jj) = std(tempEst);
                bootstat = bootstrp(1000, @std, tempEst);
                stdSEMIncorrect(ii, jj) = std(bootstat);         
            end
        end   
        
        % Calculate bias for all trials
        tempAngleDiffEstimateCCW_Correct(:, rangeCollapse) = -tempAngleDiffEstimateCCW_Correct(:, rangeCollapse);
        tempEstimateCollapse = [-fliplr(tempAngleDiffEstimateCCW_Correct); tempAngleDiffEstimateCW_Correct;...
                                -fliplr(tempAngleDiffEstimateCW_Incorrect); tempAngleDiffEstimateCCW_Incorrect];
        for jj = 1 : rangeCollapse
            estimateCorrectData{ii,jj} = tempEstimateCollapse(:,jj+rangeCollapse-1);
            tempEst = tempEstimateCollapse(:,jj+rangeCollapse-1);
            tempEst(isnan(tempEst)) = [];
            biasMeanAll(ii, jj) = nanmean(tempEst) - angleDiff(jj+rangeCollapse-1);
            biasSEMAll(ii, jj) = nanstd(tempEst) / sqrt(length(tempEst));
            stdMeanAll(ii, jj) = std(tempEst);
            bootstat = bootstrp(1000, @std, tempEst);
            stdSEMAll(ii, jj) = std(bootstat);
        end
        
    end
    legend(hLegend, legendName, 'Location', 'NorthWest')
    
    % Plot the bias 
    hBiasCorrect = figure;
    figPos = [400, -40, 800, 800];
    set(hBiasCorrect,'Units','pixels','Position',figPos)
    hold on
    for ii = 1 : length(stdNoiseLevel)
        % The mean bias
        subplot(1, 2, 1)
        hold on
        set(gca,'FontSize',fontSize)
        hShade = shadedErrorBar(angleDiff(rangeCollapse:end), biasMeanCorrect(ii,:), biasSEMCorrect(ii,:),... 
                             {'Color', colorName{ii}, 'LineWidth', lineWidth},0,1,0);        
        hLegend(ii) = hShade.mainLine;
        legendName{ii} = ['noise = ' num2str(stdNoiseLevel(ii))];
        xlabel('True angle (degree)')
        ylabel('Bias (degree)') 
%         title(['Subject ' upper(subjectID) ]);
        xlim([0 maxXplot])
        ylim([-1 12])
        set(gca, 'YTick', [0:2:12])            
        plot([0 maxXplot], [0 0], 'k--', 'LineWidth', lineWidth)
        
        % The mean std
        subplot(1, 2, 2)
        hold on
        set(gca,'FontSize',fontSize)
        hShade = shadedErrorBar(angleDiff(rangeCollapse:end), stdMeanCorrect(ii,:), stdSEMCorrect(ii,:),... 
                             {'Color', colorName{ii}, 'LineWidth', lineWidth},0,0,0);        
        hLegend(ii) = hShade.mainLine;
        legendName{ii} = ['noise = ' num2str(stdNoiseLevel(ii))];
        xlabel('True angle (degree)')
        ylabel('Standard deviation (degree)') 
        xlim([0 maxXplot])
        ylim([0 8])        
    end
%     legend(hLegend, legendName, 'Location', 'NorthEast')

%     hBiasAll = figure;
%     figPos = [400, -40, 800, 800];
%     set(hBiasAll,'Units','pixels','Position',figPos)
%     hold on
%     for ii = 1 : length(stdNoiseLevel)
%         % The mean bias
%         subplot(1, 2, 1)
%         hold on
%         set(gca,'FontSize',fontSize)
%         hShade = shadedErrorBar(angleDiff(rangeCollapse:end), biasMeanAll(ii,:), biasSEMAll(ii,:),... 
%                              {'Color', colorName{ii}, 'LineWidth', lineWidth},0,0,0);        
%         hLegend(ii) = hShade.mainLine;
%         legendName{ii} = ['noise = ' num2str(stdNoiseLevel(ii))];
%         xlabel('True angle (degree)')
%         ylabel('Bias (degree)') 
% %         title(['Subject ' upper(subjectID) ]);
%         xlim([0 maxXplot])
%         ylim([-3 4])
%         plot([0 maxXplot], [0 0], 'k--', 'LineWidth', lineWidth)
%         
%         % The mean std
%         subplot(1, 2, 2)
%         hold on
%         set(gca,'FontSize',fontSize)
%         hShade = shadedErrorBar(angleDiff(rangeCollapse:end), stdMeanAll(ii,:), stdSEMAll(ii,:),... 
%                              {'Color', colorName{ii}, 'LineWidth', lineWidth},0,0,0);        
%         hLegend(ii) = hShade.mainLine;
%         legendName{ii} = ['noise = ' num2str(stdNoiseLevel(ii))];
%         xlabel('True angle (degree)')
%         ylabel('Standard deviation (degree)') 
%         xlim([0 maxXplot])
%         ylim([0 15])        
%     end
%     legend(hLegend, legendName, 'Location', 'NorthEast')
    
%     hBiasError = figure;
%     figPos = [400, -40, 800, 800];
%     set(hBiasError,'Units','pixels','Position',figPos)
%     hold on
%     for ii = 1 : length(stdNoiseLevel)
%         hold on
%         set(gca,'FontSize',fontSize)
%         hShade = shadedErrorBar(angleDiff(rangeCollapse:end), abs(biasMeanError(ii,:)), biasSEMError(ii,:),... 
%                              {'Color', colorName{ii}, 'LineWidth', lineWidth});        
%         hLegend(ii) = hShade.mainLine;
%         legendName{ii} = ['noise = ' num2str(stdNoiseLevel(ii))];
%         xlabel('Absolute true angle (degree)')
%         ylabel('Bias (degree)') 
% %         title(['Subject ' upper(subjectID) ]);
%         xlim([0 maxXplot])
%         if strcmp(experiment, 'ControlReplication')
%             ylim([-4 20])
%             set(gca, 'XTick', [-4:4:20])
%         else
%             ylim([0 40])
%             set(gca, 'YTick', [0:5:40])            
%         end
%         plot([0 maxXplot], [0 0], 'k--', 'LineWidth', lineWidth)
%     end
%     legend(hLegend, legendName, 'Location', 'NorthEast')
    if strcmp(experiment, 'TrainArray')
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
        xlim([-maxXplot maxXplot])
        ylim([-maxYplot maxYplot])
        xlabel('\theta_{true} (degree)')
        ylabel('\theta_{estimate} (degree)')  
        title('Motor noise');
        plot([angleDiff(1) angleDiff(end)],  [angleDiff(1) angleDiff(end)], '--r', 'LineWidth', 2)
    end
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
    biasMotor = meanMotor - angleDiff;
    xlim([-maxXplot maxXplot])
    ylim([-maxYplot maxYplot])
    xlabel('\theta_{true} (degree)')
    ylabel('\theta_{estimate} (degree)')  
    title('Motor noise');
    plot([angleDiff(1) angleDiff(end)],  [angleDiff(1) angleDiff(end)], '--r', 'LineWidth', 2)
    
    % Collapse and sort the bias in trialOrder
    lengthBias = max(indexTrialOrder);
    biasSorted = NaN(1, lengthBias);
    for ii = 1 : lengthBias
        biasSorted(ii) = nanmean(biasAll(trialOrder(indexTrialOrder==ii)));
    end
    biasSorted(isnan(biasSorted)) = [];
    
    % Compute the std for all, first and second half of trials
    nDataPoints = length(biasSorted);
    stdMotor = nanstd(biasAll);
    indexFirstHalf = trialOrder(indexTrialOrder < round(nDataPoints/2));
    stdMotorFirstHalf = nanstd(biasAll(indexFirstHalf));
    indexSecondHalf = trialOrder(indexTrialOrder >= round(nDataPoints/2));
    stdMotorSecondHalf = nanstd(biasAll(indexSecondHalf));
    display(['Std for all, first and second half: ' num2str([stdMotor stdMotorFirstHalf stdMotorSecondHalf])])
    
    % Compute the error across trials    
    error = abs(biasSorted);
    errorSmooth = smooth(error, 1);
    figure
    plot(errorSmooth, 'o-')
    xlabel('Trial bin')
    ylabel('Error')
elseif strcmp(experiment, 'ControlImplicitDecision')
    stdMotor = NaN(length(stdNoiseLevel), length(angleDiff));
    meanMotor = NaN(length(stdNoiseLevel), length(angleDiff));
    hold on
    for ii = 1 : length(stdNoiseLevel)
        subplot(1,3,ii)
        set(gca,'FontSize',fontSize)
        hold on
        for jj = 1 : length(angleDiff)
            tempEst = [angleDiffEstimate{ii,jj}];
            meanMotor(ii,jj) = nanmean(tempEst);
            stdMotor(ii,jj) = nanstd(tempEst);
            errorbar(angleDiff(jj), ...
                    meanMotor(ii,jj), 2*stdMotor(ii,jj)/sqrt(sum(~isnan(tempEst))),...
                    'bo', 'MarkerFaceColor','b', 'MarkerSize', 9)
        end
        xlim([-maxXplot maxXplot])
        ylim([-maxYplot maxYplot])
        xlabel('True orientation (degree)')
        ylabel('Orientation estimate (degree)')  
        title(['Stimulus noise = ' num2str(stdNoiseLevel(ii))])
        plot([-maxXplot maxXplot],  [-maxYplot maxYplot], '--r', 'LineWidth', 2)
    end    
    tightfig
end
% optBinsCorrect = [];
% save(estimateDataName, 'estimateCorrectData', 'estimateIncorrectData', 'angleDiff', 'optBinsCorrect')

%% Calculate percent clockwise
if ~strcmp(experimentType, 'MotorNoise')  
    percentCW = NaN(length(stdNoiseLevel),length(params.barAngleDiff));
    nTrialsPerCondition = NaN(length(stdNoiseLevel),length(params.barAngleDiff));
    if strcmp(experiment, 'DecisionGiven')
        binaryDecisionAll = dataAll(:, 8);
    else
        binaryDecisionAll = dataAll(:, 6);
    end    
    for ii = 1 : length(stdNoiseLevel)
        for jj = 1 : length(percentCW)
            indexAll = dataAll(:,1) == params.barAngleDiff(jj) ...
                        & dataAll(:,4) == stdNoiseLevel(ii)...
                        & indicatorConsistent;
            indexCW = dataAll(:,1) == params.barAngleDiff(jj) ...
                        & dataAll(:,4) == stdNoiseLevel(ii)...
                        & indicatorConsistent ...
                        & binaryDecisionAll == 1;
            nTrialsPerCondition(ii,jj) = nansum((indexAll));
            percentCW(ii,jj) = 100*nansum(indexCW)/nansum(indexAll);
        end
    end
end
if extractModelParam || strcmp(experimentType, 'PerceptNoise')
    if isfield(params, 'staircaseTrack')
        nTrialPerBin = 20;
        % Bin the staircase values
        meanValues = cell(1,length(stdNoiseLevel));
        nCW = cell(1,length(stdNoiseLevel));
        percentCW = cell(1,length(stdNoiseLevel));
        for ii = 1 : length(stdNoiseLevel)
            columnIndex = dataAll(:,4)== stdNoiseLevel(ii);
            values = dataAll(columnIndex,1);
            response = dataAll(columnIndex,6);
            [meanValues{ii},nCW{ii}, nTrials] = GetAggregatedStairTrials(values, response, nTrialPerBin);
            percentCW{ii} = 100*nCW{ii}/nTrialPerBin;
        end
    end
    
    % Plot the fit curve and data
    lineWidth = 3;
    fontSize = 20;  
    figure
    hold on
    set(gca, 'FontSize', fontSize)
    colorName = {'g','r', 'b', 'cyan', 'magenta', 'y'};
    hLegend = NaN(1, length(stdNoiseLevelAll));    
    angleDiffResampled = linspace(angleDiff(1),angleDiff(end),1000);
    fitPercentCW = 100*BayesFractionCW(priorRangeAll, smoothFactorAll, stdNoiseLevelAll, lapseRate, angleDiffResampled, windowPrior);
    
    % Save the psychometric curve
    psychCurveExp2.angleDiff = angleDiffResampled;
    psychCurveExp2.percentCW = fitPercentCW;

    legendName = cell(1,length(stdNoiseLevelAll));
    for ii = 1 : length(stdNoiseLevelAll)
        hLegend(ii) = plot(angleDiffResampled,fitPercentCW(ii,:),'Color',...
            colorName{ii},'LineWidth',3);
        plot(angleDiff, percentCW(ii,:), [colorName{ii} 'o'], 'MarkerSize',13,...
            'MarkerFaceColor', colorName{ii}, 'LineWidth', lineWidth);
        legendName{ii} = ['noise = ' num2str(roundn(stdNoiseLevelAll(ii),-1))];
    end
    xlim([min(angleDiff) max(angleDiff)])
    set(gca,'YTick', 0:20:100)
end
% if ~strcmp(experiment, 'DecisionGiven')
%     save(discriminationDataName, 'binaryDecision', 'percentCW', 'nTrialsPerCondition')
% end

%% Plot the smoothed raw data
if collapseScatterPlot
    hScatter = figure;
    figPos = [0.1, 0.2, 0.8, 0.6];
    set(hScatter,'Units','normalized','Position',figPos)
    hold on
    maxYplot = 35;
    histEstimate = zeros(length(-maxYplot:maxYplot), length(-22:22));
    for ii = 1 : length(stdNoiseLevel)
        for jj = 1 :length(angleDiff)
            tempAngleDiffEst = angleDiffEstimate{ii,jj};
            tempAngleDiffEst(abs(tempAngleDiffEst)>35) = [];
            
            for kk = 1 : length(tempAngleDiffEst)
                if ~isnan(tempAngleDiffEst(kk)) && (round(tempAngleDiffEst(kk))+1<maxYplot)
                histEstimate(maxYplot-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23) = ...
                    histEstimate(maxYplot-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23) + 1;
                end
            end
        end
    end
    for jj = 1 :length(angleDiff)
        histEstimate(:, 3*jj-2) = histEstimate(:, 3*jj-1);
        histEstimate(:, 3*jj) = histEstimate(:, 3*jj-1);
    end
    myfilter = fspecial('gaussian', [1 1], 1);
    smoothImage = imfilter(histEstimate, myfilter, 'replicate');
    histEstimate= round(smoothImage*255/max(smoothImage(:)));

    % Plot
    figure
    hold on
    tempImage = uint8(histEstimate);
    [height, width] = size(tempImage);
    widthPlot = round(1*width);
    tempImage = imresize(tempImage, [height widthPlot], 'bicubic');
    imshow(tempImage)
    axis on
    hold on
    set(gca,'FontSize',fontSize)
    set(gca, 'ylim', [1 height], 'xlim', [1 widthPlot], ...
        'XTick', round(linspace(1,widthPlot,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
        'YTick', round(linspace(1,height,5)), 'YTickLabel', num2cell(round(linspace(maxYplot,-maxYplot,5))))
    plot([1 widthPlot], [round(length(-maxYplot:maxYplot)/2)+22 round(length(-maxYplot:maxYplot)/2)-22], '--b', 'LineWidth', 1.1)
    xlabel('True angle (degree)')
    ylabel('Estimated angle (degree)')
else
    hScatter = figure;
    figPos = [0.1, 0.2, 0.8, 0.6];
    set(hScatter,'Units','normalized','Position',figPos)
    hold on
    maxYplot = 35;
    histEstimate = zeros(length(-maxYplot:maxYplot), length(-22:22), length(stdNoiseLevel));
    for ii = 1 : length(stdNoiseLevel)
        for jj = 1 :length(angleDiff)
            tempAngleDiffEst = angleDiffEstimate{ii,jj};
            tempAngleDiffEst(abs(tempAngleDiffEst)>35) = [];
            for kk = 1 : length(tempAngleDiffEst)
                if ~isnan(tempAngleDiffEst(kk)) && (round(tempAngleDiffEst(kk))+1<maxYplot)
                histEstimate(maxYplot-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) = ...
                    histEstimate(maxYplot-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) + 1;
                end
            end
        end
        for jj = 1 :length(angleDiff)
            histEstimate(:, 3*jj-2, ii) = histEstimate(:, 3*jj-1, ii);
            histEstimate(:, 3*jj, ii) = histEstimate(:, 3*jj-1, ii);
        end
        myfilter = fspecial('gaussian', [1 1], 1);
        smoothImage = imfilter(histEstimate(:,:,ii), myfilter, 'replicate');
        histEstimate(:,:,ii) = round(smoothImage*255/max(smoothImage(:)));

        subplot(1,length(stdNoiseLevel),ii)
        hold on
        tempImage = uint8(histEstimate(:,:,ii));
        [height, width] = size(tempImage);
        widthPlot = round(1*width);
        tempImage = imresize(tempImage, [height widthPlot], 'bicubic');
        imshow(tempImage)
    %     sumImage = sumImage + tempImage;
        axis on
        hold on
        set(gca,'FontSize',fontSize)
        set(gca, 'ylim', [1 height], 'xlim', [1 widthPlot], ...
            'XTick', round(linspace(1,widthPlot,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
            'YTick', round(linspace(1,height,11)), 'YTickLabel', num2cell(round(linspace(maxYplot,-maxYplot,11))))
        plot([1 widthPlot], [round(length(-maxYplot:maxYplot)/2)+22 round(length(-maxYplot:maxYplot)/2)-22], '--b', 'LineWidth', 1.1)
    %     plot([width+1 width+1], [1 2*maxYplot], '--b', 'LineWidth', 1.1)    
        if strcmp(experimentType, 'MainExperiment')
            title(['Stimulus noise = ' num2str(stdNoiseLevel(ii))])
        else
            title('Motor noise')
        end
        xlabel('True angle (degree)')
        if ii == 1
            ylabel('Estimated angle (degree)')
        end
    end
end

%% Plot the response-divided version of smoothed raw data (red for CW and green for CCW)
% hScatterColor = figure;
% figPos = [0.1, 0.2, 0.8, 0.6];
% set(hScatterColor,'Units','normalized','Position',figPos)
% hold on
% maxYplotColor = 60;
% imageDataCW = zeros(length(-maxYplotColor:maxYplotColor), length(-22:22), length(stdNoiseLevel));
% imageDataCCW = zeros(length(-maxYplotColor:maxYplotColor), length(-22:22), length(stdNoiseLevel));
% imageDataBlue = zeros(length(-maxYplotColor:maxYplotColor), length(-22:22));
% for ii = 1 : length(stdNoiseLevel)
%     for jj = 1 :length(angleDiff)
%         tempAngleDiffEst = angleDiffEstimate{ii,jj};
%         tempBinaryDecision = binaryDecision{ii,jj};
%         for kk = 1 : length(tempAngleDiffEst)
%             if ~isnan(tempAngleDiffEst(kk)) && (round(tempAngleDiffEst(kk))+1<maxYplot)
%             if tempBinaryDecision(kk) == -1;
%                 imageDataCCW(maxYplotColor-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) = ...
%                     imageDataCCW(maxYplotColor-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) + 1;
%             elseif tempBinaryDecision(kk) == 1;
%                 imageDataCW(maxYplotColor-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) = ...
%                     imageDataCW(maxYplotColor-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) + 1;
%             end
%             end
%         end
%     end
%     imageData = cat(3, imageDataCW(:,:,ii), imageDataCCW(:,:,ii), imageDataBlue);
%     myfilter = fspecial('gaussian', [3 2], 1);
%     smoothImage = imfilter(imageData, myfilter, 'replicate');
%     smoothImage(:,:,1) = smoothImage(:,:,1)*255/max(max(smoothImage(:,:,1)));
%     smoothImage(:,:,2) = smoothImage(:,:,2)*255/max(max(smoothImage(:,:,2)));    
%     imageData = round(smoothImage);
%     
%     subplot_tight(1,length(stdNoiseLevel),ii,[0.03 0.03])
%     hold on
%     tempImage = uint8(imageData);
%     [height, width, channel] = size(tempImage);
%     widthPlot = 2*width;
%     tempImage = imresize(tempImage, [height widthPlot], 'bicubic');
%     imshow(tempImage)
%     [height,~,~] = size(tempImage);
%     axis on
%     hold on
%     set(gca,'FontSize',fontSize)
%     set(gca, 'ylim', [1 height], 'xlim', [1 widthPlot], ...
%         'XTick', round(linspace(3,widthPlot-2,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
%         'YTick', round(linspace(1,height,7)), 'YTickLabel', num2cell(round(linspace(maxYplotColor,-maxYplotColor,7))))
%     if ii==2 || ii==3
%         set(gca, 'YTickLabel',{''})
%     end
%     plot([1 widthPlot], [round(length(-maxYplotColor:maxYplotColor)/2)+22 round(length(-maxYplotColor:maxYplotColor)/2)-22], '--b', 'LineWidth', lineWidth)
% %     plot([width+1 width+1], [1 2*maxYplot], '--b', 'LineWidth', 1.1)  
% %     xlabel('True angle (degree)')
%     if ii == 1
%         ylabel('Estimated angle (degree)', 'FontSize', 23)
%     end
%     title(['Stimulus noise = ' num2str(stdNoiseLevel(ii))], 'FontSize', 23)
% end
% tightfig


%% Check the color discrimination task and consistence index
angleDiffEstimateError = cell(length(stdNoiseLevel), length(angleDiff));
if strcmp(experiment, 'DecisionGiven')
    load psychCurveExp2
    % Color task performance
    colorDetection = dataAll(:, 6);
    colorErrorPercent = NaN(1, length(stdNoiseLevel));
    for ii = 1 : length(stdNoiseLevel)
        indexColorError = dataAll(:,4)==stdNoiseLevel(ii) & (~colorDetection);
        colorErrorPercent(ii) = 100 * nansum(indexColorError) / nansum(dataAll(:,4)==stdNoiseLevel(ii));
    end
    percentCorrectColor = 100 * nansum(colorDetection) / sum(~isnan(colorDetection));
    disp(['Percent error color detection: ' num2str(100-percentCorrectColor)])
    
    % Deviation of estimate from decision
    indexInconsistenceTotal = ~(sign(angleDiffEst) == sign(dataAll(:,8)));
    numInconsistentTotal = nansum(indexInconsistenceTotal);
    numErrorColor = sum(~dataAll(:, 6));
    numErrorColorAndInconsistent = nansum(~dataAll(:, 6) & indexInconsistenceTotal);
    percentErrorTotal = 100*numInconsistentTotal / length(angleDiffEst);
    numComply = NaN(length(stdNoiseLevel), length(angleDiff));
    for ii = 1 : length(stdNoiseLevel)
        for jj = 1 : length(angleDiff)
            indexCondition = dataAll(:,4)==stdNoiseLevel(ii) & dataAll(:,1)==angleDiff(jj);
            indexCorrect = sign(angleDiffEst(indexCondition)) == sign(dataAll(indexCondition,8));
            numComply(ii,jj) = nansum(indexCorrect); 
            tempAngleDiffEst = angleDiffEst(indexCondition);
            angleDiffEstimateError{ii,jj} = tempAngleDiffEst(~indexCorrect);
        end
    end
    percentError = 100 - 100*numComply / sum(indexCondition);
    disp(['Percent inconsistence error: ' num2str(percentErrorTotal)])
    disp(['Percent inconsistence error given color detection error: ' num2str(100*numErrorColorAndInconsistent/numErrorColor)])
%     [corr, pvalue] = corr(indexInconsistenceTotal, indexErrorColor)
    
else
    % Deviation of estimate from decision
    indexInconsistenceTotal = ~indicatorConsistentKeep;
    numInconsistentTotal = nansum(indexInconsistenceTotal);   
    percentErrorTotal = 100*numInconsistentTotal / length(angleDiffEst);
    numComply = NaN(length(stdNoiseLevel), length(angleDiff));
    for ii = 1 : length(stdNoiseLevel)
        for jj = 1 : length(angleDiff)
            indexCondition = dataAll(:,4)==stdNoiseLevel(ii) & dataAll(:,1)==angleDiff(jj);
            indexCorrect = sign(angleDiffEst(indexCondition)) == sign(dataAll(indexCondition,6));
            numComply(ii,jj) = nansum(indexCorrect); 
            tempAngleDiffEst = angleDiffEst(indexCondition);
            angleDiffEstimateError{ii,jj} = tempAngleDiffEst(~indexCorrect);
        end
    end
    percentError = 100 - 100*numComply / sum(indexCondition);
    disp(['Percent inconsistence error: ' num2str(percentErrorTotal)])
end

% Get the hypothesized error if subjects did implicit decision
angleDiffResampled = psychCurveExp2.angleDiff;
errorTheoryImplicit = psychCurveExp2.percentCW;
errorTheoryImplicit(:,angleDiffResampled>0) = fliplr(errorTheoryImplicit(:,angleDiffResampled<0));

% Plot the inconsistence error
colorName = {'g', 'r', 'b', 'cyan', 'magenta', 'y'};
hLegend = NaN(1, length(stdNoiseLevel));
hError = figure;
figPos = [400, -40, 800, 800];
set(hError,'Units','pixels','Position',figPos)
lineWidth = 4;
legendName = {'Implicit', 'Data', 'No implicit'};
areaError = hLegend;
areaErrorTheoryImplicit = hLegend;
lapseRateColor = 1.5;
for ii = 1 : length(stdNoiseLevel)
    subplot_tight(1,3,ii,[0.1 0.04])
    hold on
    set(gca, 'FontSize', fontSize)
    hLegend(2) = plot(angleDiff, percentError(ii,:), 'b-o', 'LineWidth', lineWidth,...
        'MarkerSize', 13, 'MarkerFaceColor', 'b');
    if strcmp(experiment, 'DecisionGiven')
        %Plot hypothesized error
        hLegend(1) = plot(angleDiffResampled, errorTheoryImplicit(ii,:), 'r-.', 'LineWidth', lineWidth);
        hLegend(3) = plot([-21 21], [lapseRateColor lapseRateColor], 'Color', [0.5 0.5 0.5],...
            'LineStyle', '-.', 'LineWidth', lineWidth);
    end  
    xlim([-22 22])
    ylim([0 50])
    xlabel('True orientation')
    if ii == 1
        ylabel('Percent error (%)')
    end
    title(['noise = ' num2str(stdNoiseLevel(ii))])  
%     legend(hLegend, legendName, 'Location', 'NorthEast')
    box on
    areaError(ii) = mean(percentError(ii,:));
    areaErrorTheoryImplicit(ii) = mean(errorTheoryImplicit(ii,:));   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fit to all data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Fixed motor noise, varied lapse
% modelInconsistentError =      [3.6 4.6 1];
% modelInconsistentErrorLapse = [2.6 4.3 0];
% modelInconsistentErrorStack = modelInconsistentError - modelInconsistentErrorLapse;

% % Varied noise, varied lapse (forget stimulus)
% modelInconsistentError =      [3.9 5.9 3.7];
% modelInconsistentErrorLapse = [2.54 3.9  0];
% modelInconsistentErrorStack = modelInconsistentError - modelInconsistentErrorLapse;

% % Varied noise, varied lapse (remember stimulus)
% modelInconsistentError =      [3.3 5.6 2.6];
% modelInconsistentErrorLapse = [2.6 4.1  0];
% modelInconsistentErrorStack = modelInconsistentError - modelInconsistentErrorLapse;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fit to congruent data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Fit congruent data to get motor noise, fit all discrimination data to get
% % lapse rate, then predict incongruent error
% modelInconsistentError =      [3.5 5.4 2];
% modelInconsistentErrorLapse = [2.3 4.9 1.7];
% modelInconsistentErrorStack = modelInconsistentError - modelInconsistentErrorLapse;

% % Fit congruent data for individual subject with fixed motor noise, get the average parameter,
% % then predict incongruent error
% modelInconsistentError =      [4.9 4.9 1.7];
% modelInconsistentErrorLapse = [2.14 4.5 1.4];
% modelInconsistentErrorStack = modelInconsistentError - modelInconsistentErrorLapse;

% % Fit congruent data for average subject with fixed motor noise,
% % then predict incongruent error
% modelInconsistentError =      [3.9 5.2 1.53];
% modelInconsistentErrorLapse = [2.3 4.9 1.23];
% modelInconsistentErrorStack = modelInconsistentError - modelInconsistentErrorLapse;

% hBarAve = figure;
% load areaErrorExp2
% figPos = [400, -40, 800, 800];
% set(hBarAve,'Units','pixels','Position',figPos)
% hold on
% barValue = [mean(areaErrorTheoryImplicit) mean(areaError) modelInconsistentErrorLapse(3) mean(areaErrorExp2)  modelInconsistentErrorLapse(2)  3.0222       modelInconsistentErrorLapse(1)      mean(colorErrorPercent);
%                         0                           0     modelInconsistentErrorStack(3)          0        modelInconsistentErrorStack(2)  0            modelInconsistentErrorStack(1)        0]';
% h = bar(barValue, 'stack');
% set(gca,'XTickLabel',{'Implicit',          'Data (Exp3)',    'Model (Exp3)',             'Data (Exp2)',         'Model (Exp2)',         'Data (Exp1)',      'Model (Exp1)',                  'Color Error'},...
%          'XTick', 1:size(barValue, 1), 'FontSize', fontSize) 
% xlabel('Noise level')
% ylabel('Percent error (%)')
% ylim([0 1.1*mean(areaErrorTheoryImplicit)])

% Fit congruent data for individual subject with fixed motor noise, then
% predict incongruent error for each subject and average the errors
areaErrorTheoryImplicit = [8.85  15.86  11.22   10.04   11.34];
modelInconsistentError =      [1.44     1.6     2.2;
                               1        8.1     2.1;
                               5.6      14.1    4.14;
                               16.7     1.6     2.73;
                               4.1      3.2     0.5];
                           
modelInconsistentErrorLapse = [0        0       2.2;
                               0        7.9     1;
                               5.5      12      4.14;
                               1.8      0.5     2.5;
                               2.9      3.2     0];
                           
errorIncongruentData =         [0.11    0       0;
                                0.78    3.28  3.83;
                                5.78    14.44  6.56;
                                4.61    2.83  5.00;
                                3.83    4.33  2.61];
modelInconsistentError =      mean(modelInconsistentError, 1);
modelInconsistentErrorLapse = mean(modelInconsistentErrorLapse, 1);
modelInconsistentErrorStack = modelInconsistentError - modelInconsistentErrorLapse;
                            
errorColor = [0.06  5.67    0.22    1.89  2.1];

hBarAve = figure;
load areaErrorExp2
figPos = [400, -40, 800, 800];
set(hBarAve,'Units','pixels','Position',figPos)
hold on
barValue = [mean(areaErrorTheoryImplicit) mean(errorIncongruentData(:,3)) modelInconsistentErrorLapse(3) mean(errorIncongruentData(:,2))  modelInconsistentErrorLapse(2)  mean(errorIncongruentData(:,1))       modelInconsistentErrorLapse(1)      mean(errorColor);
                        0                           0     modelInconsistentErrorStack(3)          0        modelInconsistentErrorStack(2)  0            modelInconsistentErrorStack(1)        0]';
h = bar(barValue, 'stack');
set(gca,'XTickLabel',{'Implicit',          'Data (Exp3)',    'Model (Exp3)',             'Data (Exp2)',         'Model (Exp2)',         'Data (Exp1)',      'Model (Exp1)',                  'Color Error'},...
         'XTick', 1:size(barValue, 1), 'FontSize', fontSize) 
xlabel('Noise level')
ylabel('Percent error (%)')
ylim([0 1.1*mean(areaErrorTheoryImplicit)])

% hBar = figure;
% load areaError2
% figPos = [400, -40, 800, 800];
% set(hBar,'Units','pixels','Position',figPos)
% hold on
% set(gca, 'FontSize', fontSize)
% h = bar([1:3], [areaErrorTheoryImplicit' areaError2' areaError' colorErrorPercent']);
% set(h(1), 'FaceColor', 'red', 'EdgeColor', 'red')
% set(h(2), 'FaceColor', 'blue', 'EdgeColor', 'blue')
% set(h(3), 'FaceColor', [0 0 0.7], 'EdgeColor', 'blue')
% set(h(4), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5])
% set(gca,'XTickLabel',{'0','6','13'}, 'XTick', [1:3]) 
% xlabel('Noise level')
% ylabel('Percent error (%)')
% ylim([0 1.1*max(areaErrorTheoryImplicit)])
% legend('Implicit', 'Data (Exp 2)', 'Data (Exp3)', 'Color error', 'Location', 'NorthWest')


% Plot the smoothed raw data of error trials
hScatter = figure;
figPos = [0.1, 0.2, 0.8, 0.6];
set(hScatter,'Units','normalized','Position',figPos)
hold on
maxYplot = 35;
histEstimate = zeros(length(-maxYplot:maxYplot), length(-22:22), length(stdNoiseLevel));
estimateErrorCW = NaN(size(angleDiffEstimateError));
estimateErrorCCW = NaN(size(angleDiffEstimateError));

for ii = 1 : length(stdNoiseLevel)
    for jj = 1 :length(angleDiff)
        tempAngleDiffEst = angleDiffEstimateError{ii,jj};
        if angleDiff(jj) > 0
            estimateErrorCW(ii,jj) = mean(tempAngleDiffEst);
        elseif angleDiff(jj) < 0
            estimateErrorCCW(ii,jj) = mean(tempAngleDiffEst);
        else
            estimateErrorCW(ii,jj) = mean(tempAngleDiffEst(tempAngleDiffEst<0));
            estimateErrorCCW(ii,jj) = mean(tempAngleDiffEst(tempAngleDiffEst>0));
        end
        for kk = 1 : length(tempAngleDiffEst)
            if ~isnan(tempAngleDiffEst(kk)) &&  (maxYplot-round(tempAngleDiffEst(kk))+1>0)
            histEstimate(maxYplot-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) = ...
                histEstimate(maxYplot-round(tempAngleDiffEst(kk))+1, round(angleDiff(jj))+23, ii) + 1;
            end
        end
    end
    myfilter = fspecial('gaussian', [3 2], 1);
    smoothImage = imfilter(histEstimate(:,:,ii), myfilter, 'replicate');
    histEstimate(:,:,ii) = round(smoothImage*255/max(smoothImage(:)));
    
    subplot(1,length(stdNoiseLevel),ii)
    hold on
    tempImage = uint8(histEstimate(:,:,ii));
    [height, width] = size(tempImage);
    widthPlot = round(width);
    tempImage = imresize(tempImage, [height widthPlot], 'bicubic');
    imshow(tempImage)
    axis on
    hold on
    set(gca,'FontSize',fontSize)
    set(gca, 'ylim', [1 height], 'xlim', [1 widthPlot], ...
        'XTick', round(linspace(1,widthPlot,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
        'YTick', round(linspace(1,height,5)), 'YTickLabel', num2cell(round(linspace(maxYplot,-maxYplot,5))))
    plot([1 widthPlot], [round(length(-maxYplot:maxYplot)/2)+22 round(length(-maxYplot:maxYplot)/2)-22], '--b', 'LineWidth', lineWidth)
    plot([width+1 width+1], [1 2*maxYplot], '--b', 'LineWidth', 1.1)    
    if strcmp(experimentType, 'MainExperiment')
        title(['Stimulus noise = ' num2str(stdNoiseLevel(ii))])
    else
        title('Motor noise')
    end
    xlabel('True angle (degree)')
    ylabel('Estimated angle (degree)')
end

hAverage = figure;
hold on
rangeCollapse = round(length(angleDiff)/2);

for ii = 1 : length(stdNoiseLevel)
    plot(angleDiff, estimateErrorCCW(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth)
    plot(angleDiff, estimateErrorCW(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth) 
    xlabel('True angle (degree)')
    ylabel('Angle estimate (degree)')
end

% Compute the optimal number of bins
dimEstimateData = size(estimateCorrectData);
optBinsCorrect = NaN(dimEstimateData);
for ii = 1 : dimEstimateData(1)
    for jj = 1 : dimEstimateData(2)
        optBinsCorrect(ii,jj) = sshist(estimateCorrectData{ii,jj});
    end
end

% % Correlation btw the color and inconsistence error
% colorErrorPercent = [0.4 0.89 1.2 5.9]';
% inconsistenceError = [0.06 0.33 0.61 8.6]';
% [r, pValue] = corr(colorErrorPercent, inconsistenceError)
% p = polyfit(colorErrorPercent,inconsistenceError,1);
% fitInconsistent = polyval(p,colorErrorPercent);
% hError = figure;
% hold on
% set(gca, 'FontSize', fontSize)
% plot(colorErrorPercent, inconsistenceError, 'or', 'MarkerSize',13, 'MarkerFaceColor', 'r')
% plot(colorErrorPercent, fitInconsistent, 'LineWidth', lineWidth)
% xlabel('Color detection error (%)')
% ylabel('Inconsistence error (%)')
% box on
% export_fig(hError,'corrError.jpeg','-transparent','-m2')

%% Plot the model's parameters
% parameters: stdNoiseLevel(1:3)  LapseRate    Prior     MemNoise   Smoothness 

%%%%%%%%% Exp 1 %%%%%%%%%
% paramModel =     [2.5611    4.8570    7.1153      0.0000    28.6016      0.8805    0.8461;  %LL
%                   4.0300    5.3505    8.1590      0.0000    41.7649   2.0957    0.0257;  %SY
%                   3.8692    6.4042   12.3882      0.0000    39.5010  14.0727    0.2322;  %CZ
%                   3.7041    4.6566   10.9916      0.0000    36.4739   0.2433    0.6481;  %VS
%                   3.9226    5.7843   10.9330      0.0000    39.0534  14.2463    0.1492  %AS
%                   3.4994    5.4317    9.5625      0.0000    39.7741   6.0913    0.2953]; %Ave

%%%%%%%%% Exp 2, 3 %%%%%%%%%              
paramModel =     [2.5611    4.8570    7.1153      0.0000    23.1044    0.8805    0.8461;  %LL
                  6.3630    8.4075   14.5464      0.0000    22.5672   5.9826    0.7652;  %XFL
                  4.6736    6.1762    7.9466      0.0000    15.8459   1.2144    0.3681;  %AJ
                  4.6585    5.5272    6.8023      0.0000    17.4750   4.1846    0.0934;  %ZW
                  4.3172    4.5344    8.4083      0.0000    33.9063  14.5556    0.5305 ;  %SKB
                  4.6025    6.1782    8.9676      0.0000    23.1643   5.8385    0.8455]; %Ave

noiseSensoryExp1 = paramModel(:, 1:3);
noiseMemoryExp1 = paramModel(:, 6);
noiseAll = [noiseSensoryExp1 noiseMemoryExp1];
priorRange = paramModel(:, 5);

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

%% Plot the bootstrapped parameters
exp = 2;
if exp == 1
    subjID = {'LL';'SY';'CZ';'VS';'AS';'Ave1'};
else
    subjID = {'LL'; 'XFL';'AJ';'ZW';'SKB';'Ave23'};        
end
confidenceIntervalAll = cell(5, length(subjID));
lowerCI = NaN(6, 5);
upperCI = NaN(6, 5);

for kk = 1 : length(subjID);
    fileName = ['FitResult-Bootstrap' subjID{kk} '.txt'];
    fileID = fopen(fileName);
    paramAll = textscan(fileID,'%f %f %f %f %f %f %f %f %f','CommentStyle','//');
    fclose(fileID);
    if strcmp(subjID{kk}, 'LL')
        if exp == 1
            selectInd = [2 3 4 8 6];
        else
            selectInd = [2 3 4 8 7];
        end
    else
        selectInd = [2 3 4 7 6];
    end
    plotFig = 0;
    if plotFig
        h = figure;
        figPos = [0.01, 0.2, 0.98, 0.6];
        set(h,'Units','normalized','Position',figPos)
        hold on
        fontSize = 20;
        set(gca, 'FontSize', fontSize)
    end
    confInterval = NaN(4, 2);
    for ii = 1:5
        % Read data
        tempParam = paramAll{selectInd(ii)};
        confInterval(ii,:) = prctile(tempParam,[2.5 97.5]);
        confidenceIntervalAll{ii, kk} = confInterval(ii,:);
        lowerCI(kk, ii) = confInterval(ii,1);
        upperCI(kk, ii) = confInterval(ii,2);
        
        % Plot histogram
        if plotFig
            subplot(1, 5, ii)
            hold on
            set(gca, 'FontSize', fontSize)
            hist(tempParam, 20)
            title(['Mean = ' num2str(round(mean(tempParam),2)) ', 95% CI = [' num2str(round(confInterval(ii,1),2)) ' - ' num2str(round(confInterval(ii,2),2)) ']'])
        end
    end
end

if exp == 1
    paramModel = [2.5611    4.8570    7.1153      0.0000    28.6016      0.8805    0.8461;  %LL
                  4.0300    5.3505    8.1590      0.0000    41.7649   2.0957    0.0257;  %SY
                  3.8692    6.4042   12.3882      0.0000    39.5010  14.0727    0.2322;  %CZ
                  3.7041    4.6566   10.9916      0.0000    36.4739   0.2433    0.6481;  %VS
                  3.9226    5.7843   10.9330      0.0000    39.0534  14.2463    0.1492  %AS
                  3.4994    5.4317    9.5625      0.0000    39.7741   6.0913    0.2953]; %Ave
else
    paramModel =     [2.5611    4.8570    7.1153      0.0000    23.1044    0.8805    0.8461;  %LL
                      6.3630    8.4075   14.5464      0.0000    22.5672   5.9826    0.7652;  %XFL
                      4.6736    6.1762    7.9466      0.0000    15.8459   1.2144    0.3681;  %AJ
                      4.6585    5.5272    6.8023      0.0000    17.4750   4.1846    0.0934;  %ZW
                      4.3172    4.5344    8.4083      0.0000    33.9063  14.5556    0.5305 ;  %SKB
                      4.6025    6.1782    8.9676      0.0000    23.1643   5.8385    0.8455]; %Ave    
end

noiseSensoryExp1 = paramModel(:, 1:3);
noiseMemoryExp1 = paramModel(:, 6);
noiseAll = [noiseSensoryExp1 noiseMemoryExp1];
priorRange = paramModel(:, 5);

% Plot the parameters with bars
fontSize = 23;
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
hFig = figure;
hAx1 = subplot(1, 2, 1);
[~, hPanel] = errorbar_groups(noiseAll', abs(noiseAll-lowerCI(:, 1:4))', abs(noiseAll-upperCI(:, 1:4))', ...
                'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx1);
set(hPanel, 'EdgeColor', 'none')
set(gca, 'FontSize', fontSize)

hAx2 = subplot(1, 2, 2);
[~, hPanel] = errorbar_groups(priorRange', abs(priorRange-lowerCI(:, 5))', abs(priorRange-upperCI(:, 5))', ...
                'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx2);
set(hPanel, 'EdgeColor', 'none')
set(gca, 'FontSize', fontSize)

%% Plot the scatter plot of motor noise
% Exp 1
stdFitExp1 =  [1.94 1.24 3.06 3.36 4.49 3.76];
stdDataExp1 = [2.76 2.32 3.87 9.29 8.37 5.89];
stdDataExp1_1 = [2.96 2.42 4.01 8.18 10.17 6.72];
stdDataExp1_2 = [2.56 2.24 3.72 10.29 6.13 4.41];

hBar = figure;
figPos = [400, -40, 800, 800];
set(hBar,'Units','pixels','Position',figPos)
hold on
set(gca, 'FontSize', 23)
bar([stdFitExp1' stdDataExp1' stdDataExp1_1' stdDataExp1_2']);
set(gca,'XTickLabel',{'1','2','3','4','5','6'}, 'XTick', 1:6) 
xlabel('Subject')
ylabel('Motor std (deg)')
legend('Fit', 'Data (all)', 'Data (first half)', 'Data (second half)', 'Location', 'NorthWest')
title('Experiment 1')

% Exp 2, 3
stdFitExp23 = [2.87 4.22 3.65 2.98 4.99 4.45];
stdFitExp2 = [3.19 3.06 3.12 2.49 2.42 3.59];
stdDataExp23 = [2.76 3.69 3.69 3.5 2.49 3.25];
stdDataExp23_1 = [2.96 4.14 3.64 3.69 2.54 3.45];
stdDataExp23_2 = [2.56 3.17 3.74 3.3 2.39 2.95];

hBar = figure;
figPos = [400, -40, 800, 800];
set(hBar,'Units','pixels','Position',figPos)
hold on
set(gca, 'FontSize', 23)
bar([stdFitExp23' stdFitExp2' stdDataExp23' stdDataExp23_1' stdDataExp23_2']);
set(gca,'XTickLabel',{'1','2','3','4','5','6'}, 'XTick', 1:6) 
xlabel('Subject')
ylabel('Motor std (deg)')
legend('Fit Exp 2,3', 'Fit Exp 2', 'Data (all)', 'Data (first half)', 'Data (second half)', 'Location', 'NorthWest')
title('Experiment 2, 3')


%% Plot Exp 1 (Old)         
totalNoiseExp1 = sqrt(paramModel(:, 1:3).^2 + repmat(paramModel(:, 6).^2, 1, 3)); 
totalNoiseExp1 = reshape (totalNoiseExp1', 1, []);
memoryNoiseExp1 = paramModel(:, 6);
memoryNoiseExp1 = reshape(repmat(memoryNoiseExp1', 3, 1), 1, []);
totalNoiseExp1New = NaN(numel(totalNoiseExp1), 2);
for ii = 1 : numel(totalNoiseExp1)
    totalNoiseExp1New(ii, 1) = memoryNoiseExp1(ii);
    totalNoiseExp1New(ii, 2) = totalNoiseExp1(ii);
end
totalNoiseExp1 = reshape([reshape(totalNoiseExp1New,3,[]); zeros(1,numel(totalNoiseExp1New)/3)],[],2);
totalNoiseExp1(:, 2) = totalNoiseExp1(:, 2) - totalNoiseExp1(:, 1);

% confInternvalNoise = confidenceIntervalAll(1:3,1:6);
% confInternvalNoise = confInternvalNoise(:);
% confInternvalNoise = cell2mat(confInternvalNoise);
% confInternvalNoise = reshape([reshape(confInternvalNoise,3,[]); zeros(1,numel(confInternvalNoise)/3)],[],2);
priorRange = paramModel(:, 5);
% confInternvalPrior = confidenceIntervalAll(4,1:6);
% confInternvalPrior = cell2mat(confInternvalPrior');
% alpha = alphaCompute(priorRange, paramModelExp1(:, 8));
% confInternvalAlpha = confidenceIntervalAll(5,1:6);
% confInternvalAlpha = cell2mat(confInternvalAlpha');

% Plot the parameters with bars
hNoiseExp1 = figure;
fontSize = 30;
hold on
set(gca, 'FontSize', fontSize)
hPanel = bar(totalNoiseExp1, 'stack');
set(hPanel, 'EdgeColor', 'none')
% errorbar(1:length(totalNoiseExp1), sum(totalNoiseExp1, 2),  sum(totalNoiseExp1, 2)-confInternvalNoise(:,1)...
%                                                             , zeros(length(totalNoiseExp1), 1), 'white');
% errorbar(1:length(totalNoiseExp1), sum(totalNoiseExp1, 2),  zeros(length(totalNoiseExp1), 1)...
%                                                             , confInternvalNoise(:, 2)-sum(totalNoiseExp1, 2), 'k');                                        
% set(gca, 'XTick', [])
% ylabel('Noise (deg)')
% ylim([0 max(confInternvalNoise(:, 2))])

hPriorExp1 = figure;
hold on
set(gca, 'FontSize', fontSize)
hPanel = bar(priorRange);
% errorbar(1:length(priorRange), priorRange,  priorRange-confInternvalPrior(:,1)...
%                                 , confInternvalPrior(:, 2)-priorRange, 'k');
set(hPanel, 'EdgeColor', 'none')
ylabel('Prior range (deg)')
set(gca, 'XTick', [])

%% Plot Exp 2 and 3         
totalNoiseExp2 = sqrt(paramModelExp2(:, 1:3).^2 + repmat(paramModelExp2(:, 6).^2, 1, 3)); 
totalNoiseExp2 = reshape (totalNoiseExp2', 1, []);
totalNoiseExp2Plot = totalNoiseExp2;
memoryNoiseExp2 = paramModelExp2(:, 6);
memoryNoiseExp2 = reshape(repmat(memoryNoiseExp2', 3, 1), 1, []);
priorRange = paramModelExp2(:, 5);
totalNoiseExp2New = NaN(numel(totalNoiseExp2), 2);
for ii = 1 : numel(totalNoiseExp2)
    totalNoiseExp2New(ii, 1) = memoryNoiseExp2(ii);
    totalNoiseExp2New(ii, 2) = totalNoiseExp2(ii);
end
totalNoiseExp2 = reshape([reshape(totalNoiseExp2New,3,[]); zeros(1,numel(totalNoiseExp2New)/3)],[],2);
totalNoiseExp2(:, 2) = totalNoiseExp2(:, 2) - totalNoiseExp2(:, 1);
% confInternvalNoise = confidenceIntervalAll(1:3,1:6);
% confInternvalNoise = confInternvalNoise(:);
% confInternvalNoise = cell2mat(confInternvalNoise);
% confInternvalNoise = reshape([reshape(confInternvalNoise,3,[]); zeros(1,numel(confInternvalNoise)/3)],[],2);
% confInternvalPrior = confidenceIntervalAll(4,1:6);
% confInternvalPrior = cell2mat(confInternvalPrior');

% Plot the parameters with bars
hNoiseExp2 = figure;
fontSize = 30;
hold on
set(gca, 'FontSize', fontSize)
hPanel = bar(totalNoiseExp2, 'stack');
set(hPanel, 'EdgeColor', 'none')
% errorbar(1:length(totalNoiseExp2), sum(totalNoiseExp2, 2),  sum(totalNoiseExp2, 2)-confInternvalNoise(:,1)...
%                                                             , zeros(length(totalNoiseExp2), 1), 'white');
% errorbar(1:length(totalNoiseExp2), sum(totalNoiseExp2, 2),  zeros(length(totalNoiseExp2), 1)...
%                                                             , confInternvalNoise(:, 2)-sum(totalNoiseExp2, 2), 'k');                                        
set(gca, 'XTick', [])
ylabel('Noise (deg)')
% ylim([0 max(confInternvalNoise(:, 2))])

hPriorExp2 = figure;
hold on
set(gca, 'FontSize', fontSize)
hPanel = bar(priorRange);
% errorbar(1:length(priorRange), sum(priorRange, 2),  sum(priorRange, 2)-confInternvalPrior(:,1)...
%                                                             , zeros(length(priorRange), 1), 'white');
% errorbar(1:length(priorRange), sum(priorRange, 2),  zeros(length(priorRange), 1)...
%                                                             , confInternvalPrior(:, 2)-sum(priorRange, 2), 'k');                                                                    
set(hPanel, 'EdgeColor', 'none')
ylabel('Prior range (deg)')
set(gca, 'XTick', [])
% ylim([0 max(confInternvalPrior(:, 2))])

% % Plot the parameters scatter plot
% colorName = {'g', 'r', 'b', 'c', 'm', 'y'};
% hNoise = figure;
% hold on
% fontSize = 20;
% for ii = 1 : 5
%     set(gca, 'FontSize', fontSize)
%     ciExp2 = [confidenceIntervalAll{:, 5+ii}]';
%     ciExp3 = [confidenceIntervalAll{:, 10+ii}]';
%     errorbarxy(totalNoise(:, 5+ii), totalNoise(:, 10+ii), ...
%                 totalNoise(:, 5+ii)-ciExp2([1; 3; 5]), ciExp2([2; 4; 6])-totalNoise(:, 5+ii),...
%                 totalNoise(:, 10+ii)-ciExp3([1; 3; 5]), ciExp3([2; 4; 6])-totalNoise(:, 10+ii),...
%                 {[colorName{ii} '.'], colorName{ii}, colorName{ii}});
%     xlabel('Noise Exp 2 (deg)')
%     ylabel('Noise Exp 3 (deg)')
%     maxX = max(totalNoise(:))+1;
%     plot([0 maxX], [0 maxX], 'k--')
%     axis equal
%     xlim([0 maxX])
%     ylim([0 maxX])
%  end
% 
% hPrior = figure;
% hold on
% for ii = 1 : 5
%     set(gca, 'FontSize', fontSize)
%     ciExp2 = [confidenceIntervalAll{:, 5+ii}]';
%     ciExp3 = [confidenceIntervalAll{:, 10+ii}]';
%     errorbarxy(priorRange(2, ii), priorRange(3, ii),...
%                 priorRange(2, ii)-ciExp2(7), ciExp2(8)-priorRange(2, ii),...
%                 priorRange(3, ii)-ciExp3(7), ciExp3(8)-priorRange(3, ii),...
%                 {[colorName{ii} '.'], colorName{ii}, colorName{ii}});
%     xlabel('Prior Exp 2 (deg)')
%     ylabel('Prior Exp 3 (deg)')
%     maxX = max(priorRange(:))+1;
%     plot([0 maxX], [0 maxX], 'k--')
%     axis equal    
%     xlim([0 maxX])
%     ylim([0 maxX])
% end

%% Plot incongruent error for combined subject
incongruentErrorExp1_ave = mean([0.1111 0.7778 5.7778 4.6111 3.8333]); % ave:  3.0222
lapseExp1_ave = mean([0 0.76 4.75 1.88 2.86]);
motorExp1_ave = mean([1.3369 0.11837 0.044024 5.0207 1.5351] .* (100 - lapseExp1_ave) / 100);
incongruentErrorExp1_predicted_ave = lapseExp1_ave + motorExp1_ave;

incongruentErrorExp2_ave = mean([0.0556 3.2778 14.4444 2.8333 4.3333]);% ave: 4.9889
lapseExp2_ave = mean([0.13 5.06 13.6 0.58 2.85]);
motorExp2_ave = mean([1.3435 0.15152 1.1067 0.40338 0.0038875] .* (100 - lapseExp2_ave) / 100);
incongruentErrorExp2_predicted_ave = lapseExp2_ave + motorExp2_ave;

incongruentErrorExp3_ave = mean([0.0556 3.8333 6.5556  5 2.6111]);% ave:  3.6111
motorExp3_ave = mean([2.0946 0.66036 1.9753 1.0927 0.0056516]);
errorColorExp3_ave = mean([0.055556 5.6667 0.22222 1.8889 2.0556]);
incongruentErrorExp3_predicted_ave = mean([2.1880    5.3250   12.6090    1.5748    2.8526]);

errorPlot = [incongruentErrorExp1_ave 0;
             motorExp1_ave lapseExp1_ave;
             incongruentErrorExp2_ave 0;
             motorExp2_ave lapseExp2_ave;
             incongruentErrorExp3_ave 0;
             motorExp3_ave incongruentErrorExp3_predicted_ave-motorExp3_ave];

figure
bar(errorPlot, 'stack')

%% Plot incongruent error for individual subject (separate for each Experiment)
incongruentErrorExp1 = [0.1111 0.7778 5.7778 4.6111 3.8333]; % ave:  3.0222
lapseExp1 = [0 0.76 4.75 1.88 2.86];
motorExp1 = [1.3369 0.11837 0.044024 5.0207 1.5351] .* (100 - lapseExp1) / 100;
incongruentErrorExp1_predicted = lapseExp1 + motorExp1;

incongruentErrorExp2 = [0.0556 3.2778 14.4444 2.8333 4.3333];% ave: 4.9889
lapseExp2 = [0.13 5.06 13.6 0.58 2.85];
motorExp2 = [1.3435 0.15152 1.1067 0.40338 0.0038875] .* (100 - lapseExp2) / 100;
incongruentErrorExp2_predicted = lapseExp2 + motorExp2;

incongruentErrorExp3 = [0.0556 3.8333 6.5556  5 2.6111];% ave:  3.6111
motorExp3 = [2.0946 0.66036 1.9753 1.0927 0.0056516];
errorColorExp3 = [0.055556 5.6667 0.22222 1.8889 2.0556];
incongruentErrorExp3_predicted = [2.1880    5.3250   12.6090    1.5748    2.8526];

incongruentErrorExp1_ave = mean([0.1111 0.7778 5.7778 4.6111 3.8333]); % ave:  3.0222
lapseExp1_ave = mean([0 0.76 4.75 1.88 2.86]);
motorExp1_ave = mean([1.3369 0.11837 0.044024 5.0207 1.5351] .* (100 - lapseExp1_ave) / 100);
incongruentErrorExp1_predicted_ave = lapseExp1_ave + motorExp1_ave;

incongruentErrorExp2_ave = mean([0.0556 3.2778 14.4444 2.8333 4.3333]);% ave: 4.9889
lapseExp2_ave = mean([0.13 5.06 13.6 0.58 2.85]);
motorExp2_ave = mean([1.3435 0.15152 1.1067 0.40338 0.0038875] .* (100 - lapseExp2_ave) / 100);
incongruentErrorExp2_predicted_ave = lapseExp2_ave + motorExp2_ave;

incongruentErrorExp3_ave = mean([0.0556 3.8333 6.5556  5 2.6111]);% ave:  3.6111
motorExp3_ave = mean([2.0946 0.66036 1.9753 1.0927 0.0056516]);
errorColorExp3_ave = mean([0.055556 5.6667 0.22222 1.8889 2.0556]);
incongruentErrorExp3_predicted_ave = mean([2.1880    5.3250   12.6090    1.5748    2.8526]);


figure;
colorName = {'SlateGray', 'Crimson', 'DarkMagenta', 'DarkOrange', 'DarkGoldenRod', 'SpringGreen', 'Teal', 'DodgerBlue', 'Navy'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

subplot(1, 2, 1)
hold on
set(gca, 'FontSize', 20)
scatter(incongruentErrorExp3_predicted, incongruentErrorExp3, 13^2, colorIndex([1 6:9], :), 'filled');
plot([motorExp3_ave incongruentErrorExp3_predicted_ave], [incongruentErrorExp3_ave incongruentErrorExp3_ave], 'kd--') 
plot([-0.5 15], [-0.5 15], '--')
xlim([-0.5 15])
ylim([-0.5 15])
xlabel('% incongruent (predicted-motor+lapseExp2)')
ylabel('% incongruent (data)')
title(['Lapse + motor, r: ', num2str(corr(incongruentErrorExp3_predicted', incongruentErrorExp3'))])

subplot(1, 2, 2)
hold on
set(gca, 'FontSize', 20)
scatter(motorExp3 + errorColorExp3, incongruentErrorExp3, 13^2, colorIndex([1 6:9], :), 'filled');
plot([motorExp3_ave motorExp3_ave + errorColorExp3_ave], [incongruentErrorExp3_ave incongruentErrorExp3_ave], 'kd--') 
plot([-0.5 15], [-0.5 15], '--')
xlim([-0.5 15])
ylim([-0.5 15])
xlabel('% incongruent (predicted-motor+lapseExp2)')
ylabel('% incongruent (data)')
title(['Color + motor, r: ', num2str(corr((motorExp3 + errorColorExp3)', incongruentErrorExp3'))])


subplot(1, 2, 2)
hold on
set(gca, 'FontSize', 20)
scatter(incongruentErrorExp1_predicted, incongruentErrorExp1, 13^2, colorIndex(1:5, :), 'filled');
plot([motorExp1_ave incongruentErrorExp1_predicted_ave], [incongruentErrorExp1_ave incongruentErrorExp1_ave], 'kd--') 
plot([-0.5 15], [-0.5 15], '--')
xlim([-0.5 15])
ylim([-0.5 15])
xlabel('% incongruent (predicted)')
ylabel('% incongruent (data)')
title(['Exp 1, r: ', num2str(corr(incongruentErrorExp1_predicted', incongruentErrorExp1'))])

subplot(1, 3, 2)
hold on
set(gca, 'FontSize', 20)
scatter(incongruentErrorExp2_predicted, incongruentErrorExp2, 13^2, colorIndex([1 6:9], :), 'filled');
plot([motorExp2_ave incongruentErrorExp2_predicted_ave], [incongruentErrorExp2_ave incongruentErrorExp2_ave], 'kd--') 
plot([-0.5 15], [-0.5 15], '--')
xlim([-0.5 15])
ylim([-0.5 15])
xlabel('% incongruent (predicted)')
ylabel('% incongruent (data)')
title(['Exp 2, r: ', num2str(corr(incongruentErrorExp2_predicted', incongruentErrorExp2'))])

subplot(1, 3, 3)
hold on
set(gca, 'FontSize', 20)
scatter(incongruentErrorExp3_predicted, incongruentErrorExp3, 13^2, colorIndex([1 6:9], :), 'filled');
plot([motorExp3_ave incongruentErrorExp3_predicted_ave], [incongruentErrorExp3_ave incongruentErrorExp3_ave], 'kd--') 
plot([-0.5 15], [-0.5 15], '--')
xlim([-0.5 15])
ylim([-0.5 15])
xlabel('% incongruent (predicted-motor+lapseExp2)')
ylabel('% incongruent (data)')
title(['Exp 3, r: ', num2str(corr(incongruentErrorExp3_predicted', incongruentErrorExp3'))])

% subplot(2, 3, 1)
% hold on
% set(gca, 'FontSize', 20)
% scatter(lapseExp1, incongruentErrorExp1, 13^2, colorIndex(1:5, :), 'filled');
% plot([-0.5 15], [-0.5 15], '--')
% xlim([-0.5 15])
% ylim([-0.5 15])
% xlabel('Lapse rate')
% ylabel('% incongruent (data)')
% title(['Exp 1, r: ', num2str(corr(lapseExp1', incongruentErrorExp1'))])
% 
% subplot(2, 3, 2)
% hold on
% set(gca, 'FontSize', 20)
% scatter(lapseExp2, incongruentErrorExp2, 13^2, colorIndex([1 6:9], :), 'filled');
% plot([-0.5 15], [-0.5 15], '--')
% xlim([-0.5 15])
% ylim([-0.5 15])
% xlabel('Lapse rate')
% ylabel('% incongruent (data)')
% title(['Exp 2, r: ', num2str(corr(lapseExp2', incongruentErrorExp2'))])
% 
% subplot(2, 3, 3)
% hold on
% set(gca, 'FontSize', 30)
% nullX = zeros(1, 9);
% nullY = zeros(1, 9);
% h = NaN(1, 9);
% for ii = 1 : 9
%     h(ii) = scatter(nullX(ii), nullY(ii), 1, colorIndex(ii, :), 'filled');
% end
% [~, iconHandle]  = legend(h, 'Subject 1', 'Subject 2', 'Subject 3', 'Subject 4', 'Subject 5', 'Subject 6', 'Subject 7', 'Subject 8', 'Subject 9', 'Location', 'NorthWest');
% legend('boxoff')
% for ii = 1 : 9
%     iconHandle(ii+9).Children.MarkerSize = 15;
% end
% set(gca, 'FontSize', 20)
% scatter(motorExp3, incongruentErrorExp3, 13^2, colorIndex([1 6:9], :), 'filled');
% maxPlot = max([motorExp3 incongruentErrorExp3]);
% plot([-0.5 15], [-0.5 15], '--')
% xlim([-0.5 15])
% ylim([-0.5 15])
% xlabel('% incongruent (predicted-motor)')
% ylabel('% incongruent (data)')
% title(['Exp 3, r: ', num2str(corr(motorExp3', incongruentErrorExp3'))])


% figure
% hold on
% set(gca, 'FontSize', 20)
% scatter(errorColorExp3+motorExp3, incongruentErrorExp3, 13^2, colorIndex([1 6:9], :), 'filled');
% maxPlot = max([errorColorExp3 incongruentErrorExp3]);
% plot([-0.5 15], [-0.5 15], '--')
% xlim([-0.5 15])
% ylim([-0.5 15])
% xlabel('% incongruent (predicted)')
% ylabel('% incongruent (data)')
% title(['Exp 3, r: ', num2str(corr(errorColorExp3'+motorExp3', incongruentErrorExp3'))])

% figure
% bar(inErrorExp1)
% figure
% bar(inErrorExp2)
% figure
% bar(inErrorExp3)

%% Plot incongruent error for individual subject (combined for all experiments)
incongruentErrorExp1 = [0.1111 0.7778 5.7778 4.6111 3.8333]; % ave:  3.0222
lapseExp1 = [0 0.76 4.75 1.88 2.86];
motorExp1 = [1.3369 0.11837 0.044024 5.0207 1.5351] .* (100 - lapseExp1) / 100;
incongruentErrorExp1_predicted = lapseExp1 + motorExp1;

incongruentErrorExp2 = [0.0556 3.2778 14.4444 2.8333 4.3333];% ave: 4.9889
lapseExp2 = [0.13 5.06 13.6 0.58 2.85];
motorExp2 = [1.3435 0.15152 1.1067 0.40338 0.0038875] .* (100 - lapseExp2) / 100;
incongruentErrorExp2_predicted = lapseExp2 + motorExp2;

incongruentErrorExp3 = [0.0556 3.8333 6.5556  5 2.6111];% ave:  3.6111
motorExp3 = [2.0946 0.66036 1.9753 1.0927 0.0056516];
errorColorExp3 = [0.055556 5.6667 0.22222 1.8889 2.0556];
incongruentErrorExp3_predicted = [2.1880    5.3250   12.6090    1.5748    2.8526];

incongruentErrorExp1_ave = mean([0.1111 0.7778 5.7778 4.6111 3.8333]); % ave:  3.0222
lapseExp1_ave = mean([0 0.76 4.75 1.88 2.86]);
motorExp1_ave = mean([1.3369 0.11837 0.044024 5.0207 1.5351] .* (100 - lapseExp1_ave) / 100);
incongruentErrorExp1_predicted_ave = lapseExp1_ave + motorExp1_ave;

incongruentErrorExp2_ave = mean([0.0556 3.2778 14.4444 2.8333 4.3333]);% ave: 4.9889
lapseExp2_ave = mean([0.13 5.06 13.6 0.58 2.85]);
motorExp2_ave = mean([1.3435 0.15152 1.1067 0.40338 0.0038875] .* (100 - lapseExp2_ave) / 100);
incongruentErrorExp2_predicted_ave = lapseExp2_ave + motorExp2_ave;

incongruentErrorExp3_ave = mean([0.0556 3.8333 6.5556  5 2.6111]);% ave:  3.6111
motorExp3_ave = mean([2.0946 0.66036 1.9753 1.0927 0.0056516]);
errorColorExp3_ave = mean([0.055556 5.6667 0.22222 1.8889 2.0556]);
incongruentErrorExp3_predicted_ave = mean([2.1880    5.3250   12.6090    1.5748    2.8526]);

incongruentErrorAll = [incongruentErrorExp1 incongruentErrorExp2 incongruentErrorExp3];
incongruentErrorAll_predicted = [incongruentErrorExp1_predicted incongruentErrorExp2_predicted incongruentErrorExp3_predicted];


figure;
colorName = {'SlateGray', 'Crimson', 'DarkMagenta', 'DarkOrange', 'DarkGoldenRod', 'SpringGreen', 'Teal', 'DodgerBlue', 'Navy'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

hold on
set(gca, 'FontSize', 20)
scatter(incongruentErrorExp1_predicted, incongruentErrorExp1, 13^2, colorIndex(1:5, :), 'filled');
% plot([motorExp1_ave incongruentErrorExp1_predicted_ave], [incongruentErrorExp1_ave incongruentErrorExp1_ave], 'kd--') 

scatter(incongruentErrorExp2_predicted, incongruentErrorExp2, 13^2, colorIndex([1 6:9], :), 'square', 'filled');
% plot([motorExp2_ave incongruentErrorExp2_predicted_ave], [incongruentErrorExp2_ave incongruentErrorExp2_ave], 'ks--') 

scatter(incongruentErrorExp3_predicted, incongruentErrorExp3, 13^2, colorIndex([1 6:9], :), 'd', 'filled');
% plot([motorExp3_ave incongruentErrorExp3_predicted_ave], [S1incongruentErrorExp3_ave incongruentErrorExp3_ave], 'kx--') 
plot([-0.5 15], [-0.5 15], '--')
xlim([-0.5 15])
ylim([-0.5 15])
xlabel('% incongruent (predicted)')
ylabel('% incongruent (data)')
title(['All experiment, r: ', num2str(corr(incongruentErrorAll', incongruentErrorAll_predicted'))])

%% Compute the average params value
subjID = {'Ave23'}; % {'ll', 'sy', 'cz', 'vs', 'as',  'average', 'll', 'xfl', 'aj', 'zw', 'skb', 'average'};
paramsAll = cell(1, length(subjID));
if strcmp(subjID{1}, 'LL')
    paramsMean = NaN(length(subjID), 10);
else
    paramsMean = NaN(length(subjID), 9);
end
figure;
hold on
for kk = 1 : length(subjID);
    fileName = ['FitResult-' subjID{kk} '.txt'];
    fileID = fopen(fileName);
    if strcmp(subjID{kk}, 'LL')
        paramSubj = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f','CommentStyle','//');
    else
        paramSubj = textscan(fileID,'%f %f %f %f %f %f %f %f %f','CommentStyle','//');
    end
    paramsAll{kk} = cell2mat(paramSubj);
    paramsMean(kk, :) = nanmean(paramsAll{kk}, 1);
    paramSubjMat = cell2mat(paramSubj);
    [~, ind] = sort(paramSubjMat(:, 1), 'ascend');
    paramSort = paramSubjMat(ind, :);
    fclose(fileID);

    if strcmp(subjID{kk}, 'LL')
        plot(kk*ones(length(paramSubj{8}), 1), paramSubj{8}, 'o');    
        plot(kk, nanmean(paramSubj{8}), '*', 'MarkerSize', 10); 
        fprintf('%9.2f %9.4f %9.4f %9.4f  %10.4f %10.4f %8.4f %9.4f %9.4f %9.4f\r\n', ...
            paramsMean(1), paramsMean(2), paramsMean(3), paramsMean(4), paramsMean(5), paramsMean(6), paramsMean(7), paramsMean(8), paramsMean(9), paramsMean(10));            
    else        
        plot(kk*ones(length(paramSubj{7}), 1), paramSubj{7}, 'o');    
        plot(kk, nanmean(paramSubj{7}), '*', 'MarkerSize', 10); 
        fprintf('%9.2f %9.4f %9.4f %9.4f  %10.4f %10.4f %8.4f %9.4f %9.4f\r\n', ...
            paramsMean(1), paramsMean(2), paramsMean(3), paramsMean(4), paramsMean(5), paramsMean(6), paramsMean(7), paramsMean(8), paramsMean(9));            
    end
end
xlim([0 length(subjID)+2])
paramSort(1,2:end)

% subjID = {'LL' };
% paramsAll = cell(1, length(subjID));
% paramsMean = NaN(length(subjID), 9);
% figure;
% hold on
% for kk = 1 : length(subjID);
%     fileName = ['FitResult-' subjID{kk} '.txt'];
%     fileID = fopen(fileName);
%     myFile = textscan(fileID,'%s','delimiter','\n');
%     myFile = myFile{1};
%     saveNextLine = 0;
%     counter = 1;
%     if strcmp(subjID{kk}, 'LL')
%         paramsEachSubj = NaN(30, 10);
%     else
%         paramsEachSubj = NaN(30, 9);
%     end
%     for ii = 1 : length(myFile)
%         tempLine = myFile{ii};
%         if saveNextLine
%             paramsEachSubj(counter, :) = str2num(tempLine);
%             counter = counter + 1;
%             saveNextLine = 0;
%         else
%             indMatch = strfind(tempLine,'Iteration');
%             if ~isempty(indMatch)
%                 saveNextLine = 1;
%             end
%         end
%     end
%     paramsAll{kk} = paramsEachSubj;
%     paramsMean(kk, :) = nanmean(paramsAll{kk}, 1);    
%     
%     plot(kk*ones(length(paramsEachSubj(:, 7)), 1), paramsEachSubj(:, 7), 'o');    
%     plot(kk, mean(paramsEachSubj(:, 7)), '*', 'MarkerSize', 10); 
%     fprintf('%9.2f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f %9.4f\r\n', ...
%         paramsMean(1), paramsMean(2), paramsMean(3), paramsMean(4), paramsMean(5), paramsMean(6), paramsMean(7), paramsMean(8), paramsMean(9));            
% end
% xlim([0 length(subjID)+2])

%% Plot the motor noise for individual subject
motorNoise = [2.63 1.84 2.91 5.4 5.89 2.93 2.99 2.53 2.05 3.6 2.5];
figure
bar(motorNoise)
xlabel('Subject')
ylabel('Standard deviation of motor noise (deg)')
box off