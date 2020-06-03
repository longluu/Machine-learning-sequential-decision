%%%%%%%% Data analysis of conditional decision making experiment %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
% 
subjectIDAll = {'Sc1', 'Sc2', 'll', 'sy', 'cz', 'vs', 'as', 'll', 'xfl', 'aj', 'zw', 'skb'};
percentInconsistentSplitByNoise = NaN(length(subjectIDAll), 3);
colorName = {'Black', 'Brown', 'RosyBrown', 'SlateGray', 'Crimson', 'DarkMagenta', 'DarkOrange', 'DarkGoldenRod', 'SpringGreen', 'Teal', 'DodgerBlue', 'Navy'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
h1 = figure;
for kk = 1 : length(subjectIDAll)
    percentInconsistentSplitByAngle = NaN(8, 3);
    subjectID = subjectIDAll{kk};
    experimentNumber = 1;
    if kk <8
        experiment = 'ControlReplication';
    else
        experiment = 'Original';
    end
    if strcmp(subjectIDAll{kk}, 'Sc1')        
        subjectID = 'average';
        experimentNumber = 1:5;
    elseif strcmp(subjectIDAll{kk}, 'Sc2') 
        subjectID = 'average';
        experimentNumber = 1:5;
        experiment = 'Original';
    end
    experimentType = 'MainExperiment';
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
    indicatorInconsistent = sign(angleDiffEst)~=dataAll(:,6);    
    for ii = 1 : length(stdNoiseLevel)
        indexSelect = dataAll(:,4)==stdNoiseLevel(ii);
        indexInconsistent = indicatorInconsistent(indexSelect, :);  
        percentInconsistentSplitByNoise(kk, ii) = 100 * nansum(indexInconsistent) / sum(~isnan(indexInconsistent));
    end
    rangeCollapse = round(length(angleDiff)/2);
    for ii = rangeCollapse : length(angleDiff)
        for jj = 1 : length(stdNoiseLevel)
            indexSelect = abs(dataAll(:,1))==angleDiff(ii) & dataAll(:,4)==stdNoiseLevel(jj);
            indexInconsistent = indicatorInconsistent(indexSelect, :);  
            percentInconsistentSplitByAngle(ii-rangeCollapse+1, jj) = 100 * nansum(indexInconsistent) / sum(~isnan(indexInconsistent));
        end
    end
    
    % Plot
    if kk < 3
        figure(h1)
        indexPlot = angleDiff >= 0;
        for jj = 1 : length(stdNoiseLevel)
            subplot(1, 3, jj)
            hold on
            plot(angleDiff(indexPlot), percentInconsistentSplitByAngle(:, jj), 'Color', colorIndex(kk, :))
        end
        xlim([min(angleDiff(indexPlot)) max(angleDiff(indexPlot))])
        xlabel('Stimulus orientation (deg)')
        ylabel('Percent inconsistent (%)')
    end
end

colorName = {'Black', 'Brown', 'RosyBrown', 'SlateGray', 'Crimson', 'DarkMagenta', 'DarkOrange', 'DarkGoldenRod', 'SpringGreen', 'Teal', 'DodgerBlue', 'Navy'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

hFig = figure;
hAx1 = gca;
[~, hPanel] = errorbar_groups(percentInconsistentSplitByNoise', zeros(size(percentInconsistentSplitByNoise')), zeros(size(percentInconsistentSplitByNoise')), ...
                'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx1);
set(hPanel, 'EdgeColor', 'none')
set(gca, 'FontSize', fontSize)
ylim([0 25])
xlabel('Subject')
ylabel('Percent inconsistent (%)')

% Plot the memory noise vs. inconsistency
paramsAllSubject = [3.4994    5.4317    9.5625      0.0000    25        6.0913    0.2953    0.01;
             4.6025    6.1782    8.9676      0.0000    23.1643   5.8385    0.8455    0.01;
             2.5611    4.8570    7.1153      0.0000    28.6016      0.8805    0.8461   0.01;
             4.0300    5.3505    8.1590      0.0000    41.7649   2.0957    0.0257    0.01;
             3.8692    6.4042   12.3882      0.0000    39.5010  14.0727    0.2322    0.01;
             3.7041    4.6566   10.9916      0.0000    36.4739   0.2433    0.6481    0.01;
             3.9226    5.7843   10.9330      0.0000    39.0534  14.2463    0.1492    0.01;
             2.5611    4.8570    7.1153      0.0000    23.1044    0.8805    0.8461    0.01;
             6.3630    8.4075   14.5464      0.0000    22.5672   5.9826    0.7652    0.01;
             4.6736    6.1762    7.9466      0.0000    15.8459   1.2144    0.3681    0.01;
             4.6585    5.5272    6.8023      0.0000    17.4750   4.1846    0.0934    0.01;
             4.3172    4.5344    8.4083      0.0000    33.9063  14.5556    0.5305    0.01];

memNoise = paramsAllSubject(3:end, 6);
percentInconsistentAll = mean(percentInconsistentSplitByNoise(3:end, :), 2);

figure;
colorName = {'SlateGray', 'Crimson', 'DarkMagenta', 'DarkOrange', 'DarkGoldenRod', 'SlateGray', 'SpringGreen', 'Teal', 'DodgerBlue', 'Navy'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
regressSlope = regress(percentInconsistentAll, [ones(length(memNoise), 1) memNoise]);
Xregress = [min(memNoise) max(memNoise)];
Yregress = Xregress * regressSlope(2) + regressSlope(1);

hold on
set(gca, 'FontSize', 20)
scatter(memNoise, percentInconsistentAll, 13^2, colorIndex, 'filled');
plot(Xregress, Yregress)
xlabel('Memory noise (deg)')
ylabel('Percent inconsistent (predicted)')
title([' r: ', num2str(corr(memNoise, percentInconsistentAll))])
