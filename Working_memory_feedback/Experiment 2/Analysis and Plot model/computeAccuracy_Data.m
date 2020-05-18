%%%%%%%%%%%%%%%%%%%%%%% Compute the LLH of the models %%%%%%%%%%%%%%%%%%%%%%

%% Compute the subject LLH
subjectIDAll = {'ll', 'pw', 'eh', 'bh', 'ln', 'at', 'dh'}; % 
errorDecision = NaN(1, length(subjectIDAll));
errorEstimate = NaN(1, length(subjectIDAll));

for nn = 1 : length(subjectIDAll)
    subjectID = subjectIDAll{nn};
    if strcmp(subjectID, 'average')
        experimentNumber = 1:7;
    else
        experimentNumber = 1;    
    end
    experimentType = 'MainExperiment';
    experiment = 'Original';
    session = 1;
    dataAll = [];
    fontSize = 20;
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

    % Extract data
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
    binaryFeedbackAll = dataAll(:,8);  
    binaryDecisionAll = dataAll(:,6);  
    
    %% Compute the decision error
    errorDecision(nn) = nansum(binaryFeedbackAll ~= binaryDecisionAll) / length(binaryFeedbackAll);
    errorEstimate(nn) = nansum(abs(angleTrue - angleEstimate)) / length(angleTrue);

end

errorAll = [errorDecision; errorEstimate; errorDecision+errorEstimate];

%% Plot the LLH
colorName = {'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

hFig = figure;
[~, hPanel] = errorbar_groups(errorAll, zeros(size(errorAll)), zeros(size(errorAll)), ...
                'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig);
ylabel('Error')
xlabel('Subject') 
set(gca, 'FontSize', 20, 'XTickLabel', {'1'; '2'; '3'; '4'; '5'; '6'; '7'})
title('Incorrect trials')
legend('Decision', 'Estimate', 'Total', 'Location', 'NorthWest')