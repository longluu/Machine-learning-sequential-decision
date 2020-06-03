% Extract the estimation data
subjectID = 'average';
experimentNumber = 1:5;
experimentType = 'MainExperiment';
experiment = 'ControlReplication';
session = 1;
includeIncongruentTrials = 0;
dataAll = [];
fontSize = 20;
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
%[4.5068 2.8065 1.4884 -0.3927 -0.6928 -1.4532 -1.1594 0.2879 1.3119 0.3440 0.3761 -0.1966 -1.8427 -3.5352 -5.1920]; %tah
%[4.7655 3.3250 1.8600 1.1295 0.4500 -0.2327 -0.9759 0.4126 1.4950 0.8355 0.2337 -1.2497 -2.8566 -3.5273 -5.4698]; %rer
%[5.9333 5.2061 2.2371 0.3238 1.7823 -2.7771 -2.1744 -0.1529 2.7935 1.9921 1.3045 -2.1480 -3.1835 -4.2605 -5.5590]; %hp
%[0.6040 0.1296 -0.2591 -1.2071 -0.0587 -0.4571 -0.7466 -0.1348 0.4727 0.0008 -0.1369 0.5066 0.1502 -0.2272 -0.2773]; %kn
%[3.4593 2.6668 0.6858 -0.0477 -0.8808 -0.9526 0.0765 0.2923 0.7379 1.0920 -0.3460 -0.0342 -1.6139 -1.8800 -3.3114]; %cz 
if strcmp(experimentType, 'MotorNoise')
    a = dataAll(:,7);
    dataAll(isnan(a),:) = [];    
end

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
dataZeroDiff = dataAll(dataAll(:,1) == 0,:);
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
biasAll = angleTrue - angleEstimate;
if strcmp(experiment, 'DecisionGiven')
    binaryDecisionAll = dataAll(:,8);
    indicatorConsistent = sign(angleDiffEst)==dataAll(:,8);
elseif strcmp(experimentType, 'MotorNoise') || strcmp(experiment, 'TrainArray')
    indicatorConsistent = ones(length(biasAll));
    binaryDecisionAll = dataAll(:,6);    
else
    binaryDecisionAll = dataAll(:,6);
    indicatorConsistent = sign(angleDiffEst)==dataAll(:,6);    
end
if includeIncongruentTrials == 1
    indicatorConsistentKeep = indicatorConsistent;
    indicatorConsistent = ones(length(biasAll), 1);
end
angleDiffEstimate = cell(length(stdNoiseLevel), length(angleDiff));
biasEstimate = cell(length(stdNoiseLevel), length(angleDiff));
binaryDecisionSelf = cell(length(stdNoiseLevel), length(angleDiff));
binaryDecisionGiven = binaryDecisionSelf;

%% Extract data into angleDiffEstimate (each cell is estimate data for a condition) including all trials
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
        else
            tempDecision = binaryDecisionAll(indexSelect, :);
            tempDecision(~indexConsistent) = NaN;   
            binaryDecisionSelf{ii,jj} = tempDecision;
            binaryDecision = binaryDecisionSelf;                            
        end
    end
end

% Collapse the data into estimateData
rangeCollapse = round(length(angleDiff)/2);
estimateData = cell(length(stdNoiseLevel), rangeCollapse);

for ii = 1 : length(stdNoiseLevel)
    estimateData{ii, 1} = angleDiffEstimate{ii, rangeCollapse};
    for jj =  2 : rangeCollapse
        tempEstimateCCW = -angleDiffEstimate{ii, rangeCollapse-jj+1};
        tempEstimateCW = angleDiffEstimate{ii, rangeCollapse+jj-1};
        estimateData{ii, jj} = [tempEstimateCCW; tempEstimateCW];
    end
end
optBinsCorrect = [];
save(estimateDataName, 'estimateData', 'angleDiff', 'optBinsCorrect')

%% Extract data into angleDiffEstimate, just incongruent trials
% for ii = 1 : length(stdNoiseLevel)
%     for jj = 1 : length(angleDiff)
%         indexSelect = dataAll(:,4)==stdNoiseLevel(ii) & dataAll(:,1)==angleDiff(jj);
%         tempBias = biasAll(indexSelect, :);
%         biasEstimate{ii,jj} = tempBias(:);                       
%         indexInconsistent = ~indicatorConsistentKeep(indexSelect, :);                                
%         tempEst = angleDiffEst(indexSelect);
%         tempEst(~indexInconsistent) = NaN;                                        
%         angleDiffEstimate{ii,jj}  = tempEst;
%         if strcmp(experiment, 'DecisionGiven')
%             binaryDecisionSelf{ii,jj} = sign(angleDiffEstimate{ii,jj});
%             binaryDecisionGiven{ii,jj} = binaryDecisionAll(indexSelect, :);     
%             binaryDecision = binaryDecisionGiven;                            
%         else
%             tempDecision = binaryDecisionAll(indexSelect, :);
%             tempDecision(~indexInconsistent) = NaN;   
%             binaryDecisionSelf{ii,jj} = tempDecision;
%             binaryDecision = binaryDecisionSelf;                            
%         end
%     end
% end
% 
% % Collapse the data into estimateDataIncongruentOnly
% rangeCollapse = round(length(angleDiff)/2);
% estimateDataIncongurent = cell(length(stdNoiseLevel), rangeCollapse);
% 
% for ii = 1 : length(stdNoiseLevel)
%     estimateDataIncongurent{ii, 1} = angleDiffEstimate{ii, rangeCollapse};
%     for jj =  2 : rangeCollapse
%         tempEstimateCCW = -angleDiffEstimate{ii, rangeCollapse-jj+1};
%         tempEstimateCW = angleDiffEstimate{ii, rangeCollapse+jj-1};
%         estimateDataIncongurent{ii, jj} = [tempEstimateCCW; tempEstimateCW];
%     end
% end
% optBinsCorrect = [];
% save('estimateDataExp2IncongruentOnly', 'estimateDataIncongurent', 'angleDiff', 'optBinsCorrect')

% Get the discrimination data
percentCW = NaN(length(stdNoiseLevel),length(params.barAngleDiff));
nTrialsPerCondition = NaN(length(stdNoiseLevel),length(params.barAngleDiff));
if strcmp(experiment, 'DecisionGiven')
    binaryDecisionAll = dataAll(:, 8);
else
    binaryDecisionAll = dataAll(:, 6);
end    
for ii = 1 : length(stdNoiseLevel)
    for jj = 1 : length(percentCW)
        indexAll = dataAll(:,1) == angleDiff(jj) ...
                    & dataAll(:,4) == stdNoiseLevel(ii)...
                    & indicatorConsistent;
        indexCW = dataAll(:,1) == angleDiff(jj) ...
                    & dataAll(:,4) == stdNoiseLevel(ii)...
                    & indicatorConsistent ...
                    & binaryDecisionAll == 1;
        nTrialsPerCondition(ii,jj) = nansum((indexAll));
        percentCW(ii,jj) = 100*nansum(indexCW)/nansum(indexAll);
    end
end
save(discriminationDataName, 'binaryDecision', 'percentCW', 'nTrialsPerCondition')
