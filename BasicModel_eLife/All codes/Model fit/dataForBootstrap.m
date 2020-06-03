function [binaryDecision, percentCW, nTrialsPerCondition, estimateData, angleDiff, stdMotor] = dataForBootstrap(subjectID, expNumber, includeIncongruentTrials)
%% Extract decision and estimation data in main experiment
if strcmp(subjectID, 'average')
    experimentNumber = 1:5;
else
    experimentNumber = 1;
end
experimentType = 'MainExperiment';
if expNumber == 1
    experiment = 'ControlReplication';
    session = 1;
elseif expNumber == 2
    experiment = 'Original';
    if strcmp(subjectID, 'average')
        session = 2;
    else
        session = 1;
    end
elseif expNumber == 3
    experiment = 'DecisionGiven';
    if strcmp(subjectID, 'average')
        session = 2;
    else
        session = 1;
    end
end    
dataAll = [];
for ii = 1 : length(experimentNumber)
    dataFullPath = fullfile('Data', subjectID, experimentType, [experiment num2str(session)], ...
                                [experiment '-' num2str(experimentNumber(ii))]);          
    load(dataFullPath{1});
    if ii == 1
        dataAll = dataTotal;
    else
        dataAll = [dataAll; dataTotal];
    end
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

% Average across conditions
angleTrue = dataAll(:,5) - dataAll(:,1);
angleTrue(angleTrue>180 | angleTrue<0) = angleTrue(angleTrue>180 | angleTrue<0)...
                                            - 180 * sign(angleTrue(angleTrue>180 | angleTrue<0));
indexCorrect = find(abs(angleTrue-angleEstimate) > 70);
angleEstimate(indexCorrect) = angleEstimate(indexCorrect) + sign(angleTrue(indexCorrect)-angleEstimate(indexCorrect))*180;    
angleDiffEst = dataAll(:, 5) - angleEstimate;
indexAdjust = abs(angleDiffEst)>90;
angleDiffEst(indexAdjust) = angleDiffEst(indexAdjust) - 180 * sign(angleDiffEst(indexAdjust));
if strcmp(experiment, 'DecisionGiven')
    binaryDecisionAll = dataAll(:,8);
    indicatorConsistent = sign(angleDiffEst)==dataAll(:,8);
else
    binaryDecisionAll = dataAll(:,6);
    indicatorConsistent = sign(angleDiffEst)==dataAll(:,6);    
end
if ~isempty(includeIncongruentTrials)
    indicatorConsistent = ones(length(angleDiffEst), 1);
end
angleDiffEstimate = cell(length(stdNoiseLevel), length(angleDiff));
binaryDecisionSelf = cell(length(stdNoiseLevel), length(angleDiff));
binaryDecisionGiven = binaryDecisionSelf;

%%Extract data into angleDiffEstimate (each cell is estimate data for a condition)
for ii = 1 : length(stdNoiseLevel)
    for jj = 1 : length(angleDiff)
        indexSelect = dataAll(:,4)==stdNoiseLevel(ii) & dataAll(:,1)==angleDiff(jj);
        indexSelectLinear = find(indexSelect & indicatorConsistent);
        indexResample = indexSelectLinear(randi(length(indexSelectLinear),[1 length(indexSelectLinear)]));
        tempEst = angleDiffEst(indexResample);
        angleDiffEstimate{ii,jj}  = tempEst;
        if strcmp(experiment, 'DecisionGiven')
            binaryDecisionSelf{ii,jj} = sign(angleDiffEstimate{ii,jj});
            binaryDecisionGiven{ii,jj} = binaryDecisionAll(indexResample, :);     
            binaryDecision = binaryDecisionGiven;                            
        else
            tempDecision = binaryDecisionAll(indexResample, :);
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

%% Extract motor noise in motor experiment
if strcmp(subjectID, 'average')
    experimentNumber = 1:5;
elseif strcmp(subjectID, 'll')
    experimentNumber = 1:2;
end
experimentType = 'MotorNoise';
if expNumber == 1 || ~strcmp(subjectID, 'average')
    session = 1;
else
    session = 2;
end    
stdMotorAll = NaN(1, length(experimentNumber));

for ii = 1 : length(experimentNumber)
    dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                ['MotorNoise-' num2str(experimentNumber(ii))]);
    load(dataFullPath{1});

    % Extract data
    angleEstimate = dataTotal(:, 7);
%     nEstimate = sum(~isnan(angleEstimate));
%     params.trialOrder(nEstimate+1:end) = [];
%     angleEstimate(params.trialOrder(1:round(nEstimate/2))) = NaN;

    % Average across conditions
    angleTrue = dataTotal(:,5) - dataTotal(:,1);
    angleTrue(angleTrue>180 | angleTrue<0) = angleTrue(angleTrue>180 | angleTrue<0)...
                                                - 180 * sign(angleTrue(angleTrue>180 | angleTrue<0));
    indexCorrect = find(abs(angleTrue-angleEstimate) > 70);
    angleEstimate(indexCorrect) = angleEstimate(indexCorrect) + sign(angleTrue(indexCorrect)-angleEstimate(indexCorrect))*180;    
    angleDiffEst = dataTotal(:, 5) - angleEstimate;
    indexAdjust = abs(angleDiffEst)>90;
    angleDiffEst(indexAdjust) = angleDiffEst(indexAdjust) - 180 * sign(angleDiffEst(indexAdjust));
    biasAll = angleDiffEst - dataTotal(:, 1);
    stdRobust = NaN(1, length(angleDiff));
    for jj = 1 : length(angleDiff)
        indexSelect = dataTotal(:,1)==angleDiff(jj);
        tempBias = biasAll(indexSelect, :);
%         robustFit = mcdcov(tempBias, 'plots', 0, 'alpha', 0.8);
%         stdRobust(jj) = sqrt(robustFit.cov);
        % delete outliers that are 3 MAD away from the median MAD = median(|x_i - median(x)|)
        MAD = median(abs(tempBias - median(tempBias)));
        indDelete = abs(tempBias - median(tempBias)) > 3.5 * MAD;
        tempBias(indDelete) = [];
        stdRobust(jj) = nanstd(tempBias);
    end    

    stdMotorAll(ii) = mean(stdRobust);
end
stdMotor = mean(stdMotorAll);
end