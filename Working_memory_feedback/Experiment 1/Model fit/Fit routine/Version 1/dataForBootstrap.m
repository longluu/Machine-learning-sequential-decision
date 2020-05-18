function [binaryDecision, percentCW, nTrialsPerCondition, estimateData, angleDiff, stdMotor] = dataForBootstrap(subjectID, includeIncongruentTrials, correctIndex)
%% Extract decision and estimation data in main experiment
if strcmp(subjectID, 'average')
    experimentNumber = 1:5;
else
    experimentNumber = 1;
end
experimentType = 'MainExperiment';
experiment = 'Original';
session = 1;
dataAll = [];
for ii = 1 : length(experimentNumber)
    dataFullPath = fullfile('Data', subjectID, experimentType, [experiment num2str(session)], ...
                                [experiment '-' num2str(experimentNumber(ii))]);          
    if iscell(dataFullPath)                       
        load(dataFullPath{1});
    else
        load(dataFullPath);
    end
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
indexAdjust = find(abs(angleTrue-angleEstimate) > 70);
angleEstimate(indexAdjust) = angleEstimate(indexAdjust) + sign(angleTrue(indexAdjust)-angleEstimate(indexAdjust))*180;    
angleDiffEst = dataAll(:, 5) - angleEstimate;
indexAdjust = abs(angleDiffEst)>90;
angleDiffEst(indexAdjust) = angleDiffEst(indexAdjust) - 180 * sign(angleDiffEst(indexAdjust));
binaryDecisionAll = dataAll(:,8);
indicatorConsistent = sign(angleDiffEst)==dataAll(:,8);
if includeIncongruentTrials == 1
    indicatorConsistent = ones(length(angleDiffEst), 1);
end
if correctIndex == 1
    indicatorCorrect = dataAll(:,6) == dataAll(:,8);
else
    indicatorCorrect = dataAll(:,6) ~= dataAll(:,8);
end
angleDiffEstimate = cell(length(stdNoiseLevel), length(angleDiff));
binaryDecisionSelf = cell(length(stdNoiseLevel), length(angleDiff));
binaryDecisionGiven = binaryDecisionSelf;

%%Extract data into angleDiffEstimate (each cell is estimate data for a condition)
for ii = 1 : length(stdNoiseLevel)
    for jj = 1 : length(angleDiff)
        indexSelect = dataAll(:,4)==stdNoiseLevel(ii) & dataAll(:,1)==angleDiff(jj); 
        indexSelectLinear = find(indexSelect & indicatorConsistent & indicatorCorrect);
        if ~isempty(indexSelectLinear)                
            indexResample = indexSelectLinear(randi(length(indexSelectLinear),[1 length(indexSelectLinear)]));
            angleDiffEstimate{ii,jj}  = angleDiffEst(indexResample);        
            binaryDecisionSelf{ii,jj} = sign(angleDiffEstimate{ii,jj});
            binaryDecisionGiven{ii,jj} = binaryDecisionAll(indexResample, :);     
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
binaryDecisionAll = dataAll(:, 6);
for ii = 1 : length(stdNoiseLevel)
    for jj = 1 : length(percentCW)
        indexSelect = dataAll(:,1) == angleDiff(jj) ...
                    & dataAll(:,4) == stdNoiseLevel(ii)...
                    & indicatorConsistent;
        if sum(indexSelect) > 0  
            indexSelectCW = dataAll(:,1) == angleDiff(jj) ...
                        & dataAll(:,4) == stdNoiseLevel(ii)...
                        & indicatorConsistent ...
                        & binaryDecisionAll==1;     
            indexSelectCW(indexSelect == 0) = [];            
            indexResample = indexSelectCW(randi(length(indexSelectCW),[1 length(indexSelectCW)]));                        
            nTrialsPerCondition(ii,jj) = nansum((indexSelect));
            percentCW(ii,jj) = 100*nansum(indexResample)/nTrialsPerCondition(ii,jj);
        end
    end
end

%% Extract motor noise in motor experiment
if strcmp(subjectID, 'average')
    experimentNumber = 1:5;
end
experimentType = 'MotorNoise';
session = 1;
stdMotorAll = NaN(1, length(experimentNumber));

for ii = 1 : length(experimentNumber)
    dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                ['MotorNoise-' num2str(experimentNumber(ii))]);
    if iscell(dataFullPath)                       
        load(dataFullPath{1});
    else
        load(dataFullPath);
    end

    % Extract data
    angleEstimate = dataTotal(:, 7);
%     nEstimate = sum(~isnan(angleEstimate));
%     params.trialOrder(nEstimate+1:end) = [];
%     angleEstimate(params.trialOrder(1:round(nEstimate/2))) = NaN;

    % Average across conditions
    angleTrue = dataTotal(:,5) - dataTotal(:,1);
    angleTrue(angleTrue>180 | angleTrue<0) = angleTrue(angleTrue>180 | angleTrue<0)...
                                                - 180 * sign(angleTrue(angleTrue>180 | angleTrue<0));
    indexAdjust = find(abs(angleTrue-angleEstimate) > 70);
    angleEstimate(indexAdjust) = angleEstimate(indexAdjust) + sign(angleTrue(indexAdjust)-angleEstimate(indexAdjust))*180;    
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