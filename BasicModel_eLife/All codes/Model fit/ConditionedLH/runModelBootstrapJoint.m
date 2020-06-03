%%%%%%%%%%%%% Run the full model fitting program %%%%%%%%%%%%%
clear 
nResample = 50;
negLLHMinLoop = [];
fitParameterAll = cell(nResample, 1);
negLLH = NaN(nResample, 1);

%% Choose the parameter for the model fit
% Model type 1: Self-consistent Bayes
%            2: Optimal Bayes
%            3: Gamma distribution
% Optimization algorithm 1: Nelder-Mead simplex
%                        2: Sequential quadratic programming 
%                        3: Simulated annneaing
%                        4: Genetic algorithm
% Set starting point for optimization 0: chooose randomly
%                                     1: set starting point
% Use binned or unbinned likelihood: 0: unbinned
%                                    1: binned
modelType = 1;
optimizationAlgorithm = 1;
SetStartPoint = 1;
initialValueAll = [4.7523    6.2508    8.9348           0          22.4177    5.3987    0.9076    2.4973];
                %[stdNoiseInitial                lapseRateInitial  priorRangeInitial stdDelayInitial smoothFactorInitial stdMotor]
binnedLH = 0;
plotFitProgress = 0;
plotEstFit = 0;
subjectID = {'average'};
expNumber = [2 3];
includeIncongruentTrials = ''; % empty if not include incongruent trials Incongruent
fixMotorNoise = 1;
fixLapseRate = 1;

for ii = 1 : length(subjectID)
    initialValue = initialValueAll(ii, :);        

    %% Start the fitting
    currentFolder = pwd;
    fileNameRoot = 'FitResult';
    fileNumber = GetNextDataFileNumber(currentFolder, '.txt');
    fileID = fopen([fileNameRoot '-' num2str(fileNumber) '.txt'],'w');
    if isempty(subjectID{ii})
        subjectID{ii} = 'Average';
    end
    subjectName = ['Subject: ' subjectID{ii}];
    expID = ['Exp: ' num2str(expNumber)];
    fprintf(fileID, '%2s %27s %27s \r\n', '//', subjectName, expID); 
    fprintf(fileID,'%11s %8s %9s %27s %20s  %8s %12s %12s %12s \r\n', '//-LLHTotal', '-LLH1', '-LLH2', ...
                    'stdNoiseLevel', 'LapseRate', 'Prior', 'MemNoise', 'Smoothness', 'MotorNoise');
    counter = 1;
    for kk = 1 : nResample
        % Resample the data  
        [binaryDecision, percentCW, nTrialsPerCondition, estimateDataExp2, angleDiff, stdMotor] = dataForBootstrap(subjectID(ii), expNumber(1), includeIncongruentTrials);
        [~, ~, ~, estimateDataExp3, ~] = dataForBootstrap(subjectID(ii), expNumber(2), includeIncongruentTrials);
        if fixMotorNoise
            initialValue(end) = mean(stdMotor);
        end

        save('dataAll', 'angleDiff', 'estimateDataExp2', 'estimateDataExp3',...
                 'stdMotor', 'nTrialsPerCondition', 'percentCW')
        
        % Fit the resample
        [tempFitParameter, tempNegLLH] = modelFitBayesJoint_New(optimizationAlgorithm, SetStartPoint, initialValue, fileID, expNumber, plotFitProgress, modelType, '', fixMotorNoise, includeIncongruentTrials, fixLapseRate);
        fitParameterAll{kk} = tempFitParameter;
        negLLH(kk) = tempNegLLH;

        % Save the result
        fitParameter = fitParameterAll{kk};
        fprintf(fileID, ['//Iteration' '-' num2str(counter) ' \r\n']);
        fprintf(fileID,'%9.1f %9.4f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f\r\n', negLLH(kk),...
            fitParameter(1), fitParameter(2), fitParameter(3), fitParameter(4), fitParameter(5), fitParameter(6), fitParameter(7), fitParameter(8));
        counter = counter + 1;
    end

    % Print the best params
    negLLH = negLLH(:);
    [negLLHSort,indSort] = sort(negLLH, 'ascend')
    fitParameterAll = fitParameterAll(:);
    fitParameter = fitParameterAll{indSort(1)};
    fprintf(fileID, ['//Best params' ' \r\n']);
    fprintf(fileID,'%2s %9.1f %9.4f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f\r\n','//', negLLH(indSort(1)),...
        fitParameter(1), fitParameter(2), fitParameter(3), fitParameter(4), fitParameter(5), fitParameter(6), fitParameter(7), fitParameter(8));
    fclose(fileID);
end