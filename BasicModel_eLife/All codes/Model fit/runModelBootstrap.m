%%%%%%%%%%%%% Run the full model fitting program %%%%%%%%%%%%%
clear 
nResample = 100;
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
initialValueAll = [4.0300    5.3505    8.1590      0.0000    41.7649   2.0957    0.0257    1.8431;
                   3.8692    6.4042   12.3882      0.0000    39.5010  14.0727    0.2322    2.9146;
                   3.7041    4.6566   10.9916      0.0000    36.4739   0.2433    0.6481    5.4045;
                   3.9226    5.7843   10.9330      0.0000    39.0534  14.2463    0.1492    5.8917;
                   3.4994    5.4317    9.5625      0.0000    39.7741   6.0913    0.2953    3.5999];
                %[stdNoiseInitial                lapseRateInitial  priorRangeInitial stdDelayInitial smoothFactorInitial stdMotor]
binnedLH = 0;
plotFitProgress = 0;
plotEstFit = 0;
subjectID = {'sy', 'vs'};
expNumber = 1;
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
        [binaryDecision, percentCW, nTrialsPerCondition, estimateData, angleDiff, stdMotor] = dataForBootstrap(subjectID(ii), expNumber, includeIncongruentTrials);
        if fixMotorNoise
            initialValue(end) = mean(stdMotor);
        end
        
        save('dataAll', 'angleDiff', 'estimateData',...
         'stdMotor', 'nTrialsPerCondition', 'percentCW')

        
        % Fit the resample
        [tempFitParameter, tempNegLLH] = modelFitBayes_New(optimizationAlgorithm, SetStartPoint, initialValue, fileID, expNumber, plotFitProgress, modelType, '', fixMotorNoise, includeIncongruentTrials, fixLapseRate);
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