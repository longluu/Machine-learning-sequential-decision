%%%%%%%%%%%%% Run the full model fitting program %%%%%%%%%%%%%
clear 
nResample = 81;
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
% initialValueAll = [2.6500    6.0895           0.0000     22.2852     1.6506   0.9414    2.0976  0;
%                     3.0023    9.7384           0.0000     34.4053     0.0615   0.9480    3.1069 0;
%                     4.6136   10.4165           0.0000     29.8375     0.1325   0.9940    3.8106 0;
%                     7.7094   11.9114           0.0000     55.7419     0.0083   0.2850    3.8551 0;
%                     5.1033   10.3703           0.0000     46.6421     4.7921   0.8187    3.3313 0;
%                     5.2681   10.4547           0.0000     50.3229     0.2596   0.6929    3.3234 0];
%                 %[stdNoiseInitial                lapseRateInitial  priorRangeInitial stdDelayInitial smoothFactorInitial stdMotor]

initialValueAll = [7.7094   11.9114           0.0000     55.7419     0.0083   0.2850    3.8551 0];

binnedLH = 0;
plotFitProgress = 0;
plotEstFit = 0;
subjectID = {'jp'};% 'll', 'an', 'ep', , 'kc' , 'average'
expNumber = 1;
includeIncongruentTrials = 0; 
correctTrialOnly = 1; 
fixMotorNoise = 1;
fixLapseRate = 1;
fixBoundaryCutoff = 1;

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
    fprintf(fileID,'%11s %8s %9s %18s %20s  %8s %12s %12s %12s \r\n', '//-LLHTotal', '-LLH1', '-LLH2', ...
                    'stdNoiseLevel', 'LapseRate', 'Prior', 'MemNoise', 'Smoothness', 'MotorNoise');
    counter = 1;
    for kk = 1 : nResample
        % Resample the data  
        [binaryDecision, percentCW, nTrialsPerCondition, estimateData, angleDiff, stdMotor] = dataForBootstrap(subjectID(ii), includeIncongruentTrials, correctTrialOnly);

        save('dataAll', 'angleDiff', 'estimateData',...
                 'stdMotor', 'nTrialsPerCondition', 'percentCW')
        
        
        % Fit the resample
        [tempFitParameter, tempNegLLH] = modelFitBayes_01(optimizationAlgorithm, SetStartPoint, initialValue, fileID, expNumber, plotFitProgress, modelType, '', fixMotorNoise, includeIncongruentTrials, fixLapseRate, fixBoundaryCutoff);
        fitParameterAll{kk} = tempFitParameter;
        negLLH(kk) = tempNegLLH;

        % Save the result
        fitParameter = fitParameterAll{kk};
        fprintf(fileID, ['//Iteration' '-' num2str(counter) ' \r\n']);
        fprintf(fileID,'%9.1f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f %9.4f \r\n', negLLH(kk),...
            fitParameter(1), fitParameter(2), fitParameter(3), fitParameter(4), fitParameter(5), fitParameter(6), fitParameter(7), fitParameter(8));
        counter = counter + 1;
    end

    % Print the best params
    negLLH = negLLH(:);
    [negLLHSort,indSort] = sort(negLLH, 'ascend')
    fitParameterAll = fitParameterAll(:);
    fitParameter = fitParameterAll{indSort(1)};
    fprintf(fileID, ['//Best params' ' \r\n']);
    fprintf(fileID,'%2s %9.1f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f\r\n','//', negLLH(indSort(1)),...
        fitParameter(1), fitParameter(2), fitParameter(3), fitParameter(4), fitParameter(5), fitParameter(6), fitParameter(7), fitParameter(8));
    fclose(fileID);
end