%%%%%%%%%%%%% Run the full model fitting program %%%%%%%%%%%%%
clear 
nLoops = 30;
negLLHMinLoop = [];
fitParameterAll = cell(nLoops, 1);
negLLH = NaN(nLoops, 1);

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
SetStartPoint = 0;
initialValueAll = [5.1033   10.3703           0.0000     30     3   0.8187    3.3313    0.0000];
                %[stdNoise  lapseRate  priorRange stdMem smoothFactor stdMotor boundaryCutoff]
binnedLH = 0;
plotFitProgress = 0;
plotEstFit = 0;
subjectID = {'average'};
expNumber = 1;
includeIncongruentTrials = 0; 
correctTrialOnly = 1; 
fixMotorNoise = 1;
fixLapseRate = 1;
fixBoundaryCutoff = 1;

for ll = 1 : length(subjectID)
    initialValue = initialValueAll;
    %% Extract the data
    [binaryDecision, percentCW, nTrialsPerCondition, estimateData, angleDiff, stdMotor] = dataForFitting(subjectID(ll), includeIncongruentTrials, correctTrialOnly);
        
    save('dataAll', 'angleDiff', 'estimateData',...
             'stdMotor', 'nTrialsPerCondition', 'percentCW')
    if fixMotorNoise
        initialValue(end) = mean(stdMotor);
    end

    %% Start the fitting
    currentFolder = pwd;
    fileNameRoot = 'FitResult';
    fileNumber = GetNextDataFileNumber(currentFolder, '.txt');
    fileID = fopen([fileNameRoot '-' num2str(fileNumber) '.txt'],'w');
    if isempty(subjectID{ll})
        subjectID{ll} = 'Average';
    end
    subjectName = ['Subject: ' subjectID{ll}];
    expID = ['Exp: ' num2str(expNumber)];
    fprintf(fileID, '%27s %27s \r\n', subjectName, expID); 
    fprintf(fileID,'%11s %8s %9s %18s %20s  %8s %12s %12s %12s %12s \r\n', '//-LLHTotal', '-LLH1', '-LLH2', ...
                    'stdNoiseLevel', 'LapseRate', 'Prior', 'MemNoise', 'Smoothness', 'MotorNoise', 'BoundaryCutoff');
    counter = 1;
    for kk = 1 : nLoops
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