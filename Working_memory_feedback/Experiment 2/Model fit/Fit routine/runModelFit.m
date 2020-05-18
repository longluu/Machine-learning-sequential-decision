%%%%%%%%%%%%% Run the full model fitting program %%%%%%%%%%%%%
clear 
nLoops = 50;
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
% initialValueAll = [4.3425    6.2248           0.0000     19.1377    -9.6045   3.5283    2.0902    0.9927    0.6813;
%                 8.7205    8.8878           0.0000     33.0312   -21.2586   1.2076    1.8928    0.9215    0.4348;
%                 8.3099    8.7268           0.0000     13.9033   -12.6251   0.6127    2.7094    0.9468    0.5099;
%                 6.3379    8.4823           0.0000     19.3091   -12.7690   0.9699    2.6041    0.5558    0.4201;
%                 6.4379    9.9076           0.0000     32.5355   -17.6389   3.0964    1.5830    0.2067    0.4086;
%                 9.7284   15.1847           0.0000     56.3034   -42.2707   0.4378    4.0136    0.9881    0.5523;
%                 7.8510    9.8641           0.0000     22.8922   -18.0641   0.8543    3.9069    0.9646    0.4949;
%                 6.0243    8.0650           0.0000     32.4649   -19.9032   5.3989    2.4021    0.9986    0.5303];
                %[stdNoise                  lapseRate  priorRange   stdMem  stdMotor smoothPrior    pCW]
                
initialValueAll = [6.0243    8.0650           0.0000     32.4649   -19.9032   5.3989    2.4021    0.9986    0.5303];
                
binnedLH = 0;
plotFitProgress = 0;
plotEstFit = 0;
subjectID = {'ll', 'pw', 'eh', 'bh', 'ln', 'at', 'dh', 'average'}; % 'll', 'pw', 'eh', 'bh', 'ln', 'at', 'dh', 
expNumber = 1;
includeIncongruentTrials = 0; 
correctTrialOnly = 1; 
fixMotorNoise = 1;
fixLapseRate = 1;
fixSmoothPrior = 0;
fixPcw = 0;

for ll = 1 : length(subjectID)
    initialValue = initialValueAll(ll, :);
    %% Extract the data
    [binaryDecision, percentCW, nTrialsPerCondition, estimateData, angleDiff, stdMotor] = dataForFitting(subjectID(ll), includeIncongruentTrials, correctTrialOnly);
        
    save('dataAll', 'angleDiff', 'estimateData',...
             'stdMotor', 'nTrialsPerCondition', 'percentCW')
    if fixMotorNoise
        initialValue(end-2) = mean(stdMotor);
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
    fprintf(fileID,'%11s %8s %9s %18s %20s  %20s %12s  %12s %12s %12s \r\n', '//-LLHTotal', '-LLH1', '-LLH2', ...
                    'stdNoiseLevel', 'LapseRate', 'Prior range', 'MemNoise', 'MotorNoise', 'SmoothPrior', 'pCW');
    counter = 1;
    for kk = 1 : nLoops
        [tempFitParameter, tempNegLLH] = modelFitBayes_01(optimizationAlgorithm, SetStartPoint, initialValue, fileID, expNumber, plotFitProgress, modelType, '', fixMotorNoise, includeIncongruentTrials, fixLapseRate, fixSmoothPrior, fixPcw);
        fitParameterAll{kk} = tempFitParameter;
        negLLH(kk) = tempNegLLH;

        % Save the result
        fitParameter = fitParameterAll{kk};
        fprintf(fileID, ['//Iteration' '-' num2str(counter) ' \r\n']);
        fprintf(fileID,'%9.1f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f %9.4f %9.4f \r\n', negLLH(kk),...
            fitParameter(1), fitParameter(2), fitParameter(3), fitParameter(4), fitParameter(5), fitParameter(6), fitParameter(7), fitParameter(8), fitParameter(9));
        counter = counter + 1;
        
    end

    % Print the best params
    negLLH = negLLH(:);
    [negLLHSort,indSort] = sort(negLLH, 'ascend')
    fitParameterAll = fitParameterAll(:);
    fitParameter = fitParameterAll{indSort(1)};
    fprintf(fileID, ['//Best params' ' \r\n']);
    fprintf(fileID,'%2s %9.1f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f %9.4f %9.4f \r\n','//', negLLH(indSort(1)),...
        fitParameter(1), fitParameter(2), fitParameter(3), fitParameter(4), fitParameter(5), fitParameter(6), fitParameter(7), fitParameter(8), fitParameter(9));
    fclose(fileID);
end