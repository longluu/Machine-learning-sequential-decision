%%%%%%%%%%%%% Compute the log likelihood of self-consistent,  %%%%%%%%%%%%%
%%%%%%%%%%%%% standard Bayes and oracle model %%%%%%%%%%%%%
function computeModelLLH
numBin = 31;
subjectID = {'ll', 'sy', 'cz', 'vs', 'as',  'average', 'll', 'xfl', 'aj', 'zw', 'skb', 'average'};
expNum = [1 1 1 1 1 1 2 2 2 2 2 2];
paramsSB1 = [1.7169    4.2730    7.3499           0.0000     26.1313     2.0940   1.0000    2.6327;
            4.2738    5.4032    8.4126           0.0000     49.9981     0.0226   0.1351    1.8431;
            7.5041   10.5469   17.4636           0.0000     40.9869     1.3610   0.6663    2.9146;
            4.2756    5.6413   10.9465           0.0000     32.4254     0.0094   0.8927    5.4045;
            7.1588    8.6167   14.1209           0.0000     49.8100     1.6407   0.6285    5.8917;
            5.7866    7.3978   11.6602           0.0000     41.3896     0.0105   0.9592    3.5999;
            
            1.7169    4.2730    7.3499           0.0000     22.8378     2.0940   1.0000    2.6327;
            7.2433    8.6020   12.6927           0.0000     23.2140     0.0080   0.7465    2.9350;
            6.9466    7.6698    8.7704           0.0000     17.1810     1.1578   0.0181    2.9916;
            4.7553    5.2620    5.9237           0.0000     16.4109     2.2921   0.0118    2.5348;            
            9.9503   10.7767   13.5594           0.0000     31.9482     0.7754   0.6916    2.0456;
            6.8628    7.9477   10.2439           0.0000     23.0751     0.4136   0.9996    2.4973];

paramsSC1 = [2.4953    4.7444    6.5251           0.0000     26.4013      0.0071   0.9999    2.6327;
             4.2087    5.3793    7.6340           0.0000     41.9850     0.0055   0.8409    1.8431;
             3.8788    6.3992   12.3802          0           39.8999   14.1127    0.3329    2.9146;
             3.7794    4.7956   10.2831           0.0000     29.6435     0.0100   0.2696    5.4045;
             3.9059    5.7737   10.9297          0           38.8127   14.3210    0.2936    5.8917;
             3.5066    5.4412    9.5715          0           38.4841    6.1243    0.3334    3.5999; 
             
             2.4953    4.7444    6.5251           0.0000     22.0477     0.0071   0.9999    2.6327;             
             6.3327    8.4221   14.6399           0          22.2658    6.3080    1.0000    2.9350;
             4.5991    6.0397    7.8831           0          15.9917    1.4014    0.0120    2.9916;
             4.5627    5.4642    6.7588           0          17.4179    4.3354    0.0187    2.5348;
             4.2918    4.5450    8.4029           0          32.6473   14.7274    0.7750    2.0456;
             4.5878    6.1759    9.0336           0          22.8017    5.9303    1.0000    2.4973];

paramsSC2 = [2.5499    4.7059    6.8538         0   25.7103    0.0050    0.9979    2.6327;
             4.3642    5.3680    7.7577           0.0000     41.9866     0.0050   0.6964    1.8431;
             8.4769   10.8465   15.4944         0   37.2157    0.0050    0.7944    2.9146;
             3.6885    4.6437   10.1372         0   31.8861    0.0050    0.9862    5.4045;
             6.3837    8.6987   13.2335           0.0000     41.9764     0.0050   0.8220    5.8917;
             5.3173    6.8932    9.8580         0   49.8339    0.0050    0.5805    3.5999;
             
             2.5499    4.7059    6.8538         0   22.1601    0.0050    0.9979    2.6327;
             7.8865    9.4599   16.3686         0   20.9203    0.0050    0.9998    2.9350;
             5.0086    6.4225    7.8574         0   14.9438    0.0050    0.8581    2.9916;
             5.8036    6.5628    7.9510         0   16.7294    0.0050    0.5069    2.5348;
             9.5431   11.0343   13.3800         0   38.8442    0.0050    0.8118    2.0456;
             6.1745    7.4755    9.4670         0   22.4310    0.0050    0.9999    2.4973];

% subjectID = {'ll', 'xfl', 'aj', 'zw', 'skb', 'average', 'll', 'xfl', 'aj', 'zw', 'skb', 'average'};
% expNum = [2 2 2 2 2 2 3 3 3 3 3 3];
% paramsSB1 = [1.7169    4.2730    7.3499           0.0000     22.8378     2.0940   1.0000    2.6327;
%             7.2433    8.6020   12.6927           0.0000     23.2140     0.0080   0.7465    2.9350;
%             6.9466    7.6698    8.7704           0.0000     17.1810     1.1578   0.0181    2.9916;
%             4.7553    5.2620    5.9237           0.0000     16.4109     2.2921   0.0118    2.5348;            
%             9.9503   10.7767   13.5594           0.0000     31.9482     0.7754   0.6916    2.0456;
%             6.8628    7.9477   10.2439           0.0000     23.0751     0.4136   0.9996    2.4973;
%             
%             1.7169    4.2730    7.3499           0.0000     22.8378     2.0940   1.0000    2.6327;
%             7.2433    8.6020   12.6927           0.0000     23.2140     0.0080   0.7465    2.9350;
%             6.9466    7.6698    8.7704           0.0000     17.1810     1.1578   0.0181    2.9916;
%             4.7553    5.2620    5.9237           0.0000     16.4109     2.2921   0.0118    2.5348;            
%             9.9503   10.7767   13.5594           0.0000     31.9482     0.7754   0.6916    2.0456;
%             6.8628    7.9477   10.2439           0.0000     23.0751     0.4136   0.9996    2.4973];
% 
% paramsSC1 = [2.2692    4.5725    6.6427           0.0000     26.8237     1.5189   1.0000    2.6327;             
%              6.8748    8.6631   14.8639           0.0000     24.5196     6.6582   0.9644    2.9350;
%              4.6465    5.8595    7.9764           0.0000     15.6241     0.2107   0.0107    2.9916;
%              4.3887    5.7847    7.3365           0.0000     18.2779     6.8338   0.0224    2.5348;
%              4.3871    4.2896    8.9027           0.0000     32.0442    15.7666   0.9140    2.0456;
%              4.6472    5.9963    8.8842           0.0000     24.2385     6.2953   1.0000    2.4973;
%              
%              2.4187    4.6013    6.6852           0.0000     32.5288     1.1708   0.9992    2.6327;             
%              6.8748    8.6631   14.8639           0.0000     24.5196     6.6582   0.9644    2.9350;
%              4.6465    5.8595    7.9764           0.0000     15.6241     0.2107   0.0107    2.9916;
%              4.3887    5.7847    7.3365           0.0000     18.2779     6.8338   0.0224    2.5348;
%              4.3871    4.2896    8.9027           0.0000     32.0442    15.7666   0.9140    2.0456;
%              4.6472    5.9963    8.8842           0.0000     24.2385     6.2953   1.0000    2.4973;];
% 
% paramsSC2 = [2.5499    4.7059    6.8538         0   22.1601    0.0050    0.9979    2.6327;
%              7.8865    9.4599   16.3686         0   20.9203    0.0050    0.9998    2.9350;
%              5.0086    6.4225    7.8574         0   14.9438    0.0050    0.8581    2.9916;
%              5.8036    6.5628    7.9510         0   16.7294    0.0050    0.5069    2.5348;
%              9.5431   11.0343   13.3800         0   38.8442    0.0050    0.8118    2.0456;
%              6.1745    7.4755    9.4670         0   22.4310    0.0050    0.9999    2.4973;
%              
%              2.5499    4.7059    6.8538         0   22.1601    0.0050    0.9979    2.6327;
%              7.8865    9.4599   16.3686         0   20.9203    0.0050    0.9998    2.9350;
%              5.0086    6.4225    7.8574         0   14.9438    0.0050    0.8581    2.9916;
%              5.8036    6.5628    7.9510         0   16.7294    0.0050    0.5069    2.5348;
%              9.5431   11.0343   13.3800         0   38.8442    0.0050    0.8118    2.0456;
%              6.1745    7.4755    9.4670         0   22.4310    0.0050    0.9999    2.4973];         

logLH_Data = NaN(3, size(paramsSB1, 1));
logLH_Chance = NaN(3, size(paramsSB1, 1));
logLH_StandardBayes = NaN(3, size(paramsSB1, 1));
logLH_SelfConsistentMemory = NaN(3, size(paramsSB1, 1));
logLH_SelfConsistentNoMemory = NaN(3, size(paramsSB1, 1));

for ii = 1 : size(paramsSB1, 1)
    [logLH_Data(:, ii), logLH_StandardBayes(:, ii), logLH_SelfConsistentMemory(:, ii), logLH_SelfConsistentNoMemory(:, ii), logLH_Chance(:, ii)] = computeLLH(subjectID(ii), expNum(ii), paramsSB1(ii, :),...
                                                        paramsSC1(ii, :), paramsSC2(ii, :), numBin);
end


%% Plot the LLH
logLH_ChanceNorm = (logLH_Chance - logLH_StandardBayes) ./ (logLH_Data - logLH_StandardBayes);
logLH_SelfConsistentNorm1 = (logLH_SelfConsistentNoMemory - logLH_StandardBayes) ./ (logLH_Data - logLH_StandardBayes);
logLH_SelfConsistentNorm2 = (logLH_SelfConsistentMemory - logLH_StandardBayes) ./ (logLH_Data - logLH_StandardBayes);

% Plot the self-consistent model with and without memory noise
logLH_ChanceAve = mean(logLH_ChanceNorm, 1);
logLH_SelfConsistentAve1 = mean(logLH_SelfConsistentNorm1, 1);
logLH_SelfConsistentAve2 = mean(logLH_SelfConsistentNorm2, 1);
logLH_ModelAll = [logLH_SelfConsistentAve1; logLH_SelfConsistentAve2];

% hFig = figure;
% fontSize = 23;
% colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
% colorIndex = NaN(length(colorName), 3);
% for ii = 1 : length(colorName)
%     colorIndex(ii, :) = rgb(colorName{ii});
% end
% 
% hAx1 = gca;
% [~, hPanel] = errorbar_groups(logLH_ModelAll, zeros(size(logLH_ModelAll)), ...
%                                     zeros(size(logLH_ModelAll)), ...
%                             'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx1);
% set(hPanel, 'EdgeColor', 'none')
% set(gca, 'FontSize', fontSize)
% ylim([-0.1 1])
% ylabel('Normalized log likelihood')
% title('Experiment 1')


hFig = figure;
fontSize = 23;
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
hAx1 = subplot(1, 2, 1);
[~, hPanel] = errorbar_groups(logLH_ModelAll(:, 1:6), zeros(size(logLH_ModelAll(:, 1:6))), ...
                                    zeros(size(logLH_ModelAll(:, 1:6))), ...
                            'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx1);
set(hPanel, 'EdgeColor', 'none')
set(gca, 'FontSize', fontSize)
ylim([0 1])
ylabel('Normalized log likelihood')
title('Experiment 1')

hAx2 = subplot(1, 2, 2);
[~, hPanel] = errorbar_groups(logLH_ModelAll(:, 7:end), zeros(size(logLH_ModelAll(:, 7:end))), ...
                                    zeros(size(logLH_ModelAll(:, 7:end))), ...
                            'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx2);
set(hPanel, 'EdgeColor', 'none')
set(gca, 'FontSize', fontSize)
ylim([min(logLH_ModelAll(:))-0.1 1])
ylabel('Normalized log likelihood')
title('Experiment 2')

hFig = figure;
hAx1 = subplot(1, 2, 1);
[~, hPanel] = errorbar_groups(logLH_ChanceAve(:, 1:6), zeros(size(logLH_ChanceAve(:, 1:6))), ...
                                    zeros(size(logLH_ChanceAve(:, 1:6))), ...
                            'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx1);
set(hPanel, 'EdgeColor', 'none')
set(gca, 'FontSize', fontSize)
ylim([min(logLH_ChanceAve)-0.1 1])
text(1:length(logLH_ChanceAve(:, 1:6)), logLH_ChanceAve(:, 1:6), num2str(round(logLH_ChanceAve(:, 1:6))'),'vert','top','horiz','center'); 
ylabel('Normalized log likelihood')
title('Experiment 1')

hAx2 = subplot(1, 2, 2);
[~, hPanel] = errorbar_groups(logLH_ChanceAve(:, 7:end), zeros(size(logLH_ChanceAve(:, 7:end))), ...
                                    zeros(size(logLH_ChanceAve(:, 7:end))), ...
                            'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx2);
set(hPanel, 'EdgeColor', 'none')
set(gca, 'FontSize', fontSize)
ylim([min(logLH_ChanceAve)-0.1 1])
text(1:length(logLH_ChanceAve(:, 7:end)), logLH_ChanceAve(:, 7:end), num2str(round(logLH_ChanceAve(:, 7:end))'),'vert','top','horiz','center'); 
ylabel('Normalized log likelihood')
title('Experiment 2')

% % Plot the self-consistent model with memory noise
% hFig = figure;
% hAx1 = subplot(1, 2, 1);
% [~, hPanel] = errorbar_groups(logLH_SelfConsistentNorm1(:, 1:6), zeros(size(logLH_SelfConsistentNorm1(:, 1:6))), ...
%                                     zeros(size(logLH_SelfConsistentNorm1(:, 1:6))), ...
%                             'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx1);
% set(hPanel, 'EdgeColor', 'none')
% set(gca, 'FontSize', fontSize)
% ylim([-0.1 1])
% ylabel('Normalized log likelihood')
% title('Experiment 1')
% 
% hAx2 = subplot(1, 2, 2);
% [~, hPanel] = errorbar_groups(logLH_SelfConsistentNorm1(:, 7:end), zeros(size(logLH_SelfConsistentNorm1(:, 7:end))), ...
%                                     zeros(size(logLH_SelfConsistentNorm1(:, 7:end))), ...
%                             'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig, 'AxID', hAx2);
% set(hPanel, 'EdgeColor', 'none')
% set(gca, 'FontSize', fontSize)
% ylim([-0.1 1])
% ylabel('Normalized log likelihood')
% title('Experiment 2')



end
function [logLH_Data, logLH_StandardBayes, logLH_SelfConsistentMemory, logLH_SelfConsistentNoMemory, logLH_Chance] = computeLLH(subjectID, expNumber, paramsStandardBayes, paramsSelfConsistent1, paramsSelfConsistent2, numBin)
includeIncongruentTrials = ''; % empty if not include incongruent trials Incongruent

if expNumber == 1
    %##################################### Experiment 1 #####################################
    %% Extract the data
    [~, percentCW, nTrialsPerCondition, estimateData, angleDiff, ~] = dataForFitting(subjectID, expNumber, includeIncongruentTrials);

    %% Oracle model
    % Compute LLH of discrimination
    pCW = percentCW/100;
    numberCW = round(nTrialsPerCondition .* pCW);
    numberCCW = nTrialsPerCondition-numberCW;
    logLHDiscriminate = nansum(numberCW .* log(pCW) + ...
                        numberCCW .* log(1-pCW), 2);

    % Compute LLH of estimates
    logLHEstimate = zeros(3, 1);
    binCenter = linspace(-40, 40, numBin);
    deltaBin = diff(binCenter(1:2));
    dstep = 0.1;
    rangeth = [-60 60];
    th = rangeth(1):dstep:rangeth(2);
    binEdge = [min(th) binCenter(1:end-1)+deltaBin max(th)];
    binCountAll = cell(size(estimateData));
    for ii = 1 : size(estimateData, 1)
        for jj = 1 : size(estimateData, 2)
            binCount = hist(estimateData{ii, jj}, binCenter);
            pEmpirical = binCount / sum(binCount);
            tempLogLH = binCount .* log(pEmpirical);
            tempLogLH = nansum(tempLogLH);
            if ~isnan(tempLogLH)
                logLHEstimate(ii) = logLHEstimate(ii) +  tempLogLH;
            end
            binCountAll{ii, jj} = binCount;
        end
    end
    logLH_Data = logLHDiscriminate + logLHEstimate;

    %% Chance model
    % Compute LLH of discrimination
    pChanceModel = 0.5 * ones(size(percentCW));
    pCW = percentCW/100;
    numberCW = round(nTrialsPerCondition .* pCW);
    numberCCW = nTrialsPerCondition-numberCW;
    logLHDiscriminate = nansum(numberCW .* log(pChanceModel) + ...
                        numberCCW .* log(1-pChanceModel), 2);

    % Compute LLH of estimates
    logLHEstimate = zeros(3, 1);
    thFullCircle = -90:0.01:89.99;
    countEmpirical = hist(thFullCircle, binCenter);
    pEmpirical = countEmpirical / sum(countEmpirical);
    pEmpirical(pEmpirical==0) = 10^(-10);
    for ii = 1 : size(estimateData, 1)
        for jj = 1 : size(estimateData, 2)
            binCount = hist(estimateData{ii, jj}, binCenter);
            tempLogLH = binCount .* log(pEmpirical);
            tempLogLH = nansum(tempLogLH);
            if ~isnan(tempLogLH)
                logLHEstimate(ii) = logLHEstimate(ii) +  tempLogLH;
            end
        end
    end
    logLH_Chance = logLHDiscriminate + logLHEstimate;
    
    %% Standard Bayes
    % Fit paramters
    paramsAll = paramsStandardBayes;
    lapseRate = paramsAll(4);
    stdSensory = paramsAll(1:3);
    stdMemory = paramsAll(6);
    stdMotor = paramsAll(8);
    pC = [0.5, 0.5]'; % [cw ccw]
    priorRange = paramsAll(5);
    smoothFactor = paramsAll(7);

    % Run the model
    [pCGthetaAll, estimateModel] = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC,...
                                expNumber, angleDiff, stdMotor, 2, lapseRate, 'include');

    % LLH of discrimination
    logLHDiscriminate = nansum(numberCW .* log(pCGthetaAll) + ...
                        numberCCW .* log(1-pCGthetaAll), 2);

    % LLH of estimates                        
    logLH_Estimate = zeros(3, 1);
    for kk = 1 : length(stdSensory)                        
        rangeCollapse = round(length(angleDiff)/2);
        tempEstimateModelX = estimateModel.Xval{kk};
        pthhANDth_incorrect = estimateModel.Yval{kk};
        for jj = rangeCollapse : length(angleDiff)
            tempEstimateModelY = pthhANDth_incorrect(:, jj); 
            tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
            pBin = NaN(1, length(binCenter));
            for ii = 1 : length(binCenter)
                pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
            end
            binCount = binCountAll{kk,jj-rangeCollapse+1};
            pBin(pBin == 0) = NaN;
            logLH_Estimate(kk) = logLH_Estimate(kk) + nansum(binCount .* log(pBin));             
        end
    end

    logLH_StandardBayes = logLHDiscriminate + logLH_Estimate;

    %% Self-consistent Bayes (with memory noise)
    % Fit paramters
    paramsAll = paramsSelfConsistent1;
    lapseRate = paramsAll(4);
    stdSensory = paramsAll(1:3);
    stdMemory = paramsAll(6);
    stdMotor = paramsAll(8);
    pC = [0.5, 0.5]'; % [cw ccw]
    priorRange = paramsAll(5);
    smoothFactor = paramsAll(7);

    % Run the model
    [pCGthetaAll, estimateModel] = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC,...
                                expNumber, angleDiff, stdMotor, 1, lapseRate, includeIncongruentTrials);

    % LLH of discrimination
    logLHDiscriminate = nansum(numberCW .* log(pCGthetaAll) + ...
                        numberCCW .* log(1-pCGthetaAll), 2);

    % LLH of estimates                        
    logLH_Estimate = zeros(3, 1);
    for kk = 1 : length(stdSensory)                        
        rangeCollapse = round(length(angleDiff)/2);
        tempEstimateModelX = estimateModel.Xval{kk};
        pthhANDth_incorrect = estimateModel.Yval{kk};
        for jj = rangeCollapse : length(angleDiff)
            tempEstimateModelY = pthhANDth_incorrect(:, jj); 
            tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
            pBin = NaN(1, length(binCenter));
            for ii = 1 : length(binCenter)
                pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
            end
            binCount = binCountAll{kk,jj-rangeCollapse+1};
            pBin(pBin == 0) = NaN;
            logLH_Estimate(kk) = logLH_Estimate(kk) + nansum(binCount .* log(pBin));             
        end
    end

    logLH_SelfConsistentMemory = logLHDiscriminate + logLH_Estimate;
    
    %% Self-consistent Bayes (no memory noise)
    % Fit paramters
    paramsAll = paramsSelfConsistent2;
    lapseRate = paramsAll(4);
    stdSensory = paramsAll(1:3);
    stdMemory = paramsAll(6);
    stdMotor = paramsAll(8);
    pC = [0.5, 0.5]'; % [cw ccw]
    priorRange = paramsAll(5);
    smoothFactor = paramsAll(7);

    % Run the model
    [pCGthetaAll, estimateModel] = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC,...
                                expNumber, angleDiff, stdMotor, 1, lapseRate, includeIncongruentTrials);

    % LLH of discrimination
    logLHDiscriminate = nansum(numberCW .* log(pCGthetaAll) + ...
                        numberCCW .* log(1-pCGthetaAll), 2);

    % LLH of estimates                        
    logLH_Estimate = zeros(3, 1);
    for kk = 1 : length(stdSensory)                        
        rangeCollapse = round(length(angleDiff)/2);
        tempEstimateModelX = estimateModel.Xval{kk};
        pthhANDth_incorrect = estimateModel.Yval{kk};
        for jj = rangeCollapse : length(angleDiff)
            tempEstimateModelY = pthhANDth_incorrect(:, jj); 
            tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
            pBin = NaN(1, length(binCenter));
            for ii = 1 : length(binCenter)
                pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdge(ii) & tempEstimateModelX < binEdge(ii+1)));
            end
            binCount = binCountAll{kk,jj-rangeCollapse+1};
            pBin(pBin == 0) = NaN;
            logLH_Estimate(kk) = logLH_Estimate(kk) + nansum(binCount .* log(pBin));             
        end
    end

    logLH_SelfConsistentNoMemory = logLHDiscriminate + logLH_Estimate;    
elseif expNumber == 2
    %##################################### Experiment 2 #####################################

    %% Extract the data
    [~, percentCW, nTrialsPerCondition, estimateDataExp2, angleDiff, ~] = dataForFitting(subjectID, 2, includeIncongruentTrials);

    %% Oracle model
    % Compute LLH of discrimination
    pCW = percentCW/100;
    numberCW = round(nTrialsPerCondition .* pCW);
    numberCCW = nTrialsPerCondition-numberCW;
    logLHDiscriminate = nansum(numberCW .* log(pCW) + ...
                        numberCCW .* log(1-pCW), 2);

    % Compute LLH of estimates
    logLHEstimateExp2 = zeros(3, 1);
    binCenter = linspace(-40, 40, numBin);
    deltaBin = diff(binCenter(1:2));
    rangeth = [-60 60];
    binEdgeExp2 = [min(rangeth) binCenter(1:end-1)+deltaBin max(rangeth)];
    binCountAllExp2 = cell(size(estimateDataExp2));
    for ii = 1 : size(estimateDataExp2, 1)
        for jj = 1 : size(estimateDataExp2, 2)
            binCount = hist(estimateDataExp2{ii, jj}, binCenter);
            pEmpirical = binCount / sum(binCount);
            tempLogLH = binCount .* log(pEmpirical);
            tempLogLH = nansum(tempLogLH);
            if ~isnan(tempLogLH)
                logLHEstimateExp2(ii) = logLHEstimateExp2(ii) +  tempLogLH;
            end
            binCountAllExp2{ii, jj} = binCount;
        end
    end


    logLH_Data = logLHDiscriminate + logLHEstimateExp2;

    %% Chance model
    % Compute LLH of discrimination
    pChanceModel = 0.5 * ones(size(percentCW));
    pCW = percentCW/100;
    numberCW = round(nTrialsPerCondition .* pCW);
    numberCCW = nTrialsPerCondition-numberCW;
    logLHDiscriminate = nansum(numberCW .* log(pChanceModel) + ...
                        numberCCW .* log(1-pChanceModel), 2);

    % Compute LLH of estimates
    logLHEstimateExp2 = zeros(3, 1);
    thFullCircle = -90:0.01:89.99;
    countEmpirical = hist(thFullCircle, binCenter);
    pEmpirical = countEmpirical / sum(countEmpirical);
    pEmpirical(pEmpirical==0) = 10^(-10);
    
    for ii = 1 : size(estimateDataExp2, 1)
        for jj = 1 : size(estimateDataExp2, 2)
            binCount = hist(estimateDataExp2{ii, jj}, binCenter);
            tempLogLH = binCount .* log(pEmpirical);
            tempLogLH = nansum(tempLogLH);
            if ~isnan(tempLogLH)
                logLHEstimateExp2(ii) = logLHEstimateExp2(ii) +  tempLogLH;
            end
        end
    end
        
    logLH_Chance = logLHDiscriminate + logLHEstimateExp2;
    
    %% Standard Bayes
    % Fit paramters
    paramsAll = paramsStandardBayes;
    lapseRate = paramsAll(4);
    stdSensory = paramsAll(1:3);
    stdMemory = paramsAll(6);
    stdMotor = paramsAll(8);
    pC = [0.5, 0.5]'; % [cw ccw]
    priorRange = paramsAll(5);
    smoothFactor = paramsAll(7);

    % Run the model
    [pCGthetaAll, estimateModelExp2] = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC,...
                                2, angleDiff, stdMotor, 2, lapseRate, 'include');

    % LLH of discrimination
    logLHDiscriminate = nansum(numberCW .* log(pCGthetaAll) + ...
                        numberCCW .* log(1-pCGthetaAll), 2);

    % LLH of estimates                        
    logLH_EstimateExp2 = zeros(3, 1);
    for kk = 1 : length(stdSensory)                        
        rangeCollapse = round(length(angleDiff)/2);
        tempEstimateModelX = estimateModelExp2.Xval{kk};
        pthhANDth_incorrect = estimateModelExp2.Yval{kk};
        for jj = rangeCollapse : length(angleDiff)
            tempEstimateModelY = pthhANDth_incorrect(:, jj); 
            tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
            pBin = NaN(1, length(binCenter));
            for ii = 1 : length(binCenter)
                pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdgeExp2(ii) & tempEstimateModelX < binEdgeExp2(ii+1)));
            end
            binCount = binCountAllExp2{kk,jj-rangeCollapse+1};
            pBin(pBin == 0) = NaN;
            logLH_EstimateExp2(kk) = logLH_EstimateExp2(kk) + nansum(binCount .* log(pBin));             
        end
    end

    logLH_StandardBayes = logLHDiscriminate + logLH_EstimateExp2;

    %% Self-consistent Bayes (with memory noise)
    % Fit paramters
    paramsAll = paramsSelfConsistent1;
    lapseRate = paramsAll(4);
    stdSensory = paramsAll(1:3);
    stdMemory = paramsAll(6);
    stdMotor = paramsAll(8);
    pC = [0.5, 0.5]'; % [cw ccw]
    priorRange = paramsAll(5);
    smoothFactor = paramsAll(7);

    % Run the model
    [pCGthetaAll, estimateModelExp2] = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC,...
                                2, angleDiff, stdMotor, 1, lapseRate, includeIncongruentTrials);

    % LLH of discrimination
    logLHDiscriminate = nansum(numberCW .* log(pCGthetaAll) + ...
                        numberCCW .* log(1-pCGthetaAll), 2);

    % LLH of estimates                        
    logLH_EstimateExp2 = zeros(3, 1);
    for kk = 1 : length(stdSensory)                        
        rangeCollapse = round(length(angleDiff)/2);
        tempEstimateModelX = estimateModelExp2.Xval{kk};
        pthhANDth_incorrect = estimateModelExp2.Yval{kk};
        for jj = rangeCollapse : length(angleDiff)
            tempEstimateModelY = pthhANDth_incorrect(:, jj); 
            tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
            pBin = NaN(1, length(binCenter));
            for ii = 1 : length(binCenter)
                pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdgeExp2(ii) & tempEstimateModelX < binEdgeExp2(ii+1)));
            end
            binCount = binCountAllExp2{kk,jj-rangeCollapse+1};
            pBin(pBin == 0) = NaN;
            logLH_EstimateExp2(kk) = logLH_EstimateExp2(kk) + nansum(binCount .* log(pBin));             
        end
    end

    logLH_SelfConsistentMemory = logLHDiscriminate + logLH_EstimateExp2;
    
    %% Self-consistent Bayes (no memory noise)
    % Fit paramters
    paramsAll = paramsSelfConsistent2;
    lapseRate = paramsAll(4);
    stdSensory = paramsAll(1:3);
    stdMemory = paramsAll(6);
    stdMotor = paramsAll(8);
    pC = [0.5, 0.5]'; % [cw ccw]
    priorRange = paramsAll(5);
    smoothFactor = paramsAll(7);

    % Run the model
    [pCGthetaAll, estimateModelExp2] = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC,...
                                2, angleDiff, stdMotor, 1, lapseRate, includeIncongruentTrials);
                            
    % LLH of discrimination
    logLHDiscriminate = nansum(numberCW .* log(pCGthetaAll) + ...
                        numberCCW .* log(1-pCGthetaAll), 2);

    % LLH of estimates                        
    logLH_EstimateExp2 = zeros(3, 1);
    for kk = 1 : length(stdSensory)                        
        rangeCollapse = round(length(angleDiff)/2);
        tempEstimateModelX = estimateModelExp2.Xval{kk};
        pthhANDth_incorrect = estimateModelExp2.Yval{kk};
        for jj = rangeCollapse : length(angleDiff)
            tempEstimateModelY = pthhANDth_incorrect(:, jj); 
            tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
            pBin = NaN(1, length(binCenter));
            for ii = 1 : length(binCenter)
                pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdgeExp2(ii) & tempEstimateModelX < binEdgeExp2(ii+1)));
            end
            binCount = binCountAllExp2{kk,jj-rangeCollapse+1};
            pBin(pBin == 0) = NaN;
            logLH_EstimateExp2(kk) = logLH_EstimateExp2(kk) + nansum(binCount .* log(pBin));             
        end
    end

    logLH_SelfConsistentNoMemory = logLHDiscriminate + logLH_EstimateExp2;   
else
    %##################################### Experiment 3 #####################################

    %% Extract the data
    [~, ~, ~, estimateDataExp3, angleDiff] = dataForFitting(subjectID, 3, includeIncongruentTrials);

    %% Oracle model
    % Compute LLH of estimates
    logLHEstimateExp3 = zeros(3, 1);
    binCenter = linspace(-42, 42, numBin);
    deltaBin = diff(binCenter(1:2));
    rangeth = [-42 42];
    binEdgeExp3 = [min(rangeth) binCenter(1:end-1)+deltaBin max(rangeth)];
    binCountAllExp3 = cell(size(estimateDataExp3));
    for ii = 1 : size(estimateDataExp3, 1)
        for jj = 1 : size(estimateDataExp3, 2)
            binCount = hist(estimateDataExp3{ii, jj}, binCenter);
            pEmpirical = binCount / sum(binCount);
            tempLogLH = binCount .* log(pEmpirical);
            tempLogLH = nansum(tempLogLH);
            if ~isnan(tempLogLH)
                logLHEstimateExp3(ii) = logLHEstimateExp3(ii) +  tempLogLH;
            end
            binCountAllExp3{ii, jj} = binCount;
        end
    end

    logLH_Data = logLHEstimateExp3;

    %% Chance model
    % Compute LLH of estimates    
    logLHEstimateExp3 = zeros(3, 1);
    thFullCircle = -90:0.01:89.99;
    countEmpirical = hist(thFullCircle, binCenter);
    pEmpirical = countEmpirical / sum(countEmpirical);
    pEmpirical(pEmpirical==0) = 10^(-10);
    
    for ii = 1 : size(estimateDataExp3, 1)
        for jj = 1 : size(estimateDataExp3, 2)
            binCount = hist(estimateDataExp3{ii, jj}, binCenter);
            tempLogLH = binCount .* log(pEmpirical);
            tempLogLH = nansum(tempLogLH);
            if ~isnan(tempLogLH)
                logLHEstimateExp3(ii) = logLHEstimateExp3(ii) +  tempLogLH;
            end
        end
    end
    
    
    logLH_Chance = logLHEstimateExp3;
    
    %% Standard Bayes
    % Fit paramters
    paramsAll = paramsStandardBayes;
    lapseRate = paramsAll(4);
    stdSensory = paramsAll(1:3);
    stdMemory = paramsAll(6);
    stdMotor = paramsAll(8);
    pC = [0.5, 0.5]'; % [cw ccw]
    priorRange = paramsAll(5);
    smoothFactor = paramsAll(7);

    % Run the model
    [~, estimateModelExp3] = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC,...
                                3, angleDiff, stdMotor, 2, lapseRate, 'include');

    % LLH of estimates                        
    logLH_EstimateExp3 = zeros(3, 1);
    for kk = 1 : length(stdSensory)                        
        rangeCollapse = round(length(angleDiff)/2);
        tempEstimateModelX = estimateModelExp3.Xval{kk};
        pthhANDth_incorrect = estimateModelExp3.Yval{kk};
        for jj = rangeCollapse : length(angleDiff)
            tempEstimateModelY = pthhANDth_incorrect(:, jj); 
            tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
            pBin = NaN(1, length(binCenter));
            for ii = 1 : length(binCenter)
                pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdgeExp3(ii) & tempEstimateModelX < binEdgeExp3(ii+1)));
            end
            binCount = binCountAllExp3{kk,jj-rangeCollapse+1};
            pBin(pBin == 0) = NaN;
            logLH_EstimateExp3(kk) = logLH_EstimateExp3(kk) + nansum(binCount .* log(pBin));             
        end
    end

    logLH_StandardBayes = logLH_EstimateExp3;

    %% Self-consistent Bayes (with memory noise)
    % Fit paramters
    paramsAll = paramsSelfConsistent1;
    lapseRate = paramsAll(4);
    stdSensory = paramsAll(1:3);
    stdMemory = paramsAll(6);
    stdMotor = paramsAll(8);
    pC = [0.5, 0.5]'; % [cw ccw]
    priorRange = paramsAll(5);
    smoothFactor = paramsAll(7);

    % Run the model
    [~, estimateModelExp3] = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC,...
                                3, angleDiff, stdMotor, 1, lapseRate, includeIncongruentTrials);

    % LLH of estimates                        
    logLH_EstimateExp3 = zeros(3, 1);
    for kk = 1 : length(stdSensory)                        
        rangeCollapse = round(length(angleDiff)/2);
        tempEstimateModelX = estimateModelExp3.Xval{kk};
        pthhANDth_incorrect = estimateModelExp3.Yval{kk};
        for jj = rangeCollapse : length(angleDiff)
            tempEstimateModelY = pthhANDth_incorrect(:, jj); 
            tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
            pBin = NaN(1, length(binCenter));
            for ii = 1 : length(binCenter)
                pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdgeExp3(ii) & tempEstimateModelX < binEdgeExp3(ii+1)));
            end
            binCount = binCountAllExp3{kk,jj-rangeCollapse+1};
            pBin(pBin == 0) = NaN;
            logLH_EstimateExp3(kk) = logLH_EstimateExp3(kk) + nansum(binCount .* log(pBin));             
        end
    end

    logLH_SelfConsistentMemory = logLH_EstimateExp3;
    
    %% Self-consistent Bayes (no memory noise)
    % Fit paramters
    paramsAll = paramsSelfConsistent2;
    lapseRate = paramsAll(4);
    stdSensory = paramsAll(1:3);
    stdMemory = paramsAll(6);
    stdMotor = paramsAll(8);
    pC = [0.5, 0.5]'; % [cw ccw]
    priorRange = paramsAll(5);
    smoothFactor = paramsAll(7);

    % Run the model
    [~, estimateModelExp3] = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC,...
                                3, angleDiff, stdMotor, 1, lapseRate, includeIncongruentTrials);

    % LLH of estimates                        
    logLH_EstimateExp3 = zeros(3, 1);
    for kk = 1 : length(stdSensory)                        
        rangeCollapse = round(length(angleDiff)/2);
        tempEstimateModelX = estimateModelExp3.Xval{kk};
        pthhANDth_incorrect = estimateModelExp3.Yval{kk};
        for jj = rangeCollapse : length(angleDiff)
            tempEstimateModelY = pthhANDth_incorrect(:, jj); 
            tempEstimateModelY = tempEstimateModelY / sum(tempEstimateModelY);            
            pBin = NaN(1, length(binCenter));
            for ii = 1 : length(binCenter)
                pBin(ii) = sum(tempEstimateModelY(tempEstimateModelX >= binEdgeExp3(ii) & tempEstimateModelX < binEdgeExp3(ii+1)));
            end
            binCount = binCountAllExp3{kk,jj-rangeCollapse+1};
            pBin(pBin == 0) = NaN;
            logLH_EstimateExp3(kk) = logLH_EstimateExp3(kk) + nansum(binCount .* log(pBin));             
        end
    end

    logLH_SelfConsistentNoMemory = logLH_EstimateExp3;       
end


end    

function [pCGthetaAll, estimateModel] = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC, expNumber, thetaStim, stdMotor, modelType, lapseRate, includeIncongruentTrials)
if expNumber == 3
    flagDecisionGiven = 1;
else
    flagDecisionGiven = 0;
end
try
    estimateModel.Xval = cell(1, length(stdSensory));
    estimateModel.Yval = cell(1, length(stdSensory)); 
    pCGthetaAll = NaN(length(stdSensory), length(thetaStim));
    
    dstep = 0.1;
    rangeth = [-60 60];
    th = rangeth(1):dstep:rangeth(2);
    nth = length(th);

    pthGC = zeros(2,nth);
    if modelType == 1
        pthGC(1,:) = TukeyWindow([0 priorRange], 0, smoothFactor, th);
        pthGC(2,:) = TukeyWindow([-priorRange 0], 1, smoothFactor, th);
    elseif modelType == 2
        pth = (TukeyWindow([0 priorRange], 0, smoothFactor, th) + TukeyWindow([-priorRange 0], 1, smoothFactor, th))/2;
        pth(th==0) = 0;
        pth(th==0) = max(pth);
        pthGC(1,:) = pth;
        pthGC(2,:) = pth;    
    end
    pthGC_Discrimination(1,:) = TukeyWindow([0 priorRange], 0, smoothFactor, th);
    pthGC_Discrimination(2,:) = TukeyWindow([-priorRange 0], 1, smoothFactor, th);

    if flagDecisionGiven
        pthGC(1,:) = TukeyWindow([0 priorRange], 0, smoothFactor, th);
        pthGC(2,:) = TukeyWindow([-priorRange 0], 1, smoothFactor, th);        
    end
    for kk=1:length(stdSensory)  
        rangeM = [min(thetaStim)-5*stdSensory(kk) max(thetaStim)+5*stdSensory(kk)];
        if rangeM(2) < rangeth(2)
            rangeM = rangeth;
        end
        nm = 1000;
        m = linspace(rangeM(1), rangeM(2), nm);

        nmm = 1200;
        rangeMM = [min(rangeM)-6*stdMemory max(rangeM)+6*stdMemory];
        if rangeMM(2) < rangeth(2)
            rangeMM = rangeth;
        end
        mm = linspace(rangeMM(1), rangeMM(2), nmm);

        M = repmat(m',1,nth);
        MM_m = repmat(mm',1,nm);
        MM_th = repmat(mm',1,nth); 
        THm = repmat(th, nm, 1); 
        THmm = repmat(th, nmm, 1); 
        
        %% Generative (forward)
        % orientation noise
        pmGth = exp(-((M-THm).^2)./(2*stdSensory(kk)^2));
        pmGth = pmGth./(repmat(sum(pmGth,1),nm,1)); 

        %% Inference
        % 1: categorical judgment
        if ~flagDecisionGiven
            PCGm = (pthGC_Discrimination * pmGth') .* repmat(pC,1,nm);
            % fix the issue when sensory noise is too low
            indFirstNonZero = find(PCGm(2,:), 1);
            PCGm(2, 1: indFirstNonZero-1) = PCGm(2, indFirstNonZero);
            indLastNonZero = find(PCGm(1,:), 1, 'last');
            PCGm(1, indLastNonZero+1:end) = PCGm(1, indLastNonZero);
            PCGm = PCGm./(repmat(sum(PCGm,1),2,1));
            % max posterior decision
            PChGm = round(PCGm);
            % marginalization
            PChGtheta = PChGm * pmGth(:, ismember(th, thetaStim));
            PChGtheta_lapse = lapseRate + (1 - 2*lapseRate) * PChGtheta;
            PChGtheta_lapse = PChGtheta_lapse ./ repmat(sum(PChGtheta_lapse, 1), 2, 1);
        else
            PChGtheta_lapse = NaN(2, length(thetaStim));
            PChGtheta_lapse(1, thetaStim > 0) = 1;
            PChGtheta_lapse(1, thetaStim < 0) = 0;
            PChGtheta_lapse(1, thetaStim == 0) = 0.5;
            PChGtheta_lapse(2, :) = 1 - PChGtheta_lapse(1, :);
        end
        
        % 2: estimation
        pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); % p(mm|th) = N(th, sm^2 + smm^2)
        pmmGth = pmmGth./(repmat(sum(pmmGth,1),nmm,1)); 

        pthGmmChcw = (pmmGth.*repmat(pthGC(1,:),nmm,1))';
        pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
        pthGmmChcw(isnan(pthGmmChcw)) = 0;

        pthGmmChccw = (pmmGth.*repmat(pthGC(2,:),nmm,1))';
        pthGmmChccw = pthGmmChccw./repmat(sum(pthGmmChccw,1),nth,1);
        pthGmmChccw(isnan(pthGmmChccw)) = 0;            
        

        EthChcw = th * pthGmmChcw;
        EthChccw = th * pthGmmChccw;
        % discard repeating/decreasing values (required for interpolation) 
        indKeepCw = 1:length(EthChcw);
        while sum(diff(EthChcw)<=0) >0
            indDiscardCw = [false diff(EthChcw)<=0];
            EthChcw(indDiscardCw) = [];
            indKeepCw(indDiscardCw) = [];
        end
        indKeepCcw = 1:length(EthChccw);
        while sum(diff(EthChccw)<=0) >0
            indDiscardCcw = [diff(EthChccw)<=0 false];
            EthChccw(indDiscardCcw) = [];
            indKeepCcw(indDiscardCcw) = [];
        end

        a = 1./gradient(EthChcw,dstep);
        if ~flagDecisionGiven
            % memory noise
            pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemory^2)); 
            pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1));   

            % attention marginalization: compute distribution only over those ms that lead to cw decision!
            pmmGthChcw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim)));
            pmmGthChcw = pmmGthChcw ./ repmat(sum(pmmGthChcw,1),nmm,1);
            b = repmat(a',1,length(thetaStim)) .* pmmGthChcw(indKeepCw, :);        
        else
            b = repmat(a',1,length(thetaStim)) .* pmmGth(indKeepCw, ismember(th, thetaStim));   
        end
        pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChcw(pthhGthChcw < 0) = 0; 

        a = 1./gradient(EthChccw,dstep);
        if ~flagDecisionGiven
            % attention marginalization: compute distribution only over those ms that lead to cw decision!
            pmmGthChccw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim)));  
            b = repmat(a',1,length(thetaStim)) .* pmmGthChccw(indKeepCcw, :);        
        else
            b = repmat(a',1,length(thetaStim)) .* pmmGth(indKeepCcw, ismember(th, thetaStim));
        end  
        pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChccw(pthhGthChccw < 0) = 0; 

        if isempty(includeIncongruentTrials)
            pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
            pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            

            % modify psychometric curve p(Chat|theta, Congruent) ~ p(Congruent| Chat, theta) * p(Chat|Theta)
            pCongruentGcwTh = sum(pthhGthChcw(th' >= 0, :));
            pCongruentGccwTh = sum(pthhGthChccw(th' <= 0, :));
            PChGtheta_lapse = PChGtheta_lapse .* [pCongruentGcwTh; pCongruentGccwTh];
            PChGtheta_lapse = PChGtheta_lapse ./ repmat(sum(PChGtheta_lapse, 1), 2, 1);

            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'>= 0, :) = 0;
            pthhGthChcw(th'< 0, :) = 0;
        end 
        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);
        pthhGth = pthhGthChcw.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw.*repmat(PChGtheta_lapse(2,:),nth,1);

        pCGthetaAll(kk, :) = PChGtheta_lapse(1,:);
        estimateModel.Xval{kk} = th;
        estimateModel.Yval{kk} = pthhGth;
    end       
catch e
    keyboard
    rethrow(e)
end
end