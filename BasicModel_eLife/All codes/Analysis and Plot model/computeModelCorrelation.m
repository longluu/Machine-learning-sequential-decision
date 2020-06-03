%%%%%%%%%%%%% Compute the log likelihood of self-consistent,  %%%%%%%%%%%%%
%%%%%%%%%%%%% standard Bayes and oracle model %%%%%%%%%%%%%
function computeModelCorrelation
subjectID = {'ll', 'll', 'sy', 'cz', 'vs', 'as', 'xfl', 'aj', 'zw', 'skb'};
expNum = [1 2 1 1 1 1 2 2 2 2];        
paramsSC1 = [2.5611    4.8570    7.1153      0.0000    28.6016      0.8805    0.8461    2.6327;
             2.5611    4.8570    7.1153      0.0000    23.1044    0.8805    0.8461    2.6327;
             4.0196    5.3326    8.1636           0.0000     41.9756     2.2070   0.0100    1.8431;
             3.8789    6.4018   12.3890           0.0000     40.8589    14.1001   0.3315    2.9146;
             3.6687    4.6775   10.9523           0.0000     41.4502     0.1422   0.7650    5.4045;
             3.9257    5.7877   10.9306           0.0000     38.8974    14.2712   0.2928    5.8917;
             6.3340    8.4236   14.6740           0.0000     22.2427     6.2871   1.0000    2.9350;
             4.5991    6.0397    7.8831           0.0000     15.9917     1.4014   0.0120    2.9916;
             4.5848    5.4837    6.7755           0.0000     17.5182     4.3647   0.0119    2.5348;
             4.3183    4.5437    8.4461           0.0000     33.0670    14.7435   0.7733    2.0456];
         
biasDataExp1 = [];
biasModelExp1 = [];
biasDataExp2 = [];
biasModelExp2 = [];
biasDataCollapseExp2All = [];
biasDataCollapseExp3All = [];
biasModelCollapseExp2All = [];
biasModelCollapseExp3All = [];

for ii = 1 : size(paramsSC1, 1)
    [biasData, biasModel, biasDataCollapseExp2, biasDataCollapseExp3, biasModelCollapseExp2, biasModelCollapseExp3] = ...
                    computeBias(subjectID(ii), expNum(ii), paramsSC1(ii, :));
    if expNum(ii) == 1
        biasDataExp1 = [biasDataExp1; biasData];
        biasModelExp1 = [biasModelExp1; biasModel];
    else
        biasDataExp2 = [biasDataExp2; biasData];
        biasModelExp2 = [biasModelExp2; biasModel];
        biasDataCollapseExp2All = [biasDataCollapseExp2All; biasDataCollapseExp2];
        biasDataCollapseExp3All = [biasDataCollapseExp3All; biasDataCollapseExp3];
        biasModelCollapseExp2All = [biasModelCollapseExp2All; biasModelCollapseExp2];
        biasModelCollapseExp3All = [biasModelCollapseExp3All; biasModelCollapseExp3];                
    end        
end


%% Compare the bias data vs. model
figure
subplot(1, 2, 1)
hold on
mseExp1 = sqrt(mean((biasModelExp1 - biasDataExp1).^2));
minAxis = min([biasModelExp1; biasDataExp1]) - 1;
maxAxis = max([biasModelExp1; biasDataExp1]) + 1;
plot(biasModelExp1, biasDataExp1, 'o')
plot([minAxis maxAxis], [minAxis maxAxis])
ylim([minAxis maxAxis])
xlim([minAxis maxAxis])
ylabel('Mean bias (data)')
xlabel('Mean bias (self-consistent model)')
title(['Experiment 1, RMSE: ' num2str(mseExp1)])

subplot(1, 2, 2)
hold on
mseExp2 = sqrt(mean((biasModelExp2 - biasDataExp2).^2));
minAxis = min([biasModelExp2; biasDataExp2]) - 1;
maxAxis = max([biasModelExp2; biasDataExp2]) + 1;
plot(biasModelExp2, biasDataExp2, 'o')
plot([minAxis maxAxis], [minAxis maxAxis])
ylim([minAxis maxAxis])
xlim([minAxis maxAxis])
ylabel('Mean bias (data)')
xlabel('Mean bias (self-consistent model)')
title(['Experiment 2 and 3, RMSE: ' num2str(mseExp2)])

figure
hold on
diffDataModelExp1 = biasDataExp1 - biasModelExp1;
diffDataModelExp2 = biasDataExp2 - biasModelExp2;
ciExp1 = abs(bootci(2000, @mean, diffDataModelExp1) - mean(diffDataModelExp1));
ciExp2 = abs(bootci(2000, @mean, diffDataModelExp2) - mean(diffDataModelExp2));

plot(ones(1, length(diffDataModelExp1)), diffDataModelExp1, 'o')
plot(2*ones(1, length(diffDataModelExp2)), diffDataModelExp2, 'o')
errorbar([1 2], [mean(diffDataModelExp1) mean(diffDataModelExp2)], [ciExp1(1) ciExp2(1)], [ciExp1(2) ciExp2(2)], 'o')
plot([0 3], [0 0], '--')
xlim([0 3])
xlabel('Experiment')
ylabel('Data bias - Model bias (deg)')

%% Compare the bias between Exp 2 and Exp 3
% figure
% subplot(1, 2, 1)
% hold on
% minAxis = min([biasDataCollapseExp2All; biasDataCollapseExp3All]) - 1;
% maxAxis = max([biasDataCollapseExp2All; biasDataCollapseExp3All]) + 1;
% plot(biasDataCollapseExp2All, biasDataCollapseExp3All, 'o')
% plot([minAxis maxAxis], [minAxis maxAxis])
% ylim([minAxis maxAxis])
% xlim([minAxis maxAxis])
% ylabel('Mean bias (Experiment 3)')
% xlabel('Mean bias (Experiment 2)')
% 
% subplot(1, 2, 2)
% hold on
% minAxis = min([biasModelCollapseExp2All; biasModelCollapseExp3All]) - 1;
% maxAxis = max([biasModelCollapseExp2All; biasModelCollapseExp3All]) + 1;
% plot(biasModelCollapseExp2All, biasModelCollapseExp3All, 'o')
% plot([minAxis maxAxis], [minAxis maxAxis])
% ylim([minAxis maxAxis])
% xlim([minAxis maxAxis])
% ylabel('Mean bias (Experiment 3)')
% xlabel('Mean bias (Experiment 2)')


end

function [biasData, biasModel, biasDataExp2, biasDataExp3, biasModelExp2, biasModelExp3] = computeBias(subjectID, expNumber, paramsSelfConsistent)
includeIncongruentTrials = ''; % empty if not include incongruent trials Incongruent
biasDataExp2 = NaN(3, 1);
biasDataExp3 = NaN(3, 1);
biasModelExp2 = NaN(3, 1);
biasModelExp3 = NaN(3, 1);
nPointCollapse = 3;
if expNumber == 1
    %##################################### Experiment 1 #####################################    
    %% Bias of the data
    [~, ~, ~, estimateData, angleDiff, ~] = dataForFitting(subjectID, expNumber, includeIncongruentTrials);
    biasData = NaN(numel(estimateData), 1);
    counter = 1;
    angleDiffCollapse = angleDiff(angleDiff>=0);
    for ii = 1 : size(estimateData, 1)
        for jj = 1 : size(estimateData, 2)
            tempEst = estimateData{ii, jj};
            biasData(counter) = nanmean(tempEst(tempEst>0)) - angleDiffCollapse(jj);
            counter = counter + 1;
        end
    end
    
    %% Bias of Self-consistent Bayes
    % Fit paramters
    paramsAll = paramsSelfConsistent;
    lapseRate = paramsAll(4);
    stdSensory = paramsAll(1:3);
    stdMemory = paramsAll(6);
    stdMotor = paramsAll(8);
    pC = [0.5, 0.5]'; % [cw ccw]
    priorRange = paramsAll(5);
    smoothFactor = paramsAll(7);

    % Run the model
    meanEstimateModel = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC,...
                                expNumber, angleDiff, stdMotor, 1, lapseRate, includeIncongruentTrials);

    meanBias = meanEstimateModel - repmat(angleDiff, length(stdSensory), 1);
    biasModel = (meanBias(:, angleDiff >= 0))';
    biasModel = biasModel(:);
else
    %##################################### Experiment 2 and 3 #####################################

    %% Extract the data
    [~, ~, ~, estimateDataExp2, angleDiff, ~] = dataForFitting(subjectID, 2, includeIncongruentTrials);
    [~, ~, ~, estimateDataExp3, ~] = dataForFitting(subjectID, 3, includeIncongruentTrials);
    biasData = NaN(2*numel(estimateDataExp2), 1);
    counter = 1;
    angleDiffCollapse = angleDiff(angleDiff>=0);
    tempBias = [];
    for ii = 1 : size(estimateDataExp2, 1)
        for jj = 1 : size(estimateDataExp2, 2)
            tempEst = estimateDataExp2{ii, jj};
            biasData(counter) = nanmean(tempEst(tempEst>0)) - angleDiffCollapse(jj);
            if jj <= nPointCollapse
                tempBias = [tempBias; biasData(counter)];
            end
            counter = counter + 1;            
        end
    end
    tempBias = reshape(tempBias, [], size(estimateDataExp2, 1));
    biasDataExp2 = mean(tempBias', 2);
    
    tempBias = [];
    for ii = 1 : size(estimateDataExp3, 1)
        for jj = 1 : size(estimateDataExp3, 2)
            tempEst = estimateDataExp3{ii, jj};
            biasData(counter) = nanmean(tempEst(tempEst>0)) - angleDiffCollapse(jj);
            if jj <= nPointCollapse
                tempBias = [tempBias; biasData(counter)];
            end
            counter = counter + 1;            
        end
    end
    tempBias = reshape(tempBias, [], size(estimateDataExp2, 1));
    biasDataExp3 = mean(tempBias', 2);
    

    %% Self-consistent Bayes
    % Fit paramters
    paramsAll = paramsSelfConsistent;
    lapseRate = paramsAll(4);
    stdSensory = paramsAll(1:3);
    stdMemory = paramsAll(6);
    stdMotor = paramsAll(8);
    pC = [0.5, 0.5]'; % [cw ccw]
    priorRange = paramsAll(5);
    smoothFactor = paramsAll(7);

    % Run the model
    meanEstimateModelExp2 = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC,...
                                2, angleDiff, stdMotor, 1, lapseRate, includeIncongruentTrials);
    meanEstimateModelExp3 = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC,...
                                3, angleDiff, stdMotor, 1, lapseRate, includeIncongruentTrials);

    meanBias2 = meanEstimateModelExp2 - repmat(angleDiff, length(stdSensory), 1);
    biasModel2 = (meanBias2(:, angleDiff >= 0))';
    biasModelExp2 = (mean(biasModel2(1:nPointCollapse, :), 1))';
    biasModel2 = biasModel2(:);
    meanBias3 = meanEstimateModelExp3 - repmat(angleDiff, length(stdSensory), 1);
    biasModel3 = (meanBias3(:, angleDiff >= 0))';
    biasModelExp3 = (mean(biasModel3(1:nPointCollapse, :), 1))';    
    biasModel3 = biasModel3(:);
    biasModel = [biasModel2; biasModel3];
end
end    

function [meanEstimate] = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC, expNumber, thetaStim, stdMotor, modelType, lapseRate, includeIncongruentTrials)
if expNumber == 3
    flagDecisionGiven = 1;
else
    flagDecisionGiven = 0;
end
try
    meanEstimate = NaN(length(stdSensory), length(thetaStim));
    
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
            PCGm = (pthGC * pmGth') .* repmat(pC,1,nm);
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

        if isempty(includeIncongruentTrials)
            pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    

            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChcw(th'< 0, :) = 0;
        end 
        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not
        meanEstimate(kk, :) = th * pthhGthChcw;
    end       
catch e
    keyboard
    rethrow(e)
end
end