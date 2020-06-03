function [fitParameter, negLogLH] = modelFitBayesJoint12(optimizationAlgorithm, SetStartPoint, initialValue, fileID, expNumber, plotFitProgress, modelType, dataName, fixMotorNoise, includeIncongruentTrials, fixLapseRate)
%%%%%%%%%%%% Fitting Bayesian model to data (including the incongruent trials %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%
global hOptim
global negLLH

if plotFitProgress
    hOptim = figure;
end
negLLH = [];


%% Initialize some variable
if ~isempty(dataName)
    load(dataName)
else
    load dataAll
end
maxFunEval = 350;
tolFun = 1e-1;
tolX = 1e-2;
minStd = 1.2;
maxStd = 30;
if fixMotorNoise
    minStdMotor = stdMotor;
    maxStdMotor = stdMotor;
else
    minStdMotor = 0.01;
    maxStdMotor = 20;
end
minStdMemory = 0.005;
maxStdMemory = 20;
maxLapse = 0.15;
if fixLapseRate
    maxLapse = 0;
end
A_Inequality = [1 -1 zeros(1,5) 0;0 1 -1 zeros(1,4) 0];
b_Inequality = [0; 0];
vlb = [minStd minStd minStd  0            10  10  minStdMemory 1e-2 minStdMotor];
vub = [maxStd maxStd maxStd  maxLapse     42  42  maxStdMemory 1   maxStdMotor];
if SetStartPoint
    searchParameter0 = initialValue;
else
    stdSensoryInitial = initialValue(1:3) + 5*(rand(1,3)-0.5);
    lapseRateInitial =  initialValue(4) + 0.01*(rand-0.5);
    priorRangeInitial1 = initialValue(5) + 6*(rand-0.5);
    priorRangeInitial2 = initialValue(6) + 6*(rand-0.5);
    stdMemoryInitial = initialValue(7) + 4*(rand-0.5);
    smoothFactorInitial = rand; %(vub(end) -vlb(end)) * rand + vlb(end);
    stdMotorInitial = initialValue(9) + 5*(rand-0.5);
    searchParameter0 = [stdSensoryInitial lapseRateInitial priorRangeInitial1 priorRangeInitial2 stdMemoryInitial smoothFactorInitial stdMotorInitial];
end

%% Run the optimization
% Define the objective function
optimFun = @(searchParameter, varargin) computeNegLLH(searchParameter(1:3),searchParameter(4),...
     searchParameter(5), searchParameter(6), searchParameter(7), searchParameter(8), searchParameter(9), angleDiff, percentCW1, ...
    nTrialsPerCondition1, percentCW2, nTrialsPerCondition2, estimateDataExp1, estimateDataExp2, fileID, expNumber, plotFitProgress, modelType, includeIncongruentTrials);
switch optimizationAlgorithm
    case 1
        % using Nealder-Mead simplex
        searchOptions = struct(...
            'Display','Iter',...
            'TolFun',tolFun,...
            'TolX',tolX,...
            'MaxFunEvals',maxFunEval,...
            'FunValCheck','off');
        [fitParameter, negLogLH] = fminsearchbnd(optimFun, searchParameter0, vlb, vub, searchOptions,[]);
    case 2        
        % using Sequential quadratic programming 
        options = optimset('fmincon');
        options = optimset(options,'Diagnostics','off','Display','Iter','LargeScale','off',...
                    'Algorithm','sqp', 'MaxFunEvals', maxFunEval, 'TolFun',tolFun, 'TolX',tolX);
        [fitParameter, negLogLH] = fmincon(optimFun, searchParameter0,A_Inequality,b_Inequality,[],[],vlb,vub,[],options);
    case 3
        % using simulated annealing
        options = saoptimset('simulannealbnd');
        options = saoptimset(options, 'Display','Iter', 'MaxFunEvals', maxFunEval, 'TolFun',tolFun);
        [fitParameter, negLogLH] = simulannealbnd(optimFun,searchParameter0,vlb,vub,options);
    case 4
        % using genetic algorithm
        options = gaoptimset(@ga);
        options = gaoptimset(options, 'Display','iter', 'PopulationSize', 30, 'TolFun',tolFun);
        [fitParameter, negLogLH] = ga(optimFun,8,A_Inequality,b_Inequality,[],[],vlb,vub,[],options);
end

% %% Plot the fit curve and data
% if expNumber ~= 3
%     angleDiffResampled = linspace(angleDiff(1),angleDiff(end),50);
%     fitpercentCW = 100*BayesFractionCW(fitParameter(6), fitParameter(8), fitParameter(1:3), fitParameter(4), fitParameter(5), angleDiffResampled, windowPrior);
%     figure
%     legendName = cell(1,3);
%     colorName = {'g','r', 'b', 'cyan', 'magenta', 'y'};
%     lineWidth = 3;
%     fontSize = 20;
%     set(gca,'FontSize',15)
%     hold on
%     for ii = 1 : 3
%         hLegend(ii) = plot(angleDiffResampled,fitpercentCW(ii,:),'Color',...
%             colorName{ii},'LineWidth',3);
%         plot(angleDiff, percentCW(ii,:), [colorName{ii} 'o'], 'MarkerSize',13,...
%             'MarkerFaceColor', colorName{ii}, 'LineWidth', lineWidth);
%         legendName{ii} = ['noise = ' num2str(roundn(fitParameter(ii),-1))];
%     end
%     xlim([min(angleDiff) max(angleDiff)])
%     xlabel('True angle (degree)','FontSize',fontSize)
%     ylabel('Percent CW','FontSize',fontSize)
%     legend(hLegend, legendName, 'Location', 'SouthEast')
% end
end

function errorTotal = computeNegLLH(stdSensory, lapseRate, priorRange1, priorRange2, stdMemory, smoothFactor, stdMotor, angleDiff,...
                            percentCW1, nTrials1, percentCW2, nTrials2, estimateDataExp1, estimateDataExp2, fileID, expNumber, plotFitProgress, modelType, includeIncongruentTrials)
global hOptim
global negLLH

try
    %% Compute the full model prediction for Exp 1
    pC = [0.5 0.5]';
    [predictedFractionCW1, estimateModelExp1] = fullBayesian(stdSensory, stdMemory, priorRange1,...
            smoothFactor, pC, expNumber(1), angleDiff, stdMotor, modelType, lapseRate, includeIncongruentTrials);    
    
    % Error term for discrimination task
    epsZero = 10^(-10);
    epsDistance = epsZero;
    predictedFractionCW1(predictedFractionCW1 < epsZero) = epsZero;
    predictedFractionCW1(predictedFractionCW1 > 1-epsDistance) = 1-epsDistance;

    % Accumulate log likelihood for the entire data set
    pCW = percentCW1/100;
    numberCW = round(nTrials1 .* pCW);
    numberCCW = nTrials1-numberCW;
    logLikelihoodDiscriminate1 = sum(sum(numberCW .* log(predictedFractionCW1) + ...
                        numberCCW .* log(1-predictedFractionCW1)));

    % Error term for estimation task    
    logLikelihoodEstimate1 = 0;
    tempEstimateModelX = estimateModelExp1.Xval{1};
    for ii = 1 : length(estimateModelExp1.Yval)
        estimateModelY = estimateModelExp1.Yval{ii};
        estimateModelY = estimateModelY(:, angleDiff >=0);
        for jj = 1 : size(estimateModelY,2)            
            tempEstimateModelY = estimateModelY(:, jj); 
            tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
            tempestimateDataX = estimateDataExp1{ii,jj};
            tempestimateDataX(tempestimateDataX<min(tempEstimateModelX) | tempestimateDataX>max(tempEstimateModelX)) = [];                
            tempLikelihood = interp1(tempEstimateModelX, tempEstimateModelY, tempestimateDataX, 'pchip');
            tempLikelihood(tempLikelihood == 0) = epsZero; 
            tempLikelihood(isnan(tempLikelihood)) = [];            
            logLikelihoodEstimate1 = logLikelihoodEstimate1 + nansum(log(tempLikelihood)); 
            if logLikelihoodEstimate1>0
                keyboard
            end                
        end
    end
    errorTotal1 = -logLikelihoodDiscriminate1 - logLikelihoodEstimate1;
    
    %% Compute the full model prediction for both Exp 2 and Exp 3
    pC = [0.5 0.5]';
    [predictedFractionCW2, estimateModelExp2] = fullBayesian(stdSensory, stdMemory, priorRange2,...
            smoothFactor, pC, expNumber(2), angleDiff, stdMotor, modelType, lapseRate, includeIncongruentTrials);    

    % Error term for discrimination task Exp 2
    epsZero = 10^(-10);
    epsDistance = epsZero;
    predictedFractionCW2(predictedFractionCW2 < epsZero) = epsZero;
    predictedFractionCW2(predictedFractionCW2 > 1-epsDistance) = 1-epsDistance;
    pCW = percentCW2/100;
    numberCW = round(nTrials2 .* pCW);
    numberCCW = nTrials2-numberCW;
    logLikelihoodDiscriminate2 = sum(sum(numberCW .* log(predictedFractionCW2) + ...
                        numberCCW .* log(1-predictedFractionCW2)));

    % Error term for estimation task Exp 2   
    logLikelihoodEstimate2 = 0;
    tempEstimateModelX = estimateModelExp2.Xval{1};
    for ii = 1 : length(estimateModelExp2.Yval)
        estimateModelY = estimateModelExp2.Yval{ii};
        estimateModelY = estimateModelY(:, angleDiff >=0);
        for jj = 1 : size(estimateModelY,2)            
            tempEstimateModelY = estimateModelY(:, jj); 
            tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
            tempestimateDataX = estimateDataExp2{ii,jj};
            tempestimateDataX(tempestimateDataX<min(tempEstimateModelX) | tempestimateDataX>max(tempEstimateModelX)) = [];                
            tempLikelihood = interp1(tempEstimateModelX, tempEstimateModelY, tempestimateDataX, 'pchip');
            tempLikelihood(isnan(tempLikelihood)) = [];
            tempLikelihood(tempLikelihood == 0) = epsZero;            
            logLikelihoodEstimate2 = logLikelihoodEstimate2 + nansum(log(tempLikelihood)); 
            if logLikelihoodEstimate2>0
                keyboard
            end                
        end
    end
    errorTotal2 = -logLikelihoodDiscriminate2 - logLikelihoodEstimate2;
    

    errorTotal = errorTotal1 + errorTotal2;
    
    
    disp(['-logLH:' num2str(round(errorTotal)) ' ' num2str(round(-logLikelihoodDiscriminate1-logLikelihoodDiscriminate2)) ...
            ' ' num2str(round(-logLikelihoodEstimate1-logLikelihoodEstimate2)) ', Params: ' num2str([roundn(stdSensory, -2) roundn(lapseRate, -4) roundn([priorRange1 priorRange2 stdMemory smoothFactor stdMotor],-2)])])
    fprintf(fileID,'%2s %9.1f %8.1f %9.1f %9.4f %9.4f %9.4f %16.4f  %10.4f %10.4f %10.4f %8.4f %9.4f \r\n', '//', errorTotal, -logLikelihoodDiscriminate1-logLikelihoodDiscriminate2,...
        -logLikelihoodEstimate1-logLikelihoodEstimate2,...
        stdSensory(1), stdSensory(2), stdSensory(3), lapseRate, priorRange1, priorRange2, stdMemory, smoothFactor, stdMotor);
    
    %% Plot the optimization progress
    if plotFitProgress
        negLLH = [negLLH errorTotal];
        change_current_figure(hOptim)
        set(hOptim, 'Units', 'normalized', 'Position', [0.02 0.6 0.96 0.38])
        hold on
        maxY = max(negLLH);
        minY = min(negLLH);
        if minY == maxY
            maxY = minY + 1;
        end
        plot(negLLH,'-o', 'MarkerFaceColor','b')
        ylim([minY maxY])
        xlabel('Search point', 'FontSize', 25)
        ylabel('Negative Log LH', 'FontSize', 25)
        drawnow
    end
catch e
    keyboard
    rethrow(e)
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
    rangeth = [-42 42];
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
        pthhGthChcw(pthhGthChcw < 0) = 0; 
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
        pthhGthChccw(pthhGthChccw < 0) = 0; 
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

function change_current_figure(h)
set(0,'CurrentFigure',h)
end
