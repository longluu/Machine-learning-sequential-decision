function [fitParameter, negLogLH] = modelFitBayes_Resample(optimizationAlgorithm, SetStartPoint, initialValue, fileID, expNumber, plotFitProgress, modelType, dataName, fixMotorNoise, includeIncongruentTrials, fixLapseRate, fixBoundaryCutoff)
%%%%%%%%%%%% Fitting Bayesian model to data
%%%%%%%%%%%% Condition the PRIOR only
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
minStd = 1;
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
maxLapse = 0.2;

if fixLapseRate
    maxLapse = 0;
end

if fixBoundaryCutoff
    maxBoundaryCutoff = 0;
end

A_Inequality = [1 -1 zeros(1,5) 0;0 1 -1 zeros(1,4) 0];
b_Inequality = [0; 0];
vlb = [minStd minStd  0            10  minStdMemory 1e-2 minStdMotor 0];
vub = [maxStd maxStd  maxLapse     60  maxStdMemory 1   maxStdMotor  maxBoundaryCutoff];
if SetStartPoint
    searchParameter0 = initialValue;
else
    stdSensoryInitial = initialValue(1:2) + 5*(rand(1,2)-0.5);
    lapseRateInitial =  initialValue(3) + 0.01*(rand-0.5);
    priorRangeInitial = initialValue(4) + 8*(rand-0.5);
    stdMemoryInitial = initialValue(5) + 4*(rand-0.5);
    smoothFactorInitial = rand; %(vub(end) -vlb(end)) * rand + vlb(end);
    stdMotorInitial = initialValue(7) + 5*(rand-0.5);
    boundaryCutoffInitial = 2*rand;
    searchParameter0 = [stdSensoryInitial lapseRateInitial priorRangeInitial stdMemoryInitial smoothFactorInitial stdMotorInitial boundaryCutoffInitial];
end

%% Run the optimization
% Define the objective function
optimFun = @(searchParameter, varargin) computeNegLLH(searchParameter(1:2),searchParameter(3),...
     searchParameter(4), searchParameter(5), searchParameter(6), searchParameter(7), searchParameter(8), angleDiff, percentCW, ...
    nTrialsPerCondition, estimateData, fileID, expNumber, plotFitProgress, modelType, includeIncongruentTrials);
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

function errorTotal = computeNegLLH(stdSensory, lapseRate, priorRange, stdMemory, smoothFactor, stdMotor, boundaryCutoff, angleDiff,...
                            percentCW, nTrials, estimateData, fileID, expNumber, plotFitProgress, modelType, includeIncongruentTrials)
global hOptim
global negLLH

try
    %% Compute the full model prediction for both tasks
    pC = [0.5 0.5]';
    [predictedFractionCW, estimateModel] = fullBayesian(stdSensory, stdMemory, priorRange,...
            smoothFactor, boundaryCutoff, pC, expNumber, angleDiff, stdMotor, modelType, lapseRate, includeIncongruentTrials);    
    
    %% Error term for discrimination task
    epsZero = 10^(-10);
    epsDistance = epsZero;
    
    % Compute predicted Percent Correct for each data point
    predictedFractionCW(predictedFractionCW < epsZero) = epsZero;
    predictedFractionCW(predictedFractionCW > 1-epsDistance) = 1-epsDistance;

    % Accumulate log likelihood for the entire data set
    pCW = percentCW/100;
    numberCW = round(nTrials .* pCW);
    numberCCW = nTrials-numberCW;
    logLikelihoodDiscriminate = sum(sum(numberCW .* log(predictedFractionCW) + ...
                        numberCCW .* log(1-predictedFractionCW)));

    %% Error term for estimation task    
    logLikelihoodEstimate = 0;
    tempEstimateModelX = estimateModel.Xval{1};
    for ii = 1 : length(estimateModel.Yval)
        estimateModelY = estimateModel.Yval{ii};
        estimateModelY = estimateModelY(:, angleDiff >=0);
        for jj = 1 : size(estimateModelY,2)            
            tempEstimateModelY = estimateModelY(:, jj); 
            tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
            tempestimateDataX = estimateData{ii,jj};
            tempestimateDataX(tempestimateDataX < min(tempEstimateModelX)) = min(tempEstimateModelX);
            tempestimateDataX(tempestimateDataX > max(tempEstimateModelX)) = max(tempEstimateModelX);               
            tempLikelihood = interp1(tempEstimateModelX, tempEstimateModelY, tempestimateDataX, 'pchip');
            tempLikelihood(tempLikelihood == 0) = epsZero; 
            tempLikelihood(isnan(tempLikelihood)) = [];            
            logLikelihoodEstimate = logLikelihoodEstimate + nansum(log(tempLikelihood)); 
            if logLikelihoodEstimate>0
                keyboard
            end                
        end
    end
    errorTotal = -logLikelihoodDiscriminate - logLikelihoodEstimate;
    disp(['-logLH:' num2str(round(errorTotal)) ' ' num2str(round(-logLikelihoodDiscriminate)) ...
            ' ' num2str(round(-logLikelihoodEstimate)) ', Params: ' num2str([roundn(stdSensory, -2) roundn(lapseRate, -4) roundn([priorRange stdMemory smoothFactor stdMotor boundaryCutoff],-2)])])
    fprintf(fileID,'%2s %9.1f %8.1f %9.1f %9.4f %9.4f  %16.4f  %10.4f %10.4f %8.4f %9.4f %8.4f \r\n', '//', errorTotal, -logLikelihoodDiscriminate,...
        -logLikelihoodEstimate, stdSensory(1), stdSensory(2), lapseRate, priorRange, stdMemory, smoothFactor, stdMotor, boundaryCutoff);
    
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


function [pCGthetaAll, estimateModel] = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, boundaryCutoff, pC, ~, thetaStim, stdMotor, modelType, lapseRate, includeIncongruentTrials)
try
    estimateModel.Xval = cell(1, length(stdSensory));
    estimateModel.Yval = cell(1, length(stdSensory)); 
    pCGthetaAll = NaN(length(stdSensory), length(thetaStim));
    
    pthcw = priorRange;
    pthccw = -priorRange;

    dstep = 0.1;
    rangeth = [-65 65];
    th = rangeth(1):dstep:rangeth(2);
    nth = length(th);

    pthGC = zeros(2,nth);

    if modelType == 1
        pthGC(1,:) = TukeyWindowNew([0 pthcw], smoothFactor, th, boundaryCutoff);
        pthGC(2,:) = TukeyWindowNew([pthccw 0], smoothFactor, th, boundaryCutoff);
    elseif modelType == 2
        pth = (TukeyWindow([0 pthcw], 0, smoothFactor, th) + TukeyWindow([pthccw 0], 1, smoothFactor, th))/2;
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
        
        nmr = 1200;
        rangeMR = [min(rangeM)-6*stdMemory max(rangeM)+6*stdMemory];
        if rangeMR(2) < rangeth(2)
            rangeMR = rangeth;
        end        
        mr = linspace(rangeMR(1), rangeMR(2), nmr);

        M = repmat(m',1,nth);
        MR_th = repmat(mr', 1, nth);
        THm = repmat(th, nm, 1); 
        THmr = repmat(th, nmr, 1);

        %% Correct trials
        % Generative (forward)
        % orientation noise
        pmGth = exp(-((M-THm).^2)./(2*stdSensory(kk)^2));
        pmGth = pmGth./(repmat(sum(pmGth,1),nm,1)); 

        % Inference
        % 1: categorical judgment
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

        % 2: estimation
        % Likelihood function the same as correct decision p(mm|th) = N(th, sm^2 + smm^2)           
        pmrGth = exp(-((MR_th-THmr).^2)./(2*(stdSensory(kk)^2 + stdMemory^2)));
        pmrGth = pmrGth./(repmat(sum(pmrGth,1),nmr,1)); 
        pthGmrChcw = (pmrGth.*repmat(pthGC(1,:),nmr,1))';
        pthGmrChcw = pthGmrChcw./repmat(sum(pthGmrChcw,1),nth,1);
        pthGmrChcw(isnan(pthGmrChcw)) = 0;

        pthGmrChccw = (pmrGth.*repmat(pthGC(2,:),nmr,1))';
        pthGmrChccw = pthGmrChccw./repmat(sum(pthGmrChccw,1),nth,1);
        pthGmrChccw(isnan(pthGmrChccw)) = 0;

        EthChcw = th * pthGmrChcw;
        EthChccw = th * pthGmrChccw;
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

        % Resample m until we have a sample that is consistent with feedback
        % p(mr|m, theta, Chat)
        MR_m = repmat(mr', 1, nm);
        pmrGmth = exp(-((MR_m-repmat(m, nmr, 1)).^2)./(2*(stdMemory^2))); 

        pmrGmthChcw = pmrGmth;
        pmrGmthChcw(mr < 0, :) = 0;
        % put the tail with all 0 to 1 (deal with small memory noise)
        indZero = sum(pmrGmthChcw, 1) == 0;
        pmrGmthChcw(mr > 0, indZero) = 1;
        pmrGmthChcw = pmrGmthChcw./(repmat(sum(pmrGmthChcw,1),nmr,1));
        pmrGmthChcw(mr > 0, indZero) = 1e-50;

        pmrGmthChccw = pmrGmth;
        pmrGmthChccw(mr > 0, :) = 0;
        % put the tail with all 0 to 1 (deal with small memory noise)
        indZero = sum(pmrGmthChccw, 1) == 0;
        pmrGmthChccw(mr < 0, indZero) = 1;        
        pmrGmthChccw = pmrGmthChccw./(repmat(sum(pmrGmthChccw,1),nmr,1));
        pmrGmthChccw(mr < 0, indZero) = 1e-50;

        % Marginalize over m that lead to cw/ccw decision to compute likelihood p(mr|theta, Chat)
        pmGthChcw = pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim));
        pmrGthChcw = pmrGmthChcw * pmGthChcw;   
        pmrGthChcw = pmrGthChcw ./ (repmat(sum(pmrGthChcw,1),nmr,1)); 
        pmrGthChcw(isnan(pmrGthChcw)) = 0;

        pmGthChccw = pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim));
        pmrGthChccw = pmrGmthChccw * pmGthChccw;
        pmrGthChccw = pmrGthChccw ./ (repmat(sum(pmrGthChccw,1),nmr,1)); 
        pmrGthChccw(isnan(pmrGthChccw)) = 0;

        a = 1./gradient(EthChcw,dstep);
        b = repmat(a',1,length(thetaStim)) .* pmrGthChcw(indKeepCw, :);        
        pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChcw(pthhGthChcw < 0) = 0; 

        a = 1./gradient(EthChccw,dstep);
        b = repmat(a',1,length(thetaStim)) .* pmrGthChccw(indKeepCcw, :);        
        pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChccw(pthhGthChccw < 0) = 0; 

        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);             

        if includeIncongruentTrials == 0
            % modify psychometric curve p(Chat|theta, Congruent) ~ p(Congruent| Chat, theta) * p(Chat|Theta)
            pCongruentGcwTh = sum(pthhGthChcw(th' >= 0, :));
            pCongruentGccwTh = sum(pthhGthChccw(th' <= 0, :));
            PChGtheta_lapse_new = PChGtheta_lapse .* [pCongruentGcwTh; pCongruentGccwTh];
            PChGtheta_lapse_new = PChGtheta_lapse_new ./ repmat(sum(PChGtheta_lapse_new, 1), 2, 1);

            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'>= 0, :) = 0;
            pthhGthChcw(th'< 0, :) = 0;
        else
            PChGtheta_lapse_new = PChGtheta_lapse;
        end

        % remove incorrect trials
        pthhGthChcw(:, thetaStim < 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          hcw(:, thetaStim < 0) = 0;
        pthhGthChccw(:, thetaStim > 0) = 0;

        pthhGth_correct = pthhGthChcw.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw.*repmat(PChGtheta_lapse(2,:),nth,1);
    
        pCGthetaAll(kk, :) = PChGtheta_lapse_new(1,:);
        estimateModel.Xval{kk} = th;
        estimateModel.Yval{kk} = pthhGth_correct;    
    end
    
catch e
    keyboard
    rethrow(e)
end
end

function change_current_figure(h)
set(0,'CurrentFigure',h)
end
