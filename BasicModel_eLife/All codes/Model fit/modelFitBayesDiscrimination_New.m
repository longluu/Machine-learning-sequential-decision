function [fitParameter, negLogLH] = modelFitBayesDiscrimination_New(optimizationAlgorithm, SetStartPoint, initialValue, fileID, expNumber, plotFitProgress, modelType, dataName, fixMotorNoise, includeIncongruentTrials, fixLapseRate)
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
minStdMotor = 0.001;
maxStdMotor = 0.001;
minStdMemory = 0.001;
maxStdMemory = 0.001;
maxLapse = 0.2;
if fixLapseRate
    maxLapse = 0;
end
A_Inequality = [1 -1 zeros(1,5) 0;0 1 -1 zeros(1,4) 0];
b_Inequality = [0; 0];
vlb = [initialValue(1) initialValue(2) initialValue(3)  0            initialValue(5)  minStdMemory initialValue(7) minStdMotor];
vub = [initialValue(1) initialValue(2) initialValue(3)  maxLapse     initialValue(5)  maxStdMemory initialValue(7) maxStdMotor];

if SetStartPoint
    searchParameter0 = initialValue;
else
    stdSensoryInitial = initialValue(1:3) + 5*(rand(1,3)-0.5);
    lapseRateInitial =  initialValue(4) + 0.01*(rand-0.5);
    priorRangeInitial = initialValue(5) + 6*(rand-0.5);
    stdMemoryInitial = initialValue(6) + 4*(rand-0.5);
    smoothFactorInitial = rand; %(vub(end) -vlb(end)) * rand + vlb(end);
    stdMotorInitial = initialValue(8) + 5*(rand-0.5);
    searchParameter0 = [stdSensoryInitial lapseRateInitial priorRangeInitial stdMemoryInitial smoothFactorInitial stdMotorInitial];
end

%% Run the optimization
% Define the objective function
optimFun = @(searchParameter, varargin) computeNegLLH(searchParameter(1:3),searchParameter(4),...
     searchParameter(5), searchParameter(6), searchParameter(7), searchParameter(8), angleDiff, percentCW, ...
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

function errorTotal = computeNegLLH(stdSensory, lapseRate, priorRange, stdMemory, smoothFactor, stdMotor, angleDiff,...
                            percentCW, nTrials, estimateData, fileID, expNumber, plotFitProgress, modelType, includeIncongruentTrials)
global hOptim
global negLLH

try
    %% Compute the full model prediction for both tasks
    pC = [0.5 0.5]';
    [predictedFractionCW] = fullBayesian(stdSensory, stdMemory, priorRange,...
            smoothFactor, pC, expNumber, angleDiff, stdMotor, modelType, lapseRate, includeIncongruentTrials);    
    
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
    
    errorTotal = -logLikelihoodDiscriminate - logLikelihoodEstimate;
    disp(['-logLH:' num2str(round(errorTotal)) ' ' num2str(round(-logLikelihoodDiscriminate)) ...
            ' ' num2str(round(-logLikelihoodEstimate)) ', Params: ' num2str([roundn(stdSensory, -2) roundn(lapseRate, -4) roundn([priorRange stdMemory smoothFactor stdMotor],-2)])])
    fprintf(fileID,'%9.1f %8.1f %9.1f %9.2f %9.2f %9.2f %16.4f  %10.2f %10.2f %8.2f %9.2f \r\n', errorTotal, -logLikelihoodDiscriminate,...
        -logLikelihoodEstimate, stdSensory(1), stdSensory(2), stdSensory(3), lapseRate, priorRange, stdMemory, smoothFactor, stdMotor);
    
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


function [pCGthetaAll] = fullBayesian(stdSensory, stdMemory, priorRange, smoothFactor, pC, expNumber, thetaStim, stdMotor, modelType, lapseRate, includeIncongruentTrials)
if expNumber == 3
    flagDecisionGiven = 1;
else
    flagDecisionGiven = 0;
end
try
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

    for kk=1:length(stdSensory)  
        rangeM = [min(thetaStim)-5*stdSensory(kk) max(thetaStim)+5*stdSensory(kk)];
        if rangeM(2) < rangeth(2)
            rangeM = rangeth;
        end
        nm = 1000;
        m = linspace(rangeM(1), rangeM(2), nm);

        M = repmat(m',1,nth);
        THm = repmat(th, nm, 1); 
        
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
        

        pCGthetaAll(kk, :) = PChGtheta_lapse(1,:);
    end       
catch e
    keyboard
    rethrow(e)
end
end

function change_current_figure(h)
set(0,'CurrentFigure',h)
end
