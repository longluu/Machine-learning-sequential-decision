function [fitParameter, negLogLH, hLegend] = FullBayesFit(percentCW,nTrialsPerCondition,angleDiff, optimizationAlgorithm, windowPrior, SetStartPoint, initialValue)
global hOptim
global negLLH
hOptim = figure;
negLLH = [];

%% Initialize some variable
maxFunEval = 500;
tolFun = 1e-2;
tolX = 1e-3;
minStd = 1;
maxStd = 20;
maxLapse = 0.02;
A_Inequality = [1 -1 zeros(1,5);0 1 -1 zeros(1,4)];
b_Inequality = [0; 0];
if windowPrior == 1
    vlb = [minStd minStd minStd  0     0     1  1e-2];
    vub = [maxStd maxStd maxStd  maxLapse   maxLapse   42 1];
else
    vlb = [minStd minStd minStd  0     0     1  1e-2];
    vub = [maxStd maxStd maxStd  maxLapse   maxLapse   42 20];    
end
if SetStartPoint
    searchParameter0 = initialValue;
else
    stdNoiseInitial = [3.6 5.5 9.3]+2*(rand(1,3)-0.5);
    lapseRateInitial = maxLapse * rand;
    guessRateInitial = maxLapse * rand;
    priorRangeInitial = 22+6*(rand-0.5);
    smoothFactorInitial = vub(end) * rand;
    searchParameter0 = [stdNoiseInitial lapseRateInitial guessRateInitial priorRangeInitial smoothFactorInitial];
end

%% Run the optimization
% Define the objective function
optimFun = @(searchParameter, varargin) computeNegLLH(searchParameter(1:3),searchParameter(4),...
    searchParameter(5), searchParameter(6), searchParameter(7), angleDiff, percentCW, nTrialsPerCondition, windowPrior);
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

%% Plot the fit curve and data
angleDiffResampled = linspace(angleDiff(1),angleDiff(end),50);
fitpercentCW = 100*BayesFractionCW(fitParameter(6), fitParameter(7), fitParameter(1:3), fitParameter(4), fitParameter(5), angleDiffResampled, windowPrior);
figure
legendName = cell(1,3);
colorName = {'g','r', 'b', 'cyan', 'magenta', 'y'};
lineWidth = 3;
fontSize = 20;
set(gca,'FontSize',15)
hold on
for ii = 1 : 3
    hLegend(ii) = plot(angleDiffResampled,fitpercentCW(ii,:),'Color',...
        colorName{ii},'LineWidth',3);
    plot(angleDiff, percentCW(ii,:), [colorName{ii} 'o'], 'MarkerSize',13,...
        'MarkerFaceColor', colorName{ii}, 'LineWidth', lineWidth);
    legendName{ii} = ['noise = ' num2str(roundn(fitParameter(ii),-1))];
end
xlim([min(angleDiff) max(angleDiff)])
xlabel('True angle (degree)','FontSize',fontSize)
ylabel('Percent CW','FontSize',fontSize)
legend(hLegend, legendName, 'Location', 'SouthEast')

end

function negLogLH = computeNegLLH(stdNoise, lapseRate, guessRate, priorRange, smoothFactor, angleDiff, percentCW, nTrials, windowPrior)
global hOptim
global negLLH

% Compute predicted Percent Correct for each data point
predictedFractionCW = BayesFractionCW(priorRange, smoothFactor, stdNoise, lapseRate, guessRate, angleDiff, windowPrior);
predictedFractionCW(predictedFractionCW < 0.0001) = 0.0001;
predictedFractionCW(predictedFractionCW > 0.9999) = 0.9999;

% Accumulate log likelihood for the entire data set by adding up the
% log likelihoods for each trial.
FractionCW = percentCW/100;
numberCorrect = round(nTrials.*FractionCW);
numberNotCorrect = nTrials-numberCorrect;
negLogLH = -sum(sum(numberCorrect .* log(predictedFractionCW) + ...
                    numberNotCorrect .* log(1-predictedFractionCW)));
                
% Plot the optimization progress
negLLH = [negLLH negLogLH];
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
xlabel('Iteration', 'FontSize', 25)
ylabel('Negative Log LH', 'FontSize', 25)
drawnow

disp(['-logLH:' num2str(round(negLogLH)) ' '  ...
         ', Params: ' num2str([roundn(stdNoise, -2) roundn([lapseRate guessRate], -4) roundn([priorRange smoothFactor],-2)])])

end

function change_current_figure(h)
set(0,'CurrentFigure',h)
end