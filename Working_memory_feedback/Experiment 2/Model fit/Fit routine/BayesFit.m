function [fitParameter, negLogLH, psychCurveExp2] = BayesFit(percentCW,nTrials,angleDiff,plotFigure,stimDuration,colorName,hLegend, fixedPSE)
% We optimize two parameters: the detection threshold and its standard
% deviation of noise because we use the normal curve for fitting
minStd = 1;
maxStd = 20;
if fixedPSE
    vlb = [mean(angleDiff) minStd minStd 0 ];
    vub = [mean(angleDiff) maxStd maxStd 0.2];
else
    vlb = [-50 10   minStd minStd 0];
    vub = [-5  60   maxStd maxStd 0.2];
end
thresholdInitial = [-12 30];
stdNoiseInitial = [2 4] + 2*(rand(1,length(stimDuration))-0.5);
lapseRateInitial = 0;
searchParameter0 = [thresholdInitial stdNoiseInitial lapseRateInitial];

% Optimization begins
% % using fmincon
% options = optimset('fmincon');
% options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
% [fitParameter, negLogLH] = fmincon(@(searchParameter)normalFitErrorFunction(searchParameter(1),searchParameter(2:4),...
%         searchParameter(5),searchParameter(6),angleDiff,percentCW,nTrials),searchParameter0,[],[],[],[],vlb,vub,[],options);

% using fminsearchbnd
searchOptions = struct('Display','none', 'TolX',1e-6, 'FunValCheck','off');
optimFun = @(searchParameter, varargin) normalFitErrorFunction(searchParameter(1:2),searchParameter(3:4),...
        searchParameter(5),angleDiff,percentCW,nTrials);
[fitParameter, negLogLH] = fminsearchbnd(optimFun, searchParameter0, vlb, vub, searchOptions,[]);

% Plot the fit curve and data
angleDiffResampled = linspace(angleDiff(1),angleDiff(end),1000);
fitpercentCW = 100*BayesFractionCW(fitParameter(1), fitParameter(2:3), fitParameter(4), angleDiffResampled);
figure(plotFigure)
legendName = cell(1,length(stdNoiseInitial));
lineWidth = 3;
fontSize = 20;
set(gca,'FontSize',15)
hold on
for ii = 1 : length(stdNoiseInitial)
    hLegend(ii) = plot(angleDiffResampled,fitpercentCW(ii,:),'Color',...
        colorName{ii},'LineWidth',3);
    plot(angleDiff, percentCW(ii,:), [colorName{ii} 'o'], 'MarkerSize',13,...
        'MarkerFaceColor', colorName{ii}, 'LineWidth', lineWidth);
    legendName{ii} = ['Noise : ' num2str(roundn(fitParameter(ii+1),-2))];
end
xlim([min(angleDiff) max(angleDiff)])
xlabel('True numerosity','FontSize',fontSize)
ylabel('Percent greater','FontSize',fontSize)
legend(hLegend, legendName, 'Location', 'SouthEast')

% Print some info
fprintf('PSE = %4.2f, ',fitParameter(1));
fprintf('sigma = %4.2f %4.2f %4.2f \n',fitParameter(2:3));
fprintf('Lapse rate: %4.3f', fitParameter(4));

% Save the psychometric curve
psychCurveExp2.angleDiff = angleDiffResampled;
psychCurveExp2.percentCW = fitpercentCW;
end

function negLogLH = normalFitErrorFunction(priorRange, stdNoise, lapseRate, angleDiff, percentCW, nTrials)
% Compute predicted Percent Correct for each data point
priorRange = real(priorRange);
epsilon = 10^(-20);
predictedFractionCW = BayesFractionCW(priorRange,stdNoise, lapseRate, angleDiff);
predictedFractionCW(predictedFractionCW < epsilon) = epsilon;
predictedFractionCW(predictedFractionCW > 1 - epsilon) = 1 - epsilon;

% Accumulate log likelihood for the entire data set by adding up the
% log likelihoods for each trial.
FractionCW = percentCW/100;
numberCorrect = round(nTrials.*FractionCW);
numberNotCorrect = nTrials-numberCorrect;
negLogLH = -sum(sum(numberCorrect .* log(predictedFractionCW) + ...
                    numberNotCorrect .* log(1-predictedFractionCW)));

end

function predictFractionCW = BayesFractionCW(priorRange, stdNoise, lapseRate, angleDiff)
predictFractionCW = NaN(length(stdNoise), length(angleDiff));

pC = [0.5, 0.5]';
pthcw = priorRange(2);
pthccw = priorRange(1); % paramsAll(4)

rangeth = [-60 60];
th = rangeth(1):dstep:rangeth(2);
th = round(th, -log10(dstep));
nth = length(th);

pthGC = zeros(2,nth);
pth = (TukeyWindow([0 pthcw], 0, smoothFactor, th) + TukeyWindow([pthccw 0], 1, smoothFactor, th))/2;
pth(th==0) = 0;
pth(th==0) = max(pth);
pthGC(1,:) = pth;
pthGC(2,:) = pth;

for kk = 1 : length(stdNoise)
    rangeM = [min(angleDiff)-5*stdSensory(kk) max(angleDiff)+5*stdSensory(kk)];
    if rangeM(2) < rangeth(2)
        rangeM = rangeth;
    end
    nm = 1000;
    m = linspace(rangeM(1), rangeM(2), nm);

    M = repmat(m',1,nth);
    THm = repmat(th, nm, 1); 
    
    pmGth = exp(-((M-THm).^2)./(2*stdNoise(kk)^2));
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
    PChGtheta = PChGm * pmGth(:, ismember(th, angleDiff));
    PChGtheta_lapse = lapseRate + (1 - 2*lapseRate) * PChGtheta;
    PChGtheta_lapse = PChGtheta_lapse ./ repmat(sum(PChGtheta_lapse, 1), 2, 1);
    
end
end