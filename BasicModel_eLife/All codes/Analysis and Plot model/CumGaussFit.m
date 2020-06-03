function [fitParameter, negLogLH, psychCurveExp2] = CumGaussFit(percentCW,nTrials,angleDiff,plotFigure,stdNoiseLevel,colorName, fixedPSE)
% We optimize two parameters: the detection threshold and its standard
% deviation of noise because we use the normal curve for fitting
minStd = 0.01;
maxStd = 20;
if fixedPSE
    vlb = [0 minStd minStd minStd 0 ];
    vub = [0 maxStd maxStd maxStd 0.2];
else
    vlb = [angleDiff(1)   angleDiff(1)   angleDiff(1)   minStd minStd minStd 0];
    vub = [angleDiff(end) angleDiff(end) angleDiff(end) maxStd maxStd maxStd 0.2];
end
thresholdInitial = zeros(1, length(stdNoiseLevel));
stdNoiseInitial = [1:length(stdNoiseLevel)] + 2*(rand(1,length(stdNoiseLevel))-0.5);
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
optimFun = @(searchParameter, varargin) normalFitErrorFunction(searchParameter(1:3),searchParameter(4:6),...
        searchParameter(7),angleDiff,percentCW,nTrials);
[fitParameter, negLogLH] = fminsearchbnd(optimFun, searchParameter0, vlb, vub, searchOptions,[]);

% Plot the fit curve and data
angleDiffResampled = linspace(angleDiff(1),angleDiff(end),1000);
fitpercentCW = 100*normalFractionCW(fitParameter(1:3), fitParameter(4:6), fitParameter(7), angleDiffResampled);
figure(plotFigure)
legendName = cell(1,length(stdNoiseInitial));
lineWidth = 3;
fontSize = 20;
set(gca,'FontSize',15)
hold on
hLegend = NaN(1, length(stdNoiseInitial));
for ii = 1 : length(stdNoiseInitial)
    hLegend(ii) = plot(angleDiffResampled,fitpercentCW(ii,:),'Color',...
        colorName{ii},'LineWidth',3);
    plot(angleDiff, percentCW(ii,:), [colorName{ii} 'o'], 'MarkerSize',13,...
        'MarkerFaceColor', colorName{ii}, 'MarkerEdgeColor', 'none', 'LineWidth', lineWidth);
    legendName{ii} = ['noise = ' num2str(roundn(fitParameter(ii+3),-1))];
end
xlim([min(angleDiff) max(angleDiff)])
xlabel('True angle (degree)','FontSize',fontSize)
ylabel('Percent CW','FontSize',fontSize)
legend(hLegend, legendName, 'Location', 'SouthEast')

% Print some info
fprintf('PSE = %4.2f %4.2f %4.2f, ',fitParameter(1:3));
fprintf('sigma = %4.2f %4.2f %4.2f \n',fitParameter(4:6));
fprintf('Lapse rate: %4.2f', fitParameter(7));

% Save the psychometric curve
psychCurveExp2.angleDiff = angleDiffResampled;
psychCurveExp2.percentCW = fitpercentCW;
end

function negLogLH = normalFitErrorFunction(detectThreshold, stdNoise, lapseRate, angleDiff, percentCW, nTrials)
% Compute predicted Percent Correct for each data point
detectThreshold = real(detectThreshold);
predictedFractionCW = normalFractionCW(detectThreshold,stdNoise, lapseRate, angleDiff);
predictedFractionCW(predictedFractionCW < 0.0001) = 0.0001;
predictedFractionCW(predictedFractionCW > 0.9999) = 0.9999;

% Accumulate log likelihood for the entire data set by adding up the
% log likelihoods for each trial.
FractionCW = percentCW/100;
numberCorrect = round(nTrials.*FractionCW);
numberNotCorrect = nTrials-numberCorrect;
negLogLH = -sum(sum(numberCorrect .* log(predictedFractionCW) + ...
                    numberNotCorrect .* log(1-predictedFractionCW)));

end

function predictFractionCW = normalFractionCW(detectThreshold, stdNoise, lapseRate, angleDiff)
predictFractionCW = NaN(length(stdNoise), length(angleDiff));
for ii = 1 : length(stdNoise)
    predictFractionCW(ii,:) = lapseRate + (1-2*lapseRate)*(1-normcdf(detectThreshold(ii),angleDiff,stdNoise(ii)));
end
end