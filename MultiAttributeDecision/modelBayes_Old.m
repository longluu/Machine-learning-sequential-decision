%%%%%%%% Sequential multi-attribute judgments %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
% 

clc; clear

%% Set the parameters of the model
nTrials = 4000;
thetaNoise = 3;
gammaNoise = 6;
    
% prior on orientation decision
pCategory = [0.5 0.5];

% conditional prior on orientation
thetaLength = 1000;
theta = linspace(-50,50,thetaLength);
priorName = 'TukeyWindow';
parameterCCW = [21 0 0 -21/2]; % range, fractionTaper, cutoffPoint, mean
parameterCW = [21 0 0 21/2]; % range, fractionTaper, cutoffPoint, mean
pTheta_CCW = priorFunction(theta, priorName, parameterCCW);
pTheta_CW = priorFunction(theta, priorName, parameterCW);

% conditional prior on numerosity
gammaLength = 1000;
gamma = linspace(-40, 40,gammaLength);
priorName = 'TukeyWindow';
parameterCCW = [20 0 0 -10]; % range, fractionTaper, cutoffPoint, mean
parameterCW = [20 0 0 10]; % range, fractionTaper, cutoffPoint, mean
pGamma_CCW = priorFunction(gamma, priorName, parameterCCW);
pGamma_CW = priorFunction(gamma, priorName, parameterCW);

%% Simulate the experiment
% Generate the stimulus
categoryTrue = [zeros(1, nTrials/2) ones(1, nTrials/2)]; % 1:CCW, 2:CW
indexSampleCCW = randp(pTheta_CCW, 1, ceil(nTrials/2));
indexSampleCW = randp(pTheta_CW, 1, nTrials-ceil(nTrials/2));
thetaTrue = [theta(indexSampleCCW) theta(indexSampleCW)]';
indexSampleCCW = randp(pGamma_CCW, 1, ceil(nTrials/2));
indexSampleCW = randp(pGamma_CW, 1, nTrials-ceil(nTrials/2));
gammaTrue = [gamma(indexSampleCCW) gamma(indexSampleCW)]';

% Generate the measurements
mTheta = normrnd(thetaTrue, thetaNoise);
stdGamma = gammaNoise*ones(nTrials, 1);
mGamma = normrnd(gammaTrue, stdGamma);

%% Standard Bayesian
integrandCCW = repmat(pTheta_CCW, nTrials, 1) ...
    .* normpdf(repmat(theta, nTrials, 1), repmat(mTheta, 1, thetaLength), thetaNoise);
wCCW = trapz(theta, integrandCCW, 2);
integrandCW = repmat(pTheta_CW, nTrials, 1) ...
    .* normpdf(repmat(theta, nTrials, 1), repmat(mTheta, 1, thetaLength), thetaNoise);
wCW = trapz(theta, integrandCW, 2);
priorCCW = repmat(wCCW*pCategory(1), 1, gammaLength) .* repmat(pGamma_CCW, nTrials, 1);
priorCW = repmat(wCW*pCategory(2), 1, gammaLength) .* repmat(pGamma_CW, nTrials, 1);
priorGamma = priorCCW + priorCW;
likelihoodGamma = normpdf(repmat(gamma, nTrials, 1), repmat(mGamma, 1, gammaLength), repmat(stdGamma, 1, gammaLength));
pGamma_m = likelihoodGamma .* priorGamma;
scaleConstant = trapz(gamma, pGamma_m, 2);
pGamma_m = pGamma_m ./ repmat(scaleConstant, 1, gammaLength);
gammaEstimate = trapz(gamma, repmat(gamma, nTrials, 1) .* pGamma_m, 2);

% Scatter plot of the result
gammaEstAve = NaN(1, 55);
imagePlot = zeros(55);
for ii = 1 : length(imagePlot)
    indexSelect = gammaTrue>=ii+4 & gammaTrue<ii+5;
    gammaEstAve(ii) = nanmean(gammaEstimate(indexSelect));
    for jj = 1 : length(imagePlot)
        imagePlot(jj, ii) = sum(gammaEstimate(indexSelect)>=jj+4 & gammaEstimate(indexSelect)<jj+5);
    end
end
imagePlot = flipud(imagePlot);
myfilter = fspecial('gaussian', [1 1], 1);
smoothImage = imfilter(imagePlot, myfilter, 'replicate');
smoothImage = imresize(smoothImage, 2);

hScatter=figure;
figPos = [0.3, 0.2, 0.5, 0.7];
imshow(smoothImage, [min(imagePlot(:)) max(imagePlot(:))])
hold on
set(gca, 'FontSize',20)
plot([length(smoothImage) 1], [1 length(smoothImage)], '--b', 'LineWidth', 1.2)
xlabel('True numerosity')
ylabel('Numerosity estimate')
axis on
set(gca, 'XTick', round(linspace(1,length(smoothImage),5)), 'XTickLabel', num2cell(round(linspace(5,60,5))), ...
    'YTick', round(linspace(1,length(smoothImage),5)), 'YTickLabel', num2cell(round(linspace(60,5,5))));
set(hScatter,'Units','normalized','Position',figPos)
title(['Noise = ' num2str(thetaNoise)])
set(hScatter,'Units','normalized','Position',figPos)
tightfig

hMean = figure;
hold on
set(gca, 'FontSize',20)
plot(5.5:1:59.5, gammaEstAve, 'LineWidth', 2)
plot(5.5:1:59.5, 5.5:1:59.5, '--k', 'LineWidth', 1.5)
xlim([5 60])
ylim([5 60])
axis square
xlabel('True numerosity')
ylabel('Numerosity estimate')

%% Self-consistent Bayesian
% Decision task
integrandCCW = repmat(pTheta_CCW, nTrials, 1) ...
    .* normpdf(repmat(theta, nTrials, 1), repmat(mTheta, 1, thetaLength), thetaNoise);
wCCW = trapz(theta, integrandCCW, 2);
integrandCW = repmat(pTheta_CW, nTrials, 1) ...
    .* normpdf(repmat(theta, nTrials, 1), repmat(mTheta, 1, thetaLength), thetaNoise);
wCW = trapz(theta, integrandCW, 2);
integrandCCW = repmat(pGamma_CCW, nTrials, 1) ...
    .* normpdf(repmat(gamma, nTrials, 1), repmat(mGamma, 1, gammaLength), gammaNoise(1));
uCCW = trapz(gamma, integrandCCW, 2);
integrandCW = repmat(pGamma_CW, nTrials, 1) ...
    .* normpdf(repmat(gamma, nTrials, 1), repmat(mGamma, 1, gammaLength), gammaNoise(2));
uCW = trapz(gamma, integrandCW, 2);
pCCW_mGamma_mTheta = pCategory(1) * wCCW .* uCCW;
pCW_mGamma_mTheta = pCategory(2) * wCW .* uCW;
decisionCW = pCW_mGamma_mTheta > pCCW_mGamma_mTheta;

% Estimation task
priorCCW = repmat(pGamma_CCW, nTrials, 1);
priorCW = repmat(pGamma_CW, nTrials, 1);
priorGamma = NaN(size(priorCW));
priorGamma(decisionCW==0, :) = priorCCW(decisionCW==0, :);
priorGamma(decisionCW==1, :) = priorCW(decisionCW==1, :);
likelihoodGamma = normpdf(repmat(gamma, nTrials, 1), repmat(mGamma, 1, gammaLength), repmat(stdGamma, 1, gammaLength));
pGamma_m = likelihoodGamma .* priorGamma;
scaleConstant = trapz(gamma, pGamma_m, 2);
pGamma_m = pGamma_m ./ repmat(scaleConstant, 1, gammaLength);
gammaEstimate = trapz(gamma, repmat(gamma, nTrials, 1) .* pGamma_m, 2);

% Plot the result
figure
gammaEstAveCCW = NaN(1, 55);
gammaEstAveCW = NaN(1, 55);
imagePlot = zeros(55);
for ii = 1 : length(imagePlot)
    indexSelect = gammaTrue>=ii+4 & gammaTrue<ii+5;
    indexCCW = indexSelect & ~decisionCW;
    gammaEstAveCCW(ii) = nanmean(gammaEstimate(indexCCW));
    indexCW = indexSelect & decisionCW;
    gammaEstAveCW(ii) = nanmean(gammaEstimate(indexCW));    
    for jj = 1 : length(imagePlot)
        imagePlot(jj, ii) = sum(gammaEstimate(indexSelect)>=jj+4 & gammaEstimate(indexSelect)<jj+5);
    end
end
imagePlot = flipud(imagePlot);
myfilter = fspecial('gaussian', [1 1], 1);
smoothImage = imfilter(imagePlot, myfilter, 'replicate');
smoothImage = imresize(smoothImage, 2);

hScatter=figure;
figPos = [0.3, 0.2, 0.5, 0.7];
imshow(smoothImage, [min(imagePlot(:)) max(imagePlot(:))])
hold on
set(gca, 'FontSize',20)
plot([length(smoothImage) 1], [1 length(smoothImage)], '--b', 'LineWidth', 1.2)
xlabel('True numerosity')
ylabel('Numerosity estimate')
axis on
set(gca, 'XTick', round(linspace(1,length(smoothImage),5)), 'XTickLabel', num2cell(round(linspace(5,60,5))), ...
    'YTick', round(linspace(1,length(smoothImage),5)), 'YTickLabel', num2cell(round(linspace(60,5,5))));
set(hScatter,'Units','normalized','Position',figPos)
tightfig

hMean = figure;
hold on
set(gca, 'FontSize',20)
plot(5.5:1:59.5, gammaEstAveCCW, 'b-', 'LineWidth', 2)
plot(5.5:1:59.5, gammaEstAveCW, 'r-', 'LineWidth', 2)
plot(5.5:1:59.5, 5.5:1:59.5, '--k', 'LineWidth', 1.5)
xlim([5 60])
ylim([5 60])
axis square
xlabel('True numerosity')
ylabel('Numerosity estimate')
