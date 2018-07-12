%%%%%%%% The sampling Bayesian model of conditioned perception %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
% 

clc; clear

%% Set the parameters of the model
nTrials = 10000;
stdNoise = 7;
    
% prior on orientation decision
pCategory = [0.5 0.5];

% conditional prior on orientation
thetaLength = 1000;
theta = linspace(-50,50,thetaLength);
priorName = 'TukeyWindow';
thetaRange = 21;
parameterCCW = [thetaRange 0 0 -thetaRange/2]; % range, fractionTaper, cutoffPoint, mean
parameterCW = [thetaRange 0 0 thetaRange/2]; % range, fractionTaper, cutoffPoint, mean
pTheta_CCW = priorFunction(theta, priorName, parameterCCW);
pTheta_CW = priorFunction(theta, priorName, parameterCW);


%% Simulate the experiment
% Generate the stimulus
categoryTrue = [zeros(1, nTrials/2) ones(1, nTrials/2)]; % 1:CCW, 2:CW
indexSampleCCW = randp(pTheta_CCW, 1, ceil(nTrials/2));
indexSampleCW = randp(pTheta_CW, 1, nTrials-ceil(nTrials/2));
thetaTrue = [theta(indexSampleCCW) theta(indexSampleCW)]';
m1 = normrnd(thetaTrue, stdNoise);
m2 = normrnd(thetaTrue, stdNoise);

%% Standard Bayesian
priorCCW = pCategory(1) .* repmat(pTheta_CCW, nTrials, 1);
priorCW = pCategory(2) .* repmat(pTheta_CW, nTrials, 1);
priorTheta = priorCCW + priorCW;
likelihoodTheta = normpdf(repmat(theta, nTrials, 1), repmat(m1, 1, thetaLength), stdNoise) ...
                .* normpdf(repmat(theta, nTrials, 1), repmat(m2, 1, thetaLength), stdNoise);
pTheta_m1_m2 = priorTheta .* likelihoodTheta;
scaleConstant = trapz(theta, pTheta_m1_m2, 2);
pTheta_m1_m2 = pTheta_m1_m2 ./ repmat(scaleConstant, 1, thetaLength);
thetaEstimateSB = trapz(theta, repmat(theta, nTrials, 1) .* pTheta_m1_m2, 2);

%% Self-consistent Bayesian
% Decision task
integrandCCW = repmat(pTheta_CCW, nTrials, 1) ...
    .* normpdf(repmat(theta, nTrials, 1), repmat(m1, 1, thetaLength), stdNoise);
integrandCW = repmat(pTheta_CW, nTrials, 1) ...
    .* normpdf(repmat(theta, nTrials, 1), repmat(m1, 1, thetaLength), stdNoise);
pCCW_m1 = pCategory(1) * trapz(theta, integrandCCW, 2);
pCW_m1 = pCategory(2) * trapz(theta, integrandCW, 2);
decisionCW = pCW_m1 > pCCW_m1;

% Estimation task
priorThetaCond = NaN(size(priorCW));
priorThetaCond(decisionCW==0, :) = priorCCW(decisionCW==0, :);
priorThetaCond(decisionCW==1, :) = priorCW(decisionCW==1, :);
pTheta_m1_m2_C = likelihoodTheta .* priorThetaCond;
scaleConstant = trapz(theta, pTheta_m1_m2_C, 2);
pTheta_m1_m2_C = pTheta_m1_m2_C ./ repmat(scaleConstant, 1, thetaLength);
thetaEstimateSC = trapz(theta, repmat(theta, nTrials, 1) .* pTheta_m1_m2_C, 2);

%% Plot the result
rangeIllustrate = [10 21];

% Standard Bayes
figure
hold on
plot(thetaTrue, thetaEstimateSB, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none')
plot([-thetaRange thetaRange], [0 0], '--b')
plot([0 0], [-thetaRange thetaRange], '--b')
plot([-thetaRange thetaRange], [-thetaRange thetaRange], '--b')
xlim([-thetaRange-1 thetaRange+1])
ylim([-thetaRange-1 thetaRange+1])
axis('square')
xlabel('True attractiveness')
ylabel('Estimated attractiveness')

% Self-consistent Bayes
figure
hold on
plot(thetaTrue, thetaEstimateSC, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none')
plot([-thetaRange thetaRange], [0 0], '--b')
plot([0 0], [-thetaRange thetaRange], '--b')
plot([-thetaRange thetaRange], [-thetaRange thetaRange], '--b')
xlim([-thetaRange-1 thetaRange+1])
ylim([-thetaRange-1 thetaRange+1])
axis('square')
xlabel('True attractiveness')
ylabel('Estimated attractiveness')

% Illustrate the distortion of additional evidence
indInclude = (thetaTrue >= rangeIllustrate(1)) & (thetaTrue <= rangeIllustrate(2) & m1 < 0);
figure
hold on
subplot(4, 1, 1)
hold on
plot(m1(indInclude), 'o-')
plot(zeros(1, sum(indInclude)), '--')
ylim([-27 27])
xlim([0 sum(indInclude)])

subplot(4, 1, 2)
hold on
plot(m2(indInclude), 'o-')
plot(zeros(1, sum(indInclude)), '--')
ylim([-27 27])
xlim([0 sum(indInclude)])

subplot(4, 1, 3)
hold on
plot(thetaEstimateSB(indInclude), 'o-')
plot(zeros(1, sum(indInclude)), '--')
ylim([-27 27])
xlim([0 sum(indInclude)])

subplot(4, 1, 4)
hold on
plot(thetaEstimateSC(indInclude), 'o-')
plot(zeros(1, sum(indInclude)), '--')
ylim([-27 27])
xlim([0 sum(indInclude)])
