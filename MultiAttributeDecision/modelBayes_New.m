%%%%%%%% Sequential multi-attribute judgments %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
% 
%% Set up parameters
nTrials = 10000;
thetaRange = 20;
gammaRange = 20;
thetaNoise = 1;
gammaNoise = 7;
pC = [0.5 0.5]; % [pCW pCCW]

% conditional prior on orientation
thetaLength = 1000;
theta = linspace(-30,30,thetaLength);
priorName = 'TukeyWindow';
parameterCCW = [thetaRange 0 0 -thetaRange/2]; % range, fractionTaper, cutoffPoint, mean
parameterCW = [thetaRange 0 0 thetaRange/2]; % range, fractionTaper, cutoffPoint, mean
pTheta_CCW = priorFunction(theta, priorName, parameterCCW);
pTheta_CW = priorFunction(theta, priorName, parameterCW);

% conditional prior on numerosity
gammaLength = 1000;
gamma = linspace(-40, 40,gammaLength);
priorName = 'TukeyWindow';
parameterCCW = [gammaRange 0 0 -gammaRange/2]; % range, fractionTaper, cutoffPoint, mean
parameterCW = [gammaRange 0 0 gammaRange/2]; % range, fractionTaper, cutoffPoint, mean
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
Mth = repmat(mTheta, 1, thetaLength);
ThM = repmat(theta, nTrials, 1);
mGamma = normrnd(gammaTrue, gammaNoise);
Mga = repmat(mGamma, 1, gammaLength);
GaM = repmat(gamma, nTrials, 1);
Chat = -ones(1, nTrials);

%% Standard Bayesian
% Decision task
pMtheta_theta = exp(-((Mth-ThM).^2)./(2*thetaNoise^2));
pMgamma_gamma = exp(-((Mga-GaM).^2)./(2*gammaNoise^2));
pMtheta_C = pMtheta_theta * [pTheta_CW' pTheta_CCW'];
pMgamma_C = pMgamma_gamma * [pGamma_CW' pGamma_CCW'];
pC_MthMga = pMtheta_C .* pMgamma_C .* repmat(pC, nTrials, 1);
Chat(pC_MthMga(:, 1) > pC_MthMga(:, 2)) = 1; 

% Estimation task
pC_gamma = [pGamma_CW; pGamma_CCW] .* repmat(pC', 1, gammaLength);
pMtheta_gamma = pMtheta_C * pC_gamma;
pGamma_Mth_Mga =  pMtheta_gamma .* pMgamma_gamma;
pGamma_Mth_Mga = bsxfun(@rdivide, pGamma_Mth_Mga, sum(pGamma_Mth_Mga, 2));
gammaHat = gamma * pGamma_Mth_Mga';

% Plot the estimate
figure
hold on
plot(gammaTrue, gammaHat, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none')
plot([-gammaRange gammaRange], [0 0], '--b')
plot([0 0], [-gammaRange gammaRange], '--b')
plot([-gammaRange gammaRange], [-gammaRange gammaRange], '--b')
axis('square')
xlabel('True attractiveness')
ylabel('Estimated attractiveness')

%% Self-consistent Bayesian
% Estimation task
pGamma_Chat = NaN(size(pMgamma_gamma));
pGamma_Chat(Chat == 1, :) = repmat(pGamma_CW, sum(Chat == 1), 1);
pGamma_Chat(Chat == -1, :) = repmat(pGamma_CCW, sum(Chat == -1), 1);
pGamma_Mth_Chat = pMgamma_gamma .* pGamma_Chat;
pGamma_Mth_Chat = bsxfun(@rdivide, pGamma_Mth_Chat, sum(pGamma_Mth_Chat, 2));
gammaHat = gamma * pGamma_Mth_Chat';

% Plot the estimate
figure
hold on
plot(gammaTrue, gammaHat, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none')
plot([-gammaRange gammaRange], [0 0], '--b')
plot([0 0], [-gammaRange gammaRange], '--b')
plot([-gammaRange gammaRange], [-gammaRange gammaRange], '--b')
axis('square')
xlabel('True attractiveness')
ylabel('Estimated attractiveness')









