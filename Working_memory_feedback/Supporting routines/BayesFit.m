function [fitParameter, negLogLH] = BayesFit(percentEllipseCW,nTrials,angleDiff,plotFigure)
% We optimize 4 parameters: 3 std of ellispe coherence and 1 prior range
% (uniform prior)
sigma_ThetaInitial = [7.8 4 2.6];
thetaPriorInitial = 20;
searchParameter0 = [sigma_ThetaInitial thetaPriorInitial];
[fitParameter, negLogLH] = fminsearch(@(searchParameter)BayesFitErrorFunction(searchParameter(1:3), searchParameter(4)...
    ,angleDiff,percentEllipseCW,nTrials),searchParameter0);

% Add the fit curve to the plot in green.
figure(plotFigure)
hold on
angleDiffResampled = linspace(angleDiff(1),angleDiff(end),30);
sigma_theta = fitParameter(1:3);
thetaPrior = fitParameter(4);
fitPercentCW = 100*BayesFractionCW(angleDiffResampled, sigma_theta, thetaPrior);
disp('******* Full Bayes fit *******')
for ii = 1 : length(sigma_ThetaInitial)
    subplot(1,3,ii);
    hold on
    plot(angleDiffResampled,fitPercentCW(ii,:),'red','LineWidth',2);
    plot(angleDiff,percentEllipseCW(ii,:),'bo','MarkerFaceColor','b', 'MarkerSize', 7);
    ylim([0 110])

    % Print some info
    [dummy, ind75] = min(abs(fitPercentCW(ii,:)-75));
    [dummy, ind50] = min(abs(fitPercentCW(ii,:)-50));
    deltaThreshold = angleDiffResampled(ind75)-angleDiffResampled(ind50);
    fprintf(['Contrast ' num2str(ii) ': '])
    fprintf('PSE = %6.4f, ',angleDiffResampled(ind50));
    fprintf('Deltax_thres = %6.4f, ',deltaThreshold);
    fprintf('sigma = %6.4f \n',sigma_theta(ii));
end
end

function f = BayesFitErrorFunction(sigma_theta,thetaPrior,angleDiff,percentEllipseCW,nTrials)
% Compute predicted Percent Correct for each data point
predictedFractionCW = BayesFractionCW(angleDiff, sigma_theta, thetaPrior);
predictedFractionCW(predictedFractionCW <= 0) = 0.0001;
predictedFractionCW(predictedFractionCW >= 1) = 0.9999;

% Accumulate log likelihood for the entire data set by adding up the
% log likelihoods for each trial.
logLikelihood = 0;
FractionCW = percentEllipseCW/100;
numberCW = round(nTrials*FractionCW);
numberCCW = nTrials-numberCW;
for kk = 1 : length(sigma_theta)
    for i = 1:length(angleDiff)
        for j = 1:numberCW(kk,i)
            logLikelihood = logLikelihood + log10(predictedFractionCW(kk,i));
        end
        for j = 1:numberCCW(kk,i)
            logLikelihood = logLikelihood + log10(1-predictedFractionCW(kk,i));
        end
    end
end
f = -logLikelihood;
disp(num2str([f sigma_theta thetaPrior]))
end

function FractionCW = BayesFractionCW(angleDiff, sigma_theta, thetaPrior)
pH = [0.5 0.5];
nSamples = 1000;
FractionCW = NaN(length(sigma_theta), length(angleDiff));

for jj = 1 : length(sigma_theta)    
    % Generate counterbalanced grid of m
    m = linspace(-40,50,nSamples);
    theta1 = linspace(-thetaPrior, 0, nSamples);
    theta2 = linspace(0, thetaPrior, nSamples);

    % Calculate the BLS estimation of theta
    pM =  pH(1) * (normcdf(0, m, sigma_theta(jj)) - normcdf(-thetaPrior, m, sigma_theta(jj)))...
                 + pH(2) * (normcdf(thetaPrior, m, sigma_theta(jj)) - normcdf(0, m, sigma_theta(jj)));
    pM = pM/thetaPrior;
    thetaHat = NaN(1, nSamples);
    for ii = 1 : nSamples
        thetaHat(ii) = (1/(pM(ii)*thetaPrior)) * (pH(1) * trapz(theta1, theta1.*normpdf(theta1, m(ii), sigma_theta(jj))) ...
                        + pH(2) * trapz(theta2, theta2 .* normpdf(theta2, m(ii), sigma_theta(jj))));
    end

    % Calculate the percept
    m(isnan(thetaHat) | isinf(thetaHat)) = [];
    thetaHat(isnan(thetaHat) | isinf(thetaHat))  = [];
    thetaHatResampled = linspace(min(thetaHat),max(thetaHat), nSamples);
    thetaHatInverse = interp1(thetaHat, m, thetaHatResampled,'spline');    
    for kk = 1 : length(angleDiff)  
        pThetaHat_Theta = normpdf(thetaHatInverse, angleDiff(kk), sigma_theta(jj));
        scaling_factor = 1/trapz(thetaHatResampled,pThetaHat_Theta);
        pThetaHat_Theta = pThetaHat_Theta * scaling_factor;
        FractionCW(jj,kk) = trapz(thetaHatResampled(thetaHatResampled>=0), pThetaHat_Theta(thetaHatResampled>=0));
    end
end
end