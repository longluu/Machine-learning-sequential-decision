function thetaMedian = medianPdf(theta, pTheta)
pTheta(isnan(theta)) = [];
theta(isnan(theta)) = [];
pThetaCum = cumtrapz(theta, pTheta);
[~, indexMean] = min(abs(pThetaCum - 0.5*max(pThetaCum)));
thetaMedian = theta(indexMean);
end