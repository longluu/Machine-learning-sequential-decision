stdNoiseLevel = [3 18];
nStimOrientation = 13;
nTrialPerCondition = 70;
nTrialTotal = length(stdNoiseLevel) * nStimOrientation * nTrialPerCondition;
sampleSize = 24;
thetaStimulus = NaN(nTrialTotal, sampleSize);
counter = 1;
for ii = 1 : nStimOrientation
    for jj = 1 : length(stdNoiseLevel)
        thetaStimulus(counter:counter+nTrialPerCondition-1, :) = normrnd(0, stdNoiseLevel(jj), [nTrialPerCondition, sampleSize]);
        counter = counter + nTrialPerCondition;
    end
end

figure; 
subplot(2, 1, 1)
plot(1:nTrialTotal, thetaStimulus, 'o')
subplot(2, 1, 2)
plot(1:nTrialTotal, mean(thetaStimulus, 2), '-')
save('thetaStimulus', 'thetaStimulus')