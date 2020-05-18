%%%%%%%%%%% Simulate the estimate for different sensory measurement
th = -80:0.01:80;
stdTotal = sqrt([6.0243    8.0650].^2 + 5.3989^2);
m = (-30:5:0)';
nTrialPerM = 500;
meanEstimate = NaN(length(stdTotal), length(m));
semEstimate = NaN(length(stdTotal), length(m));

% Prior                
pTheta = ones(size(th));
pTheta(th<0 | th>32) = 0;
pTheta = pTheta ./ sum(pTheta);

for kk = 1 : length(stdTotal)
    % Resampling distribution centered on m
    pResample = normpdf(repmat(th, length(m), 1), ...
                        repmat(m, 1, length(th)), stdTotal(kk)); 
    pResample(:,th < 0) = 0;
    
    for ii = 1 : length(m)
        mResampled = randpdf(pResample(ii, :), th, [nTrialPerM 1]);
        mResampled(isnan(mResampled)) = [];
        pMr_theta = normpdf(repmat(th, length(mResampled), 1), ...
                        repmat(mResampled, 1, length(th)), stdTotal(kk)); 
        pTheta_mr = bsxfun(@times,  pMr_theta, pTheta);
        pTheta_mr = bsxfun(@rdivide, pTheta_mr, sum(pTheta_mr, 2));
        tempEstimate = th * pTheta_mr';
        meanEstimate(kk, ii) = mean(tempEstimate);
        semEstimate(kk, ii) = std(bootstrp(1000, @(x) mean(x), tempEstimate));
    end
end

figure;
hold on
% plot(m, meanEstimate, 'o-')
for kk = 1 : length(stdTotal)
    errorbar(m, meanEstimate(kk, :), semEstimate(kk, :))
end
xlabel('Memory sample (deg)')
ylabel('Estimate (deg)')
legend('Low noise', 'High noise', 'Location', 'NorthWest')