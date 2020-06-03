function [predictFractionCW] = BayesFractionCW(priorRange, smoothFactor, stdNoiseLevel, lapseRate, angleDiff, windowPrior)
%% Compute p(CW_hat|theta) using integral(p(m|theta), 0, Inf)
% NOTE that this code is NOT general and confined to the case of a
% symmetric prior around the reference
try
    mLength = 2000;
    thetaRepmat = repmat(angleDiff, mLength, 1);
    numAngleDiff = length(angleDiff);
    predictFractionCW = NaN(length(stdNoiseLevel), numAngleDiff);
    for ii = 1 : length(stdNoiseLevel) 
        % Generate grid of m
        mOriginal = linspace(-80, 80, mLength);
        mRepmat = repmat(mOriginal', 1, length(angleDiff));
        pM_Theta = normpdf(mRepmat, thetaRepmat, stdNoiseLevel(ii));

        % Calculate the percent CW p(CW_hat|theta) = integral(p(m|theta), 0, Inf)
        indexIntegrate = mOriginal > 0;
        predictFractionCW(ii, :) = lapseRate + (1-2*lapseRate) * trapz(mOriginal(indexIntegrate), pM_Theta(indexIntegrate, :), 1);
    end  
    
%% Compute p(CW_hat|theta) using integral(p(thetaHat|theta), 0, Inf)
%     cutoffPoint = priorRange;
%     mLength = 1000;
%     if windowPrior == 1 
%         theta = linspace(-50, 50, 1000);  
%         pTheta_CCW = TukeyWindow([-cutoffPoint 0], 1, smoothFactor, theta);
%         pTheta_CW = TukeyWindow([0 cutoffPoint], 0, smoothFactor, theta);
%     else
%         theta = linspace(-75, 75, 1500);
%         pTheta_CCW = HalfGaussWindow([-cutoffPoint 0], 1, smoothFactor, theta);
%         pTheta_CW = HalfGaussWindow([0 cutoffPoint], 0, smoothFactor, theta);
%     end 
%     thetaRepmat = repmat(theta, mLength, 1);
%     pTheta = pTheta_CCW + pTheta_CW;
%     pTheta = repmat(pTheta, mLength, 1);
%     numAngleDiff = length(angleDiff);
%     rangeCollapse = round(numAngleDiff/2);
%     estimateModel.Xval = cell(length(stdNoiseLevel), rangeCollapse);
%     estimateModel.Yval = cell(length(stdNoiseLevel), rangeCollapse);
% 
%     for jj = 1 : length(stdNoiseLevel) 
%         % Generate grid of m
%         mOriginal = linspace(-80, 80, mLength);
%         mRepmat = repmat(mOriginal', 1, length(theta));
%         pM_Theta = normpdf(thetaRepmat, mRepmat, stdNoiseLevel(jj));
%         pM = trapz(theta, pTheta.*pM_Theta, 2);
% 
%         % Calculate the BLS estimation of theta
%         fM = (1./pM) .* trapz(theta, thetaRepmat.*pTheta.*pM_Theta, 2);
% 
%         for kk = 1 : rangeCollapse         
%             % Calculate the percept
%             m = mOriginal;
%             fM_new = fM;
%             m(isnan(fM_new)) = [];
%             fM_new(isnan(fM_new)) = [];
%             m(isinf(fM_new)) = [];
%             fM_new(isinf(fM_new)) = [];                    
%             thetaHatResampled = linspace(min(fM_new),max(fM_new), mLength);
%             thetaHatInverse = interp1(fM_new, m, thetaHatResampled,'pchip');
%             D_thetaHatInverse = (thetaHatInverse(3:end) - thetaHatInverse(1:end-2)) ./ (thetaHatResampled(3:end) - thetaHatResampled(1:end-2));
%             D_thetaHatInverse = [D_thetaHatInverse(1) D_thetaHatInverse D_thetaHatInverse(end)];                    
%             pThetaHat_Theta = normpdf(thetaHatInverse, angleDiff(kk), stdNoiseLevel(jj)) .* abs(D_thetaHatInverse);
%             normalizedFactor = trapz(thetaHatResampled, pThetaHat_Theta);
%             pThetaHat_Theta = pThetaHat_Theta/normalizedFactor;        
%             estimateModel.Xval{jj, kk} = thetaHatResampled;
%             estimateModel.Yval{jj, kk} = pThetaHat_Theta; 
%         end
%     end  
% 
%     % Calculate the percent CW
%     predictFractionCW = NaN(length(stdNoiseLevel), numAngleDiff);
%     for ii = 1 : length(stdNoiseLevel)
%         for jj = 1 : numAngleDiff
%             if angleDiff(jj) <= 0
%                 thetaHat = estimateModel.Xval{ii, jj};
%                 pThetaHat = estimateModel.Yval{ii, jj};
%                 indexIntegrate = thetaHat > 0;
%                 predictFractionCW(ii,jj) = lapseRate + (1-2*lapseRate)*trapz(thetaHat(indexIntegrate), pThetaHat(indexIntegrate));
%             else
%                 thetaHat = estimateModel.Xval{ii, numAngleDiff-jj+1};
%                 pThetaHat = estimateModel.Yval{ii, numAngleDiff-jj+1};
%                 indexIntegrate = thetaHat < 0;
%                 predictFractionCW(ii,jj) = lapseRate + (1-2*lapseRate)*trapz(thetaHat(indexIntegrate), pThetaHat(indexIntegrate));
%             end            
%         end
%     end
catch e
    keyboard
    rethrow(e)
end
end