function y = priorFunction(theta, priorName, params)
if strcmp(priorName, 'UniGauss')
    % params = mean, std, range
    meanGauss = params(1);
    stdGauss = params(2);
    range = params(3);
    y = normpdf(theta, meanGauss, stdGauss);
    if strcmp(choice, 'CCW')
        y(theta>0 | theta<-range) = 0;
    else
        y(theta<0 | theta>range) = 0;
    end
    
    % Normalize to make a probability distribution
    scalingFactor = trapz(theta, y);
    y = y / scalingFactor;
elseif strcmp(priorName, 'TukeyWindow')
    % params = range, fractionTaper, cutoffPoint, mean
    range = [0 params(1)];
    fractionTaper = params(2);
    cutoffPoint = params(3);
    y = TukeyWindowTruncate(range, fractionTaper, theta, cutoffPoint);
    meanY = trapz(theta, theta.*y);
    stepShift = round((params(4) - meanY)/ diff(theta(1:2)));
    y = circshift(y, [1,stepShift]);
elseif strcmp(priorName, 'Gamma')
    % params = a, b, range (mean = a*b, var = a*b^2)
    y = zeros(size(theta));
    indexSelect = theta>0 & theta<=params(3);
    y(indexSelect) = gampdf(theta(indexSelect), params(1), params(2));
    
    if strcmp(choice, 'CCW')
        stepShift = round(params(3)/diff(theta(1:2)));
        y = circshift(y, [1,-stepShift]);
    end
    
    % Normalize to make a probability distribution
    scalingFactor = trapz(theta, y);
    y = y / scalingFactor;    
end