function [y, pH] = priorFunction(theta, priorName, choice, params)
if strcmp(priorName, 'UniGauss')
    % params = mean, std, range
    meanGauss = params(1);
    stdGauss = params(2);
    range = params(3);
    y = normpdf(theta, meanGauss, stdGauss);
    pH(1) = trapz(theta(theta<0), y(theta<0));
    pH(2) = trapz(theta(theta>0), y(theta>0));
    if strcmp(choice, 'CCW')
        y(theta>0 | theta<-range) = 0;
    else
        y(theta<0 | theta>range) = 0;
    end
    
    % Normalize to make a probability distribution
    scalingFactor = trapz(theta, y);
    y = y / scalingFactor;
elseif strcmp(priorName, 'TukeyWindow')
    % params = range, fractionTaper, cutoffPoint
    if strcmp(choice, 'CCW')
        range = [-params(1) 0];
    else
        range = [0 params(1)];
    end
    fractionTaper = params(2);
    cutoffPoint = params(3);
    y = TukeyWindowNew(range, fractionTaper, theta, cutoffPoint);
end