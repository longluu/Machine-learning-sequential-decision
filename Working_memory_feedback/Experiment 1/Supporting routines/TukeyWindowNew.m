function y = TukeyWindowNew(range, fractionTaper, theta, cutoffPoint)
% Find the end points of Tukey window in true coordinate
if sum(range) > 0
    decision = 1;
    A = range(1) + cutoffPoint;
    B = range(2);
else
    decision = -1;
    A = range(1);
    B = range(2) - cutoffPoint;    
end
r = fractionTaper;
a = r/4;
b = 1 - a;
endPointNorm = [0 1];
endPointTrue = (B-A) * (endPointNorm - a) / (b-a) + A;

% Create the Tukey window
y = zeros(size(theta));
nPoints = sum(theta >= endPointTrue(1) & theta <= endPointTrue(2));
y(theta >= endPointTrue(1) & theta <= endPointTrue(2)) = tukeywin(nPoints, r);

% Create the sharp edge on one side of the window
switch decision
    case -1
        midPoint = (A - B)/2;
        y(theta > -cutoffPoint) = 0;
        y(theta >= midPoint & theta < -cutoffPoint) = max(y);
    case 1
        midPoint = (B - A)/2;        
        y(theta < cutoffPoint) = 0;
        y(theta >= cutoffPoint & theta < midPoint) = max(y);
    case 'NA'
end
    
% Normalize to make a probability distribution
scalingFactor = trapz(theta, y);
y = y / scalingFactor;
end