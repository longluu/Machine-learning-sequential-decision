function y = TukeyWindowTruncate(range, fractionTaper, theta, cutoffPoint)
% Find the end points of Tukey window in true coordinate
r = fractionTaper;
if sum(range) > 0
    decision = 1;
    A = -range(2);
    B = range(2);
else
    decision = -1;
    A = range(1);
    B = -range(1);
end
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
        y(theta > -cutoffPoint | theta < range(1)) = 0;
    case 1
        y(theta < cutoffPoint | theta > range(2)) = 0;
    case 'NA'
end
    
% Normalize to make a probability distribution
scalingFactor = trapz(theta, y);
y = y / scalingFactor;
end