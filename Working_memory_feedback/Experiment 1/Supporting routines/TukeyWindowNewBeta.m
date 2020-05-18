function y = TukeyWindowNewBeta(range, fractionTaper, theta, condPriorMean)
% Find the end points of Tukey window in true coordinate
r = fractionTaper;
A = range(1);
B = range(2);
a = r/4;
b = 1 - a;
endPointNorm = [0 1];
endPointTrue = (B-A) * (endPointNorm - a) / (b-a) + A;
if sum(range) > 0
    decision = 1;
else
    decision = -1;
end

% Create the Tukey window.
y = zeros(size(theta));
nPoints = sum(theta >= endPointTrue(1) & theta <= endPointTrue(2));
y(theta >= endPointTrue(1) & theta <= endPointTrue(2)) = tukeywin(nPoints, r);

% Create the sharp edge on one side of the window
diffTheta = condPriorMean - abs(sum(range)/2);
nShift = round(diffTheta/(theta(2) - theta(1)));
switch decision
    case -1
        y = circshift(y, [0 -nShift]);
    case 1
        y = circshift(y, [0 nShift]);
    case 'NA'
end
    
% Normalize to make a probability distribution
scalingFactor = trapz(theta, y);
y = y / scalingFactor;
end