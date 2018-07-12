function y = HalfGaussWindowNew(range, stdGauss, theta, cutoffPoint)
% Initialize the window 
y = zeros(size(theta));
if sum(range) > 0
    decision = 1;
else
    decision = -1;
end

% Create the window
switch decision
    case -1
        indexGauss = theta <= range(1); 
        y(indexGauss) = normpdf(theta(indexGauss), range(1), stdGauss);
        y(theta >=range(1) & theta <= -cutoffPoint) = max(y(indexGauss));
    case 1
        indexGauss = theta >= range(2); 
        y(indexGauss) = normpdf(theta(indexGauss), range(2), stdGauss);
        y(theta >= cutoffPoint & theta <= range(2)) = max(y(indexGauss));        
    case 'NA'
end
    
% Normalize to make a probability distribution
scalingFactor = trapz(theta, y);
y = y / scalingFactor;
end