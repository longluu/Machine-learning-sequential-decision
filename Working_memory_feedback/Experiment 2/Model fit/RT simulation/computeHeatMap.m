%%%%% Compute the heatmap  of reaction time %%%%%
% Input: estimateData (subject's estimate)
%        rtData (subject's reaction time)
%        binCenter
% Output: heatMap (a vector of heat map showing reaction time as a function
% of estimate)

function heatMap = computeHeatMap(estimateData, rtData, binCenter)
    binWidth = diff(binCenter(1:2));
    binEdge = [binCenter-binWidth/2 max(binCenter)+binWidth/2];
    heatMap = NaN(1, length(binCenter));
    for ii = 1 : length(binCenter)
        heatMap(ii) = nanmean(rtData(estimateData >= binEdge(ii) & estimateData < binEdge(ii+1)));
    end
end