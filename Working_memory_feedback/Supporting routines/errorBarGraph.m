function hBar = errorBarGraph(dataMat, lowerBound, upperBound, colorIndex)
lowerBound = lowerBound(:);
upperBound = upperBound(:);

% Creating the bar graph
hBar = bar(dataMat,'BarWidth',0.8, 'EdgeColor', 'none');
hold on;

% Set the color
for ii = 1 : size(dataMat, 2)
%     hBar(ii).FaceColor = colorIndex(ii,:);
    set(hBar(ii), 'FaceColor', colorIndex(ii,:));
end

% Finding the number of groups and the number of bars in each group
ngroups = size(dataMat, 1);
nbars = size(dataMat, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
counter = 1;
for ii = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars);
    
    % Plot the errorbar
    for jj = 1 : length(x)
        plot([x(jj) x(jj)], [lowerBound(counter) upperBound(counter)], 'k', 'LineWidth', 2)
        counter = counter + 1;
    end
end