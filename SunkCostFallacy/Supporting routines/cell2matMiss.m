function yMatrix = cell2matMiss(yCell)
% Convert vector cell to matrix, fill in NaN for missing value
nRow = 0;
for ii = 1:length(yCell)
    if nRow < size(yCell{ii}, 1)
        nRow = size(yCell{ii}, 1);
    end
end
yMatrix = NaN(nRow, size(yCell, 2));
for ii = 1:length(yCell)
    if size(yCell{ii}, 1) < nRow
        yMatrix(1:size(yCell{ii}, 1), ii) = yCell{ii};
    else
        yMatrix(:, ii) = yCell{ii};
    end
end