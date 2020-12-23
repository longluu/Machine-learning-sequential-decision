%% Display the fit parameters (9 params)
subjIDAll = {'ll', 'pw', 'eh', 'bh', 'ln', 'at', 'dh'};
for kk = 1 : length(subjIDAll)
    subjID = subjIDAll{kk};
    fileName = ['FitResult-' subjID '.txt'];
    fileID = fopen(fileName);
    myFile = textscan(fileID,'%s','delimiter','\n');
    myFile = myFile{1};
    saveNextLine = 0;
    counter = 1;
    paramsAll = NaN(30, 10);
    for ii = 1 : length(myFile)
        tempLine = myFile{ii};
        if saveNextLine
            paramsAll(counter, :) = str2num(tempLine);
            counter = counter + 1;
            saveNextLine = 0;
        else
            indMatch = strfind(tempLine,'Iteration');
            if ~isempty(indMatch)
                saveNextLine = 1;
            end
        end
    end

    [~, indexSort] = sort(paramsAll(:, 1));
    paramsSort = paramsAll(indexSort, :);
    paramsSort = [paramsSort; mean(paramsSort(1:20, :), 1)];
    currentFolder = pwd;
    fileNameRoot = ['FitResult-' subjID];
    fileID = fopen([fileNameRoot '-extracted' '.txt'],'w');
    for ii = 1 : size(paramsSort, 1)-1
        fprintf(fileID, '%9.2f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f %9.4f %9.4f\r\n', ...
            paramsSort(ii, 1), paramsSort(ii, 2), paramsSort(ii, 3), paramsSort(ii, 4),...
            paramsSort(ii, 5), paramsSort(ii, 6), paramsSort(ii, 7), paramsSort(ii, 8), paramsSort(ii, 9), paramsSort(ii, 10));     
    end
    ii = ii + 1;
    fprintf(fileID, '\n');
    fprintf(fileID, '%9.2f %9.4f %9.4f %16.4f  %10.4f %10.4f %8.4f %9.4f %9.4f %9.4f\r\n', ...
        paramsSort(ii, 1), paramsSort(ii, 2), paramsSort(ii, 3), paramsSort(ii, 4),...
        paramsSort(ii, 5), paramsSort(ii, 6), paramsSort(ii, 7), paramsSort(ii, 8), paramsSort(ii, 9), paramsSort(ii, 10));     

    fclose(fileID);
end
