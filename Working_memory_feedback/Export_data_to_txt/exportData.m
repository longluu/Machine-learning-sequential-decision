%%%%%%%%%%%%%% Export data to text file %%%%%%%%%%%%%%
subjectIDAll = {'ll', 'an', 'ep', 'jp', 'kc', 'll', 'pw', 'eh', 'bh', 'ln', 'at', 'dh'};
subjectNumber = [1 2:5 1 6:11];
experimentNumerID = [ones(1, 5) 2*ones(1, 7)];

currentFolder = pwd;
fileNameRoot = 'DataAllExperiment';
fileNumber = GetNextDataFileNumber(currentFolder, '.txt');
fileID = fopen([fileNameRoot '-' num2str(fileNumber) '.txt'],'w');
dir_root = "C:\Users\longluu\Documents\GitHub\Machine-learning-sequential-decision\Working_memory_feedback";
subjMissTrialOrder = zeros(size(subjectIDAll));
for mm = 1 : length(subjectIDAll)
    subjectID = subjectIDAll{mm};
    experimentNumber = 1;
    experimentType = 'MainExperiment';
    experiment = 'Original';
    if experimentNumerID(mm) == 1
        data_path = fullfile(dir_root, 'Experiment 1');
    elseif experimentNumerID(mm) == 2
        data_path = fullfile(dir_root, 'Experiment 2');
    end
    session = 1;

    dataFullPath = fullfile(data_path, 'Data', subjectID, experimentType, [experiment num2str(session)], ...
                                [experiment '-' num2str(experimentNumber)]);          
    load(dataFullPath);   
    
    %% Export data to txt file
    %  Column 1: true angle difference (population)
    %  Column 2: bar presentation time
    %  Column 3: SOA
    %  Column 4: bar noise level
    %  Column 5: bar reference angle
    %  Column 6: subject's categorical decision
    %  Column 7: subject's estimate
    %  Column 8: true category
    %  Column 9: true angle difference (sample) 
    %  Column 10: decision time
    %  Column 11: estimation time
    fprintf(fileID, 'Experiment: %1s,    Subject: %1s \r\n', num2str(experimentNumerID(mm)), num2str(subjectNumber(mm)));     
    if experimentNumerID(mm) == 1
        fprintf(fileID, 'AngleDiff   Noise  RefAngle   SubjDecision  SubjEstimate  Feedback   TrialOrder \n');
        for ii = 1 : size(dataTotal, 1)
            fprintf(fileID, '%6.0f %8.0f %10.2f %10.0f %16.5f %10.0f %10.0f \r\n', dataTotal(ii, 1), ...
                        dataTotal(ii, 4), dataTotal(ii, 5), dataTotal(ii, 6), dataTotal(ii, 7), dataTotal(ii, 8), params.trialOrder(ii));
        end     
    else
        fprintf(fileID, 'AngleDiff   Noise  RefAngle   SubjDecision  SubjEstimate  Feedback   DecisionTime   EstimateTime   TrialOrder \n');
        for ii = 1 : size(dataTotal, 1)
            fprintf(fileID, '%6.0f %8.0f %10.2f %10.0f %16.5f %10.0f %14.5f %14.5f %10.0f \r\n', dataTotal(ii, 1), ...
                        dataTotal(ii, 4), dataTotal(ii, 5), dataTotal(ii, 6), dataTotal(ii, 7), dataTotal(ii, 8), dataTotal(ii, 10), dataTotal(ii, 11), params.trialOrder(ii));
        end    
    end
end
fclose(fileID);



