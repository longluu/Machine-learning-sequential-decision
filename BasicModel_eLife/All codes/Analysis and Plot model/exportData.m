%%%%%%%%%%%%%% Export data to text file %%%%%%%%%%%%%%
subjectIDAll = {'ll', 'sy', 'cz', 'vs', 'as', 'll', 'xfl', 'aj', 'zw', 'skb', 'll', 'xfl', 'aj', 'zw', 'skb'};
subjectNumber = [1 2:5 1 6:9 1 6:9];
experimentNumerID = [ones(1, 5) 2*ones(1, 5) 3*ones(1, 5)];

% subjectIDAll = {'ll'};
% subjectNumber = [1];
% experimentNumerID = [1];
currentFolder = pwd;
fileNameRoot = 'DataAllExperiment';
fileNumber = GetNextDataFileNumber(currentFolder, '.txt');
fileID = fopen([fileNameRoot '-' num2str(fileNumber) '.txt'],'w');
subjMissTrialOrder = zeros(size(subjectIDAll));
for mm = 1 : length(subjectIDAll)
    subjectID = subjectIDAll{mm};
    experimentNumber = 1;
    experimentType = 'MainExperiment';
    if experimentNumerID(mm) == 1
        experiment = 'ControlReplication';
    elseif experimentNumerID(mm) == 2
        experiment = 'Original';
    else
        experiment = 'DecisionGiven';
    end
    session = 1;
    dataAll = [];
    trialOrder = [];
    indexTrialOrder = [];
    maxTrialOrder = 0;
    includeIncorrectTrial = 1;
    earlyTrial = 0;

    for ii = 1 : length(experimentNumber)
        if strcmp(experiment, 'PilotData')
            dataFullPath = fullfile('Data', subjectID, experimentType, experiment, ['Session' num2str(session)], ...
                                        ['ConditionDecisionMaking-' num2str(experimentNumber(ii))]);
        elseif strcmp(experimentType, 'MotorNoise')
            dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                        ['MotorNoise-' num2str(experimentNumber(ii))]);
        elseif strcmp(experimentType, 'PerceptNoise')
            dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                        ['PerceptNoise-' num2str(experimentNumber(ii))]);                                
        elseif strcmp(experimentType, 'TrainArray')
            dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                        ['TrainArray-' num2str(experimentNumber(ii))]);  
        else
            dataFullPath = fullfile('Data', subjectID, experimentType, [experiment num2str(session)], ...
                                        [experiment '-' num2str(experimentNumber(ii))]);          
        end
        load(dataFullPath);   
        trialOrder = [trialOrder params.trialOrder+maxTrialOrder];
        maxTrialOrder = max(trialOrder);
        indexTrialOrder = [indexTrialOrder 1:length(params.trialOrder)];
        tempTrialOrder = params.trialOrder;
    %     if earlyTrial == 1
    %         dataTotal(tempTrialOrder(length(dataTotal) - round(percentInclude * length(dataTotal)) : end), 7) = NaN;
    %     else
    %         dataTotal(tempTrialOrder(1 : round(percentInclude * length(dataTotal))), 7) = NaN;
    %     end
        if ii == 1
            dataAll = dataTotal;
            if strcmp(experimentType, 'PerceptNoise') && isfield(params, 'staircaseTrack')
                staircaseTrack = params.staircaseTrack;
            end
        else
            dataAll = [dataAll; dataTotal];
            if strcmp(experimentType, 'PerceptNoise') && isfield(params, 'staircaseTrack')
                staircaseTrack = [staircaseTrack; params.staircaseTrack];
            end
        end
    end  
    if length(params.trialOrder) == size(dataAll, 1)
        excludeIndex = isnan(dataAll(:, 6)) | isnan(dataAll(:, 7));
        dataAll(excludeIndex, :) = [];
        params.trialOrder(excludeIndex) = [];
    else
        excludeIndex = isnan(dataAll(:, 6)) | isnan(dataAll(:, 7));
        dataAll(excludeIndex, :) = [];
        subjMissTrialOrder(mm) = 1;
        params.trialOrder = zeros(size(dataAll, 1), 1);        
    end
    
    %% Export data to txt file
    %  Column 1: angle differences (theta1 - theta2)
    %  Column 4: bar Stimulus noise
    %  Column 5: bar reference angle
    %  Column 6: CW/CCW OR Red/Green
    %  Column 7: estimated angle
    %  Column 8: given CW/CCW
    fprintf(fileID, 'Experiment: %1s,    Subject: %1s \r\n', num2str(experimentNumerID(mm)), num2str(subjectNumber(mm)));     
    if experimentNumerID(mm) == 3
        for ii = 1 : size(dataAll, 1)
            fprintf(fileID, '%6.0f %6.0f %10.2f %6.0f %12.5f %6.0f %10.0f \r\n', dataAll(ii, 1), ...
                        dataAll(ii, 4), dataAll(ii, 5), dataAll(ii, 6), dataAll(ii, 7),  params.trialOrder(ii), dataAll(ii, 8));
        end     
    else
        for ii = 1 : size(dataAll, 1)
            fprintf(fileID, '%6.0f %6.0f %10.2f %6.0f %12.5f %10.0f \r\n', dataAll(ii, 1), ...
                        dataAll(ii, 4), dataAll(ii, 5), dataAll(ii, 6), dataAll(ii, 7), params.trialOrder(ii));
        end    
    end
end
fclose(fileID);
subjectIDAll(logical(subjMissTrialOrder))



