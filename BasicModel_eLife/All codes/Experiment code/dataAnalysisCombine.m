%% Combine two experimental data
clear; 

subjectID = 'vs';

% Load the 'base' data
experimentNumber = 4;
experimentType = 'MainExperiment';
experiment = 'ControlReplication';
session = 1;
extractModelParam = 1;
dataAll = [];
fontSize = 20;
if strcmp(experiment, 'PilotData')
    dataFullPath = fullfile('Data', subjectID, experimentType, experiment, ['Session' num2str(session)], ...
                                ['ConditionDecisionMaking-' num2str(experimentNumber)]);
elseif strcmp(experimentType, 'MotorNoise')
    dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                ['MotorNoise-' num2str(experimentNumber)]);
elseif strcmp(experimentType, 'PerceptNoise')
    dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                ['PerceptNoise-' num2str(experimentNumber)]);                                
else
    dataFullPath = fullfile('Data', subjectID, experimentType,[experiment num2str(session)], ...
                                [experiment '-' num2str(experimentNumber)]);
end
load(dataFullPath);
dataBase = dataTotal;
paramsBase = params;

% Load the 'substitute' data
experimentNumber = 3;
experimentType = 'MainExperiment';
experiment = 'ControlReplication';
session = 1;
extractModelParam = 1;
dataAll = [];
fontSize = 20;
if strcmp(experiment, 'PilotData')
    dataFullPath = fullfile('Data', subjectID, experimentType, experiment, ['Session' num2str(session)], ...
                                ['ConditionDecisionMaking-' num2str(experimentNumber)]);
elseif strcmp(experimentType, 'MotorNoise')
    dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                ['MotorNoise-' num2str(experimentNumber)]);
elseif strcmp(experimentType, 'PerceptNoise')
    dataFullPath = fullfile('Data', subjectID, experimentType, ['Session' num2str(session)], ...
                                ['PerceptNoise-' num2str(experimentNumber)]);                                
else
    dataFullPath = fullfile('Data', subjectID, experimentType,[experiment num2str(session)], ...
                                [experiment '-' num2str(experimentNumber)]);
end
load(dataFullPath);
stdLine = unique(dataTotal(:,4));
angleDiff = (unique(dataTotal(:,1)))';

% Replace the data in 'base' by 'substitute'
for ii = 1 : length(angleDiff)
    dataBase(dataBase(:,1) == angleDiff(ii),:) =  dataTotal(dataTotal(:,1) == angleDiff(ii),:);
end

% Save the data
dataTotal = dataBase;
params = paramsBase;
dataFolder = fullfile('Data', params.subject, 'MainExperiment', [params.experimentName num2str(params.session)]);
if ~exist(dataFolder,'dir')
    mkdir(dataFolder)
end
dataFile = sprintf('%s/%s-%d.mat', dataFolder, params.experimentName, GetNextDataFileNumber(dataFolder, '.mat'));    
save(dataFile,'dataTotal','params')
