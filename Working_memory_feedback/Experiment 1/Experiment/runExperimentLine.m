function runExperimentLine(subjectID)
% This program tests the agreement of computer vision algorithm and human
% perception in scene classification task

%% Define the parameters in the experiment
% General parameters
params.fixation1 = [0 0 0.15 15 1 1 1]; % (disk) origin, radius, numVertex, color
params.fixationDuration = 1;
params.fixationToCueTime = 0.5;
params.barDuration = 0.5;
params.SOA = 1;
params.viewDistance = 835;
params.backgroundRGB = [0.5 0.5 0.5];	% RGB of the background.  All values are in the [0,1] range.
params.intertrialInterval = [0.3 0.6];  % time btw consecutive trials (sec)
params.staircase = 0;
params.nAngle = 10; % for staircase only
params.step = 1; % for staircase only
params.enableFeedback = 1;
params.feedbackDuration = 0.5;
params.refDuration = 1;
params.apertureSize = 5;
params.proportionCategoryGiven = 0; % proportion of trial giving subject the answer
params.proportionEstimation = 1; % proportion of trial having estimation task if answer not given
params.proportionFeedback = 0;
params.timeOutDecision = 4;
params.timeWait = 0.5;
params.sigmaModel = [7.7624 4.0246 2.6583];
params.pseModel = [0.3126 -0.5959 0.2921];
params.useModel = 0; % Simulate the decision using model parameters
params.lineToWedgeInterval = 0.5; % for motor noise exp
params.instruction = 0;
params.usePriorWedge = 0;
params.useReferenceLine = 0;
params.nTrialPerCondition = 55;
params.session = 1;
params.control = 0;
params.nTrialDisplayScore = 10;
params.costFunction = 'a/(b+deltaTheta^2)';
params.costFunctionParam = [20 5];

% Stimulus parameter 
params.boundaryLine = 3; % length of boundary reference line (degree)
params.numBars = 24;
params.barRGB = [1 1 1];
params.barLength = 0.6; % the length of bars (degree visual angle)
params.barStdAngle = [3 18]; 
params.barAngleDiff = -21:3:21; % the angle difference (degree)
params.barAngleReference = linspace(0,180,params.nTrialPerCondition);
params.gapAngle = 0.5; % gap between the cue wedge and reference
params.gapLine = 2; % gap btw estimate line (or aperture) and reference
if params.usePriorWedge
    params.sweepAngle = max(params.barAngleDiff);
else
    params.sweepAngle = 0;
end    
params.cueWedgeRGB = [0 0.8 0; 0.8 0 0];
params.priorWedgeRGB =  [0.6 0.6 0.6];
params.wedgeRadius = [params.gapLine/2 params.gapLine/2+0.2];
params.cueDuration = 0.3;
params.jitterRange = 0.3;
params.lineTrueRGB = [0 1 0];
params.score = zeros(1, 1+ceil(params.nTrialPerCondition * length(params.barStdAngle) * ...
                               length(params.barAngleDiff) / params.nTrialDisplayScore));

% Compute the samples of bar angle
angleBarMean = 0;
sampleSize = params.numBars;
nSampleSearch = 50000;
sampleSelect = NaN(length(params.barStdAngle), params.nTrialPerCondition, sampleSize);
meanPopulation = angleBarMean * ones(nSampleSearch, sampleSize);
for ii = 1 : length(params.barStdAngle)
    sampleAll = normrnd(meanPopulation, params.barStdAngle(ii));
    meanSampleAll = mean(sampleAll, 2);
    [~, ind] = sort(abs(meanSampleAll));
    sampleSelect(ii,:,:) = sampleAll(ind(1:params.nTrialPerCondition),:);
end
params.sampleSelect = sampleSelect;
params.counterSample = zeros(1, length(params.barStdAngle));

% Key mappings
params.leftKey = {'a'};             % Keys accepted for left/up/absent response
params.rightKey = {'l'};            % Keys accepted for right/down/present response

% Experiment info
if nargin < 1
    subjectID = 'dummy';
end
params.subject = subjectID;            % Name of the subject.
if params.proportionCategoryGiven == 1
    params.experimentName = 'DecisionGiven'; % for save directory.
elseif params.proportionCategoryGiven == 0
    params.experimentName = 'Original';
else
    params.experimentName = 'Mixed';
end
if params.control == 1
    params.experimentName = 'ControlReplication';
elseif params.control == 2
    params.experimentName = 'ControlImplicitDecision';
end

%% Run the driver programOriginal-1
fileLoad = '';
Experiment_DriverPTB_motor(params, fileLoad)
