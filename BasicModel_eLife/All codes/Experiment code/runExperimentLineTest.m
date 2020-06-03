function runExperimentLineTest(subjectID)


%% Define the parameters in the experiment
% General parameters
params.fixation1 = [0 0 0.15 15 1 1 1]; % (disk) origin, radius, numVertex, color
params.fixationDuration = 1.3;
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
params.nTrialPerCondition = 40;

% Stimulus parameter 
params.boundaryLine = 3; % length of boundary reference line (degree)
params.numBars = 24;
params.barRGB = [1 1 1];
params.barLength = 0.6; % the length of bars (degree visual angle)
params.barStdAngle = [0 6 18]; 
params.barAngleDiff = linspace(-21,21,15); % the angle difference (degree)
params.barAngleReference = linspace(0,180,params.nTrialPerCondition);
params.gapAngle = 0.5; % gap between the cue wedge and reference
params.gapLine = 2; % gap btw estimate line (or aperture) and reference
params.sweepAngle = max(params.barAngleDiff);
params.cueWedgeRGB = [0 0.8 0; 0.8 0 0];
params.priorWedgeRGB =  [0.6 0.6 0.6];
params.wedgeRadius = [params.gapLine /2 params.gapLine/2+0.2];
params.cueDuration = 0.3;
params.jitterRange = 0.3;
params.lineTrueRGB = [0 1 0];

% Key mappings
params.leftKey = {'a'};             % Keys accepted for left/up/absent response
params.rightKey = {'l'};            % Keys accepted for right/down/present response

% Experiment info
if nargin < 1
    subjectID = 'aj';
end
params.subject = subjectID;            % Name of the subject.
if params.proportionCategoryGiven == 1
    params.experimentName = 'DecisionGiven'; % for save directory.
elseif params.proportionCategoryGiven == 0
    params.experimentName = 'Original';
else
    params.experimentName = 'Mixed';
end
params.session = 1;
params.control = 0;
if params.control == 1
    params.experimentName = 'ControlReplication';
elseif params.control == 2
    params.experimentName = 'ControlImplicitDecision';
end

%% Run the driver program
fileLoad = 'Original-1';
Experiment_DriverPTB_Main_Line(params, fileLoad)
