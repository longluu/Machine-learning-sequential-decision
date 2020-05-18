function Experiment_DriverPTB_percept(params, fileLoad)
% Set some parameters we'll use.
params.experimentName = 'PerceptNoise';
screenParams = mglDescribeDisplays;   
if ~isempty(fileLoad)
    dataFile = fullfile('Data', params.subject, 'PerceptNoise',  ['Session' num2str(params.session)], fileLoad);
    load(dataFile);
end    
nAngleDiff = size(params.barAngleDiff,2);
nBarDuration = size(params.barDuration,2);
nSOA = size(params.SOA,2);
nAngleReference = size(params.barAngleReference,2);
nBarStd = size(params.barStdAngle,2);
barLength = visangle2stimsize(params.barLength,0,params.viewDistance,...
    screenParams(2).screenSizeMM(1),screenParams(2).screenSizePixel(1));
lineEstimateLength = visangle2stimsize(params.apertureSize, 0,params.viewDistance,...
    screenParams(2).screenSizeMM(1),screenParams(2).screenSizePixel(1));
lineEstimateSize = [1 lineEstimateLength];
barSize = [1 barLength];
boundarySize = visangle2stimsize(params.boundaryLine,0,params.viewDistance,...
    screenParams(2).screenSizeMM(1),screenParams(2).screenSizePixel(1)) ;
gapLine = visangle2stimsize(params.gapLine,0,params.viewDistance,...
        screenParams(2).screenSizeMM(1),screenParams(2).screenSizePixel(1));
ringRadius1 = params.apertureSize/2 - (params.barLength+params.jitterRange)/2;
ringRadius = [ringRadius1; ringRadius1 - 1];
numBarTotal = params.numBars;
numBarOuter = 16;
numBar = [numBarOuter numBarTotal-numBarOuter];
if params.instruction
    params.timeOutDecision = 25;
end

% Create the data array.  
%  Column 1: angle differences
%  Column 2: bar presentation time
%  Column 3: SOA
%  Column 4: bar aspect ratio
%  Column 5: bar reference angle
%  Column 6 to 8: data (categorical decision and estimation)
nTrials = nAngleDiff*nBarDuration*nSOA*nBarStd*nAngleReference;
if isempty(fileLoad)
    dataTotal = NaN(nTrials, 9);
    index = 1;
    for ii = 1:nAngleDiff
        for jj = 1:nBarDuration
            for ll = 1:nSOA
                for kk = 1:nBarStd	
                    for mm = 1:nAngleReference
                        dataTotal(index,1) = params.barAngleDiff(ii);
                        dataTotal(index,2) = params.barDuration(jj);
                        dataTotal(index,3) = params.SOA(ll);
                        if ~params.staircase
                            dataTotal(index,4) = params.barStdAngle(kk);
                        end
                        dataTotal(index,5) = params.barAngleReference(mm);                
                        index = index + 1;
                    end
                end
            end
        end
    end
end
% Create the bookeeping variable for staircase
if params.staircase
    nTrialPerPsychCurve = nAngleDiff*nAngleReference;
    staircaseTrack = cell(nAspectRatio, 2);
    staircaseIndex = ones(nAspectRatio, nTrialPerPsychCurve);
    responseAll = cell(nAspectRatio, 2);
    staircaseIndex(:, nTrialPerPsychCurve/2:end) = 2;
    staircaseIndex = (Shuffle(staircaseIndex'))';
    staircaseCount = ones(1, nAspectRatio);
    staircaseStep = ones(nAspectRatio, 2) * params.step;
end

% Prevent any key to Matlab terminal
if ~params.instruction
    ListenChar(2)
end

% Set the screen to bright mode
Datapixx('Open')
Datapixx('DisableVideoScanningBacklight');
Datapixx('RegWrRd')

% Open mgl window
mglOpen;

try
    % clear both buffers to gray
    mglClearScreen(0.5);mglFlush;
    mglClearScreen(0.5);mglFlush;

    % Clear the  buffer.
    mglGetKeyEvent;

    % Check gamepad availability
    numGamepads = Gamepad('GetNumGamepads');
    if numGamepads == 0
        error('abort')
    else
        gamepadIndex = 1;
    end
    numButtons = Gamepad('GetNumButtons', gamepadIndex);

    % Present the start text and wait for go signal
    mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);
    mglClearScreen(0.5);  
    mglTextDraw('Press any key to start',[0 0]);
    mglFlush
    mglGetKeyEvent;
    keepLooping = true;
    while keepLooping
        buttonStates = NaN(1,numButtons);
        for ii = 1 : numButtons
            buttonStates(ii) = Gamepad('GetButton', gamepadIndex, ii);
        end
        hatPosition = Gamepad('GetHat', gamepadIndex, 1);
        if (sum(buttonStates) > 0) || (hatPosition ~= 360)
            keepLooping = false;
        end
    end
    mglClearScreen(0.5);  
    mglFlush
    
    % Choose the unit of stimulus display
    mglVisualAngleCoordinates(params.viewDistance/10,screenParams(2).screenSizeMM/10);

    % Create the line and fixation texture
    boundaryMat = rectTextureCreate(1, boundarySize, [0 0 0], params.backgroundRGB, [], lineEstimateSize(2), gapLine, []);
    boundaryTexture = mglCreateTexture(boundaryMat);
    barMat = rectTextureCreate(barSize(1), barSize(2), params.barRGB, params.backgroundRGB, [], [], gapLine, []);
    barTexture = mglCreateTexture(barMat);
    angleVertex = linspace(0,2*pi, params.fixation1(4)+1);
    xVertex = params.fixation1(3) * cos(angleVertex) + params.fixation1(1);
    yVertex = params.fixation1(3) * sin(angleVertex) + params.fixation1(2);
    
    % Create the sound
    fs = 5000;
    fSound1 = 587.33; % D5
    fSound2 = 987.77; % B5
    durationSound = 0.1;
    t = linspace(0,durationSound,fs*durationSound);
    yLow = sin(2*pi*fSound1*t);
    yHigh = sin(2*pi*fSound2*t);
        
    % Start the experiment
    if isempty(fileLoad)
        % New session
        trialOrder = randperm(nTrials);  
        params.trialOrder = trialOrder;
        trialIndex = 1;
        while trialIndex <= nTrials 
            repeatFlag = 0;
            
            % Determine counter index for sample
            indexCounter = find(dataTotal(trialOrder(trialIndex),4) == params.barStdAngle);
            
            % boundary line  
            angleReference = dataTotal(trialOrder(trialIndex),5);
            
            % Ellipse aspect ratio
            aspectRatio = dataTotal(trialOrder(trialIndex),4);
            
            % Pick the angle difference
            if ~params.staircase
                angleDiff = dataTotal(trialOrder(trialIndex),1);
            else
                % Update the staircase
                rowIndex = find(params.ellipseAspectRatio == aspectRatio);
                columnIndex = staircaseIndex(rowIndex, staircaseCount(rowIndex)); 
                staircaseUse = staircaseTrack{rowIndex, columnIndex};
                if isempty(staircaseUse)
                    % First trial
                    if columnIndex == 1
                        staircaseUse = -12;
                    else
                        staircaseUse = 12;
                    end
                else
                    % Half staircase step at 25% and 80% 
                    staircaseLength = length(staircaseUse);
                    if staircaseLength == round(nTrialPerPsychCurve/8)
                        staircaseStep(rowIndex, columnIndex) = staircaseStep(rowIndex, columnIndex)/2;
                    elseif staircaseLength == round(nTrialPerPsychCurve/2.5)
                        staircaseStep(rowIndex, columnIndex) = staircaseStep(rowIndex, columnIndex)/2;
                    end   
                    % Update
                    responseTrack = responseAll{rowIndex, columnIndex};
                    if isempty(responseTrack)
                        downFlag = 1;
                    elseif length(responseTrack) == 1
                        downFlag = responseTrack;
                    else
                        if sum(responseTrack(end-1:end)) == 2
                            if columnIndex == 1
                                downFlag = 0;
                            else
                                downFlag = 1;
                            end
                        elseif responseTrack(end) == 0
                            if columnIndex == 1
                                downFlag = 1;
                            else
                                downFlag = 0;
                            end
                        end
                    end
                    if downFlag
                        staircaseUse(end+1) = staircaseUse(end) - staircaseStep(rowIndex, columnIndex);
                    else
                        staircaseUse(end+1) = staircaseUse(end) + staircaseStep(rowIndex, columnIndex);
                    end
                    if (columnIndex == 1) && (staircaseUse(end)>0)
                        staircaseUse(end) = -staircaseStep(rowIndex, columnIndex);
                    elseif (columnIndex == 2) && (staircaseUse(end)<0)
                        staircaseUse(end) = staircaseStep(rowIndex, columnIndex);                        
                    end
                end
                % Save and increase count
                staircaseTrack{rowIndex, columnIndex} = staircaseUse;
                staircaseCount(rowIndex) = staircaseCount(rowIndex) + 1;
                angleDiff = staircaseUse(end);
                dataTotal(trialOrder(trialIndex),1) = angleDiff;
            end 
            
            % Fixation point
            mglClearScreen(0.5); 
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)
            mglPolygon(xVertex, yVertex, params.fixation1(5:7))             
            mglFlush;    
            startTimeFixation = mglGetSecs;
            if params.instruction
                pause
            end
            
            % Present bar arrays
            angleBarMean = angleReference - angleDiff;
            xBarCenter = [];
            yBarCenter = [];
            mglClearScreen(0.5); 
            for ii = 1 : 2
                barAngleLocation = linspace(0,2*pi,numBar(ii)+1)+rand;
                barAngleLocation(end) = [];
                xBarCenter = [xBarCenter ringRadius(ii) * cos(barAngleLocation)];
                yBarCenter = [yBarCenter ringRadius(ii) * sin(barAngleLocation)];
            end
            xBarCenter = xBarCenter + params.jitterRange * (rand(1,numBarTotal)-0.5);
            yBarCenter = yBarCenter + params.jitterRange * (rand(1,numBarTotal)-0.5);
            barAngles = angleBarMean + squeeze(params.sampleSelect(indexCounter, ...
                                                mod(params.counterSample(indexCounter), params.nTrialPerCondition)+1, :));
            dataTotal(trialOrder(trialIndex),9) = mean(angleReference - barAngles);             
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)  
            mglBltTexture(repmat(barTexture,1,numBarTotal),[xBarCenter' yBarCenter'],zeros(1,numBarTotal),zeros(1,numBarTotal),barAngles');
            mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
            while (mglGetSecs - startTimeFixation) < params.fixationDuration
            end            
            mglFlush
            if params.instruction
                pause
            end            
            startTimeEllipse = mglGetSecs;            
            while (mglGetSecs - startTimeEllipse) < dataTotal(trialOrder(trialIndex),2)
            end
            mglClearScreen(0.5);
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)  
            mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
            mglFlush;
            
            % Get the categorical response
            mglGetKeyEvent;
            keepLooping = true;
            startTimeDecision = mglGetSecs;
            while keepLooping
                % Process any gamepad input.
                buttonStates = NaN(1,numButtons);
                for ii = 1 : numButtons
                    buttonStates(ii) = Gamepad('GetButton', gamepadIndex, ii);
                end
                if sum(buttonStates) == 1
                    keyPosition = find(buttonStates,1);
                    switch keyPosition
                        % Left (CCW/Green)
                        case 5
                            dataTotal(trialOrder(trialIndex),6) = -1;
                            keepLooping = false;
                        % Right (CW/Red)
                        case 6
                            dataTotal(trialOrder(trialIndex),6) = 1;
                            keepLooping = false;
                    end
                end                    
                key = mglGetKeyEvent;
                if ~isempty(key) && (key.keyCode == 13)
                    % Abort.
                    error('abort');
                end    
                if (mglGetSecs - startTimeDecision) > params.timeOutDecision 
                    keepLooping = false;
                    repeatFlag = 1;
                end
            end
            
            if repeatFlag
                temp = trialOrder(trialIndex);
                trialOrder(trialIndex) = trialOrder(end);
                trialOrder(end) = temp;
            else
                % Record response correctness
                if sign(angleDiff) ~= 0
                    responseCorrect = (sign(angleDiff) == dataTotal(trialOrder(trialIndex),6));
                else
                    if rand > 0.5
                        responseCorrect = 0;
                    else
                        responseCorrect = 1;
                    end
                end
                if params.staircase
                    responseAll{rowIndex, columnIndex}(end+1) = responseCorrect;
                end
                if params.instruction
                    pause
                end

                % Present feedback (if turned on)
                if params.enableFeedback
                    if responseCorrect
                        sound(yHigh,fs)
                    else
                        sound(yLow,fs)
                    end
                end
                
                % Increment the counters
                if responseCorrect
                    params.score(1) = params.score(1) + 1;
                    params.score(floor((trialIndex-1)/params.nTrialDisplayScore)+2) = ...
                        params.score(floor((trialIndex-1)/params.nTrialDisplayScore)+2) + 1;
                end                
                trialIndex = trialIndex + 1;
                params.currentTrialIndex = trialIndex;
                params.counterSample(indexCounter) = params.counterSample(indexCounter) + 1;
            end
            
            % Turn off the feedback text.
            mglClearScreen(0.5);
            mglFlush;

            % Wait for subject's signal for next trial
            if mod(trialIndex-1, params.nTrialDisplayScore) == 0
                startTimeWait = mglGetSecs;            
                while (mglGetSecs - startTimeWait) < params.timeWait
                end
                mglClearScreen(0.5);
                mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);
                scoreNormalized = round(100*params.score(1)/(trialIndex-1));
                mglTextDraw(['Your current score is ' num2str(scoreNormalized) '/100'],[0 0]);
                mglFlush
                keepLooping = true;
                while keepLooping
                    buttonStates = NaN(1,numButtons);
                    for ii = 1 : numButtons
                        buttonStates(ii) = Gamepad('GetButton', gamepadIndex, ii);
                    end
                    hatPosition = Gamepad('GetHat', gamepadIndex, 1);
                    if (sum(buttonStates) > 0) || (hatPosition ~= 360)
                        keepLooping = false;
                    end
                end                
            elseif repeatFlag
                startTimeWait = mglGetSecs;            
                while (mglGetSecs - startTimeWait) < params.timeWait
                end
                mglClearScreen(0.5);
                mglPolygon(xVertex, yVertex, params.fixation1(5:7))  
                mglFlush
                keepLooping = true;
                while keepLooping
                    buttonStates = NaN(1,numButtons);
                    for ii = 1 : numButtons
                        buttonStates(ii) = Gamepad('GetButton', gamepadIndex, ii);
                    end
                    hatPosition = Gamepad('GetHat', gamepadIndex, 1);
                    if (sum(buttonStates) > 0) || (hatPosition ~= 360)
                        keepLooping = false;
                    end
                end
            else
                timeToNextTrial = rand*(params.intertrialInterval(2)-params.intertrialInterval(1)) + params.intertrialInterval(1);
                timeLastTrial = mglGetSecs;
                while (mglGetSecs - timeLastTrial) < timeToNextTrial
                end
            end
            
            mglClearScreen(0.5);  
            mglFlush            
        end
    else
        % Continued session
        trialOrder = params.trialOrder;
        trialIndex = params.currentTrialIndex;       
        while trialIndex <= nTrials 
            repeatFlag = 0;
            
            % Determine counter index for sample
            indexCounter = find(dataTotal(trialOrder(trialIndex),4) == params.barStdAngle);
            
            % boundary line  
            angleReference = dataTotal(trialOrder(trialIndex),5);
            
            % Ellipse aspect ratio
            aspectRatio = dataTotal(trialOrder(trialIndex),4);
            
            % Pick the angle difference
            if ~params.staircase
                angleDiff = dataTotal(trialOrder(trialIndex),1);
            else
                % Update the staircase
                rowIndex = find(params.ellipseAspectRatio == aspectRatio);
                columnIndex = staircaseIndex(rowIndex, staircaseCount(rowIndex)); 
                staircaseUse = staircaseTrack{rowIndex, columnIndex};
                if isempty(staircaseUse)
                    % First trial
                    if columnIndex == 1
                        staircaseUse = -12;
                    else
                        staircaseUse = 12;
                    end
                else
                    % Half staircase step at 25% and 80% 
                    staircaseLength = length(staircaseUse);
                    if staircaseLength == round(nTrialPerPsychCurve/8)
                        staircaseStep(rowIndex, columnIndex) = staircaseStep(rowIndex, columnIndex)/2;
                    elseif staircaseLength == round(nTrialPerPsychCurve/2.5)
                        staircaseStep(rowIndex, columnIndex) = staircaseStep(rowIndex, columnIndex)/2;
                    end   
                    % Update
                    responseTrack = responseAll{rowIndex, columnIndex};
                    if isempty(responseTrack)
                        downFlag = 1;
                    elseif length(responseTrack) == 1
                        downFlag = responseTrack;
                    else
                        if sum(responseTrack(end-1:end)) == 2
                            if columnIndex == 1
                                downFlag = 0;
                            else
                                downFlag = 1;
                            end
                        elseif responseTrack(end) == 0
                            if columnIndex == 1
                                downFlag = 1;
                            else
                                downFlag = 0;
                            end
                        end
                    end
                    if downFlag
                        staircaseUse(end+1) = staircaseUse(end) - staircaseStep(rowIndex, columnIndex);
                    else
                        staircaseUse(end+1) = staircaseUse(end) + staircaseStep(rowIndex, columnIndex);
                    end
                    if (columnIndex == 1) && (staircaseUse(end)>0)
                        staircaseUse(end) = -staircaseStep(rowIndex, columnIndex);
                    elseif (columnIndex == 2) && (staircaseUse(end)<0)
                        staircaseUse(end) = staircaseStep(rowIndex, columnIndex);                        
                    end
                end
                % Save and increase count
                staircaseTrack{rowIndex, columnIndex} = staircaseUse;
                staircaseCount(rowIndex) = staircaseCount(rowIndex) + 1;
                angleDiff = staircaseUse(end);
                dataTotal(trialOrder(trialIndex),1) = angleDiff;
            end 
            
            % Fixation point
            mglClearScreen(0.5); 
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)
            mglPolygon(xVertex, yVertex, params.fixation1(5:7))             
            mglFlush;    
            startTimeFixation = mglGetSecs;
            if params.instruction
                pause
            end
            
            % Present bar arrays
            angleBarMean = angleReference - angleDiff;
            xBarCenter = [];
            yBarCenter = [];
            mglClearScreen(0.5); 
            for ii = 1 : 2
                barAngleLocation = linspace(0,2*pi,numBar(ii)+1)+rand;
                barAngleLocation(end) = [];
                xBarCenter = [xBarCenter ringRadius(ii) * cos(barAngleLocation)];
                yBarCenter = [yBarCenter ringRadius(ii) * sin(barAngleLocation)];
            end
            xBarCenter = xBarCenter + params.jitterRange * (rand(1,numBarTotal)-0.5);
            yBarCenter = yBarCenter + params.jitterRange * (rand(1,numBarTotal)-0.5);
            barAngles = angleBarMean + squeeze(params.sampleSelect(indexCounter, ...
                                                mod(params.counterSample(indexCounter), params.nTrialPerCondition)+1, :));
            dataTotal(trialOrder(trialIndex),9) = mean(angleReference - barAngles);             
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)  
            mglBltTexture(repmat(barTexture,1,numBarTotal),[xBarCenter' yBarCenter'],zeros(1,numBarTotal),zeros(1,numBarTotal),barAngles');
            mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
            while (mglGetSecs - startTimeFixation) < params.fixationDuration
            end            
            mglFlush
            if params.instruction
                pause
            end            
            startTimeEllipse = mglGetSecs;            
            while (mglGetSecs - startTimeEllipse) < dataTotal(trialOrder(trialIndex),2)
            end
            mglClearScreen(0.5);
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)  
            mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
            mglFlush;
            
            % Get the categorical response
            mglGetKeyEvent;
            keepLooping = true;
            startTimeDecision = mglGetSecs;
            while keepLooping
                % Process any gamepad input.
                buttonStates = NaN(1,numButtons);
                for ii = 1 : numButtons
                    buttonStates(ii) = Gamepad('GetButton', gamepadIndex, ii);
                end
                if sum(buttonStates) == 1
                    keyPosition = find(buttonStates,1);
                    switch keyPosition
                        % Left (CCW/Green)
                        case 5
                            dataTotal(trialOrder(trialIndex),6) = -1;
                            keepLooping = false;
                        % Right (CW/Red)
                        case 6
                            dataTotal(trialOrder(trialIndex),6) = 1;
                            keepLooping = false;
                    end
                end                    
                key = mglGetKeyEvent;
                if ~isempty(key) && (key.keyCode == 13)
                    % Abort.
                    error('abort');
                end    
                if (mglGetSecs - startTimeDecision) > params.timeOutDecision 
                    keepLooping = false;
                    repeatFlag = 1;
                end
            end
            
            if repeatFlag
                temp = trialOrder(trialIndex);
                trialOrder(trialIndex) = trialOrder(end);
                trialOrder(end) = temp;
            else
                % Record response correctness
                if sign(angleDiff) ~= 0
                    responseCorrect = (sign(angleDiff) == dataTotal(trialOrder(trialIndex),6));
                else
                    if rand > 0.5
                        responseCorrect = 0;
                    else
                        responseCorrect = 1;
                    end
                end
                if params.staircase
                    responseAll{rowIndex, columnIndex}(end+1) = responseCorrect;
                end
                if params.instruction
                    pause
                end

                % Present feedback (if turned on)
                if params.enableFeedback
                    if responseCorrect
                        sound(yHigh,fs)
                    else
                        sound(yLow,fs)
                    end
                end
                
                % Increment the counters
                if responseCorrect
                    params.score(1) = params.score(1) + 1;
                    params.score(floor((trialIndex-1)/params.nTrialDisplayScore)+2) = ...
                        params.score(floor((trialIndex-1)/params.nTrialDisplayScore)+2) + 1;
                end                                
                trialIndex = trialIndex + 1;
                params.currentTrialIndex = trialIndex;
                params.counterSample(indexCounter) = params.counterSample(indexCounter) + 1;
            end
            
            % Turn off the feedback text.
            mglClearScreen(0.5);
            mglFlush;

            % Wait for subject's signal for next trial
            if mod(trialIndex-1, params.nTrialDisplayScore) == 0
                startTimeWait = mglGetSecs;            
                while (mglGetSecs - startTimeWait) < params.timeWait
                end
                mglClearScreen(0.5);
                mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);
                scoreNormalized = round(100*params.score(1)/(trialIndex-1));
                mglTextDraw(['Your current score is ' num2str(scoreNormalized) '/100'],[0 0]);
                mglFlush
                keepLooping = true;
                while keepLooping
                    buttonStates = NaN(1,numButtons);
                    for ii = 1 : numButtons
                        buttonStates(ii) = Gamepad('GetButton', gamepadIndex, ii);
                    end
                    hatPosition = Gamepad('GetHat', gamepadIndex, 1);
                    if (sum(buttonStates) > 0) || (hatPosition ~= 360)
                        keepLooping = false;
                    end
                end                
            elseif repeatFlag
                startTimeWait = mglGetSecs;            
                while (mglGetSecs - startTimeWait) < params.timeWait
                end
                mglClearScreen(0.5);
                mglPolygon(xVertex, yVertex, params.fixation1(5:7))  
                mglFlush
                keepLooping = true;
                while keepLooping
                    buttonStates = NaN(1,numButtons);
                    for ii = 1 : numButtons
                        buttonStates(ii) = Gamepad('GetButton', gamepadIndex, ii);
                    end
                    hatPosition = Gamepad('GetHat', gamepadIndex, 1);
                    if (sum(buttonStates) > 0) || (hatPosition ~= 360)
                        keepLooping = false;
                    end
                end
            else
                timeToNextTrial = rand*(params.intertrialInterval(2)-params.intertrialInterval(1)) + params.intertrialInterval(1);
                timeLastTrial = mglGetSecs;
                while (mglGetSecs - timeLastTrial) < timeToNextTrial
                end
            end
            
            mglClearScreen(0.5);  
            mglFlush            
        end  
    end
    % Save the data
    if params.staircase 
        params.staircaseTrack = staircaseTrack;
    end
    if isempty(fileLoad)
        dataFolder = fullfile('Data', params.subject, 'PerceptNoise', ['Session' num2str(params.session)]);
        if ~exist(dataFolder,'dir')
            mkdir(dataFolder)
        end
        dataFile = sprintf('%s/%s-%d.mat', dataFolder, params.experimentName, GetNextDataFileNumber(dataFolder, '.mat'));    
    end
    save(dataFile,'dataTotal','params')
    
    % Close the window
    ListenChar(0)
    mglClose
catch e
    % Save the data
    if params.staircase 
        params.staircaseTrack = staircaseTrack;
    end
    if isempty(fileLoad)
        dataFolder = fullfile('Data', params.subject, 'PerceptNoise', ['Session' num2str(params.session)]);
        if ~exist(dataFolder,'dir')
            mkdir(dataFolder)
        end
        dataFile = sprintf('%s/%s-%d.mat', dataFolder, params.experimentName, GetNextDataFileNumber(dataFolder, '.mat'));    
    end
    save(dataFile,'dataTotal','params')

    % Close the window
    ListenChar(0)
    mglClose
    
    % Report error
    if strcmp(e.message, 'abort')
        fprintf('- Experiment aborted, some data saved.\n');
    else
        rethrow(e);
    end    
end
