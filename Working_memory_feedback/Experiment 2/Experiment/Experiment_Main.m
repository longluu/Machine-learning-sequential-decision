function Experiment_Main(params, fileLoad)
% Set some parameters we'll use.
load('thetaStimulus')
if ~isempty(fileLoad)
    dataFile = fullfile('Data', params.subject, 'MainExperiment', [params.experimentName num2str(params.session)], fileLoad);
    load(dataFile);
end
timeWait = params.timeWait;
screenParams = mglDescribeDisplays;   
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
gapAngle = params.gapAngle(1);
sweepAngle = params.sweepAngle;
priorWedgeRGB = params.priorWedgeRGB;
ringRadius1 = params.apertureSize/2 - params.barLength;
ringRadius = [ringRadius1; ringRadius1 - 1];
numBarTotal = params.numBars;
numBarOuter = 16;
numBar = [numBarOuter numBarTotal-numBarOuter];
if params.instruction
    params.timeOutDecision = 25;
end

% Create the data array.  
%  Column 1: true angle difference (population)
%  Column 2: bar presentation time
%  Column 3: SOA
%  Column 4: bar noise level
%  Column 5: bar reference angle
%  Column 6 to 8: data (categorical decision and estimation)
%  Column 9: true angle difference (sample) 
%  Column 10: decision time
%  Column 11: estimation time
nTrials = nAngleDiff*nBarDuration*nSOA*nBarStd*nAngleReference;
if isempty(fileLoad)
    dataTotal = NaN(nTrials, 11);
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
nCategoryGivenTrial = round(params.proportionCategoryGiven * nTrials);
nEstimationTrial = round(params.proportionEstimation * (nTrials - nCategoryGivenTrial));
params.timeWait = timeWait;

% Prevent any key to Matlab terminal
if ~params.instruction
    ListenChar(2)
else
    params.timeOutDecision = 25;
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
    lineMat = rectTextureCreate(lineEstimateSize(1), lineEstimateSize(2), params.barRGB, params.backgroundRGB, [], [], gapLine,[]);
    lineTexture = mglCreateTexture(lineMat);
    boundaryMat = rectTextureCreate(1, boundarySize, [0 0 0], params.backgroundRGB, [], lineEstimateSize(2), gapLine,[]);
    boundaryTexture = mglCreateTexture(boundaryMat);

    barMat = rectTextureCreate(barSize(1), barSize(2), params.barRGB, params.backgroundRGB, [], [], gapLine,[]);
    barTexture = mglCreateTexture(barMat);
    angleVertex = linspace(0,2*pi, params.fixation1(4)+1);
    xVertex = params.fixation1(3) * cos(angleVertex) + params.fixation1(1);
    yVertex = params.fixation1(3) * sin(angleVertex) + params.fixation1(2);
    
    % Create the sound
    fs = 5000;
    fSound1 = 587.33; % D5
    fSound2 = 987.77; % B5
    durationSound = params.feedbackDuration;
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
            
            % Dertermine which condition to use
            decisionGiven = (trialOrder(trialIndex) <= nCategoryGivenTrial);
            
            % boundary line  
            angleReference = dataTotal(trialOrder(trialIndex),5);

            % Pick the angle difference
            if ~params.staircase
                angleDiff = dataTotal(trialOrder(trialIndex),1);
            else
                % Update the value using Quest
            end 
            
            % Prior wedge
            startAngle1 = 90-angleReference+gapAngle;
            startAngle2 = 90-angleReference-gapAngle;
            startAngle3 = 270-angleReference+gapAngle;
            startAngle4 = 270-angleReference-gapAngle;
            
            % Fixation point
            mglClearScreen(0.5); 
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)
            mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                startAngle1, sweepAngle(1), priorWedgeRGB, 100, 2);
            mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                startAngle3, sweepAngle(1), priorWedgeRGB, 100, 2);                         
            mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                startAngle2, sweepAngle(2), priorWedgeRGB, 100, 2);
            mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                startAngle4, sweepAngle(2), priorWedgeRGB, 100, 2);
            mglPolygon(xVertex, yVertex, params.fixation1(5:7))             
            startTimeFixation = mglGetSecs;
            mglFlush
            if params.instruction
                pause
            end
            
            
            % Clear the prior wedge after t(FixationToCue) + t(cue)
            while (mglGetSecs - startTimeFixation) < params.fixationToCueTime + params.cueDuration
            end
            mglClearScreen(0.5);
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference) 
            mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
            mglFlush;                

            % Present bars
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
            barAngles = angleBarMean + thetaStimulus(trialOrder(trialIndex), :);
            dataTotal(trialOrder(trialIndex),9) = mean(angleReference - barAngles); 
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)  
            mglBltTexture(repmat(barTexture,1,numBarTotal),[xBarCenter' yBarCenter'],zeros(1,numBarTotal),zeros(1,numBarTotal),barAngles);
            mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
            while (mglGetSecs - startTimeFixation) < params.fixationDuration
            end            
            mglFlush
            startTimebar = mglGetSecs;            
            while (mglGetSecs - startTimebar) < dataTotal(trialOrder(trialIndex),2)
            end
            if params.instruction
                pause
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
                            if ~decisionGiven
                                dataTotal(trialOrder(trialIndex),6) = -1;
                                keepLooping = false;
                            end
                        % Right (CW/Red)
                        case 6
                            if ~decisionGiven
                                dataTotal(trialOrder(trialIndex),6) = 1;
                                keepLooping = false;
                            end
                        case 4
                            if decisionGiven
                                if wedgeColor > 0.5
                                    dataTotal(trialOrder(trialIndex),6) = 0;
                                else
                                    dataTotal(trialOrder(trialIndex),6) = 1;
                                end
                                keepLooping = false;
                            end 
                        case 2
                            if decisionGiven
                                if wedgeColor > 0.5
                                    dataTotal(trialOrder(trialIndex),6) = 1;
                                else
                                    dataTotal(trialOrder(trialIndex),6) = 0;
                                end
                                keepLooping = false;
                            end 

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
            
            % Record decision time
            dataTotal(trialOrder(trialIndex),10) = mglGetSecs - startTimeDecision;
            
            if repeatFlag
                temp = trialOrder(trialIndex);
                trialOrder(trialIndex) = trialOrder(end);
                trialOrder(end) = temp;
                params.trialOrder = trialOrder;
            else                
                % Get the response for estimation task
                if (trialOrder(trialIndex) <= nCategoryGivenTrial) ...
                        || ((trialOrder(trialIndex) >= nCategoryGivenTrial) ...
                                && (trialOrder(trialIndex) <= nCategoryGivenTrial + nEstimationTrial))
                    startTimeEstimation = mglGetSecs;        
                    mglClearScreen(0.5);
                    mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)
                    mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 

                    mglFlush
                    mglGetKeyEvent
                    keepLooping = true;
                    thetaLine = [];
                    if ~decisionGiven
                        % Present feedback (if turned on)
                        if params.enableFeedback
                            if sign(angleDiff) ~= 0
                                responseCorrect = (sign(angleDiff) == dataTotal(trialOrder(trialIndex),6));
                            else
                                if rand > 0.5
                                    responseCorrect = 0;
                                else
                                    responseCorrect = 1;
                                end
                            end
                            while (mglGetSecs - startTimeEstimation) < params.feedbackLatency
                            end                            
                            if responseCorrect
                                sound(yHigh,fs)
                            else
                                sound(yLow,fs)
                            end
                        else                
                            sound(yHigh,fs)
                        end
                        if responseCorrect
                            dataTotal(trialOrder(trialIndex),8) = dataTotal(trialOrder(trialIndex),6);
                        else
                            dataTotal(trialOrder(trialIndex),8) = -dataTotal(trialOrder(trialIndex),6);
                        end
                    elseif dataTotal(trialOrder(trialIndex),6) == 1
                        sound(yHigh,fs)
                    elseif dataTotal(trialOrder(trialIndex),6) == 0
                        sound(yLow,fs)
                    end

                    while keepLooping
                        buttonState = [Gamepad('GetButton', gamepadIndex, 5) Gamepad('GetButton', gamepadIndex, 6)];
                        switch sum(buttonState)
                            % Update line position
                            case 0
                                mglClearScreen(0.5);
                                xCoordinate1 = Gamepad('GetAxis', gamepadIndex, 3);
                                yCoordinate1 = Gamepad('GetAxis', gamepadIndex, 4);
                                xCoordinate2 = Gamepad('GetAxis', gamepadIndex, 1);
                                yCoordinate2 = Gamepad('GetAxis', gamepadIndex, 2);

                                if (sqrt(xCoordinate1^2+yCoordinate1^2) >= 0.99*AxisMax)
                                    thetaLine = atan2(-yCoordinate1, xCoordinate1);
                                    mglBltTexture(boundaryTexture,[0 0],0,0,angleReference) 
                                    mglBltTexture(lineTexture,[0 0],0,0,rad2deg(thetaLine));
                                    mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
                                    mglFlush
                                elseif (sqrt(xCoordinate2^2+yCoordinate2^2) >= 0.99*AxisMax)
                                    thetaLine = atan2(-yCoordinate2, xCoordinate2);
                                    mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)  
                                    mglBltTexture(lineTexture,[0 0],0,0,rad2deg(thetaLine));
                                    mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
                                    mglFlush
                                end
                            % Finish estimation
                            case 1
                                if ~isempty(thetaLine)
                                    thetaEstimate = wrapTo180(rad2deg(thetaLine));
                                    if (thetaEstimate > 180) || (thetaEstimate < 0) 
                                        dataTotal(trialOrder(trialIndex),7) = thetaEstimate-sign(thetaEstimate)*180;
                                    else
                                        dataTotal(trialOrder(trialIndex),7) = thetaEstimate;
                                    end                                
                                    keepLooping = false;  
                                end
                        end
                        key = mglGetKeyEvent;
                        if ~isempty(key) && (key.keyCode == 13)
                            % Quit
                            error('abort');                                
                        end 

                    end
                end
                
                % Record estimation time
                dataTotal(trialOrder(trialIndex),11) = mglGetSecs - startTimeEstimation;
                
                % Increment the trial counter and the sample counter
                if responseCorrect
                    params.score(1) = params.score(1) + 1;
                    params.score(floor((trialIndex-1)/params.nTrialDisplayScore)+2) = ...
                        params.score(floor((trialIndex-1)/params.nTrialDisplayScore)+2) + 1;
                end                                
                trialIndex = trialIndex + 1;
                params.currentTrialIndex = trialIndex;
            end
            mglClearScreen(0.5);
            mglFlush;

            % Wait for subject's signal for next trial
            if repeatFlag
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
            params.currentTrialIndex = trialIndex;
        end
    else
        % Continued session
        trialOrder = params.trialOrder;
        trialIndex = params.currentTrialIndex; 
        while trialIndex <= nTrials 
            repeatFlag = 0;
            
            % Dertermine which condition to use
            decisionGiven = (trialOrder(trialIndex) <= nCategoryGivenTrial);
            
            % boundary line  
            angleReference = dataTotal(trialOrder(trialIndex),5);

            % Pick the angle difference
            if ~params.staircase
                angleDiff = dataTotal(trialOrder(trialIndex),1);
            else
                % Update the value using Quest
            end 
            
            % Prior wedge
            startAngle1 = 90-angleReference+gapAngle;
            startAngle2 = 90-angleReference-gapAngle;
            startAngle3 = 270-angleReference+gapAngle;
            startAngle4 = 270-angleReference-gapAngle;
            
            % Fixation point
            mglClearScreen(0.5); 
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)
            mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                startAngle1, sweepAngle(1), priorWedgeRGB, 100, 2);
            mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                startAngle3, sweepAngle(1), priorWedgeRGB, 100, 2);                         
            mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                startAngle2, sweepAngle(2), priorWedgeRGB, 100, 2);
            mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                startAngle4, sweepAngle(2), priorWedgeRGB, 100, 2);
            mglPolygon(xVertex, yVertex, params.fixation1(5:7))             
            startTimeFixation = mglGetSecs;
            mglFlush
            if params.instruction
                pause
            end
            
            
            % Clear the prior wedge after t(FixationToCue) + t(cue)
            while (mglGetSecs - startTimeFixation) < params.fixationToCueTime + params.cueDuration
            end
            mglClearScreen(0.5);
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference) 
            mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
            mglFlush;                

            % Present bars
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
            barAngles = angleBarMean + thetaStimulus(trialOrder(trialIndex), :);
            dataTotal(trialOrder(trialIndex),9) = mean(angleReference - barAngles); 
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)  
            mglBltTexture(repmat(barTexture,1,numBarTotal),[xBarCenter' yBarCenter'],zeros(1,numBarTotal),zeros(1,numBarTotal),barAngles);
            mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
            while (mglGetSecs - startTimeFixation) < params.fixationDuration
            end            
            mglFlush
            startTimebar = mglGetSecs;            
            while (mglGetSecs - startTimebar) < dataTotal(trialOrder(trialIndex),2)
            end
            if params.instruction
                pause
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
                            if ~decisionGiven
                                dataTotal(trialOrder(trialIndex),6) = -1;
                                keepLooping = false;
                            end
                        % Right (CW/Red)
                        case 6
                            if ~decisionGiven
                                dataTotal(trialOrder(trialIndex),6) = 1;
                                keepLooping = false;
                            end
                        case 4
                            if decisionGiven
                                if wedgeColor > 0.5
                                    dataTotal(trialOrder(trialIndex),6) = 0;
                                else
                                    dataTotal(trialOrder(trialIndex),6) = 1;
                                end
                                keepLooping = false;
                            end 
                        case 2
                            if decisionGiven
                                if wedgeColor > 0.5
                                    dataTotal(trialOrder(trialIndex),6) = 1;
                                else
                                    dataTotal(trialOrder(trialIndex),6) = 0;
                                end
                                keepLooping = false;
                            end 

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

            % Record decision time
            dataTotal(trialOrder(trialIndex),10) = mglGetSecs - startTimeDecision;
            
            if repeatFlag
                temp = trialOrder(trialIndex);
                trialOrder(trialIndex) = trialOrder(end);
                trialOrder(end) = temp;
                params.trialOrder = trialOrder;
            else                
                % Get the response for estimation task
                if (trialOrder(trialIndex) <= nCategoryGivenTrial) ...
                        || ((trialOrder(trialIndex) >= nCategoryGivenTrial) ...
                                && (trialOrder(trialIndex) <= nCategoryGivenTrial + nEstimationTrial))
                    startTimeEstimation = mglGetSecs;        
                    mglClearScreen(0.5);
                    mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)
                    mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 

                    mglFlush
                    mglGetKeyEvent
                    keepLooping = true;
                    thetaLine = [];
                    if ~decisionGiven
                        % Present feedback (if turned on)
                        if params.enableFeedback
                            if sign(angleDiff) ~= 0
                                responseCorrect = (sign(angleDiff) == dataTotal(trialOrder(trialIndex),6));
                            else
                                if rand > 0.5
                                    responseCorrect = 0;
                                else
                                    responseCorrect = 1;
                                end
                            end
                            while (mglGetSecs - startTimeEstimation) < params.feedbackLatency
                            end                                                        
                            if responseCorrect
                                sound(yHigh,fs)
                            else
                                sound(yLow,fs)
                            end
                        else                
                            sound(yHigh,fs)
                        end
                        if responseCorrect
                            dataTotal(trialOrder(trialIndex),8) = dataTotal(trialOrder(trialIndex),6);
                        else
                            dataTotal(trialOrder(trialIndex),8) = -dataTotal(trialOrder(trialIndex),6);
                        end
                    elseif dataTotal(trialOrder(trialIndex),6) == 1
                        sound(yHigh,fs)
                    elseif dataTotal(trialOrder(trialIndex),6) == 0
                        sound(yLow,fs)
                    end

                    while keepLooping
                        buttonState = [Gamepad('GetButton', gamepadIndex, 5) Gamepad('GetButton', gamepadIndex, 6)];
                        switch sum(buttonState)
                            % Update line position
                            case 0
                                mglClearScreen(0.5);
                                xCoordinate1 = Gamepad('GetAxis', gamepadIndex, 3);
                                yCoordinate1 = Gamepad('GetAxis', gamepadIndex, 4);
                                xCoordinate2 = Gamepad('GetAxis', gamepadIndex, 1);
                                yCoordinate2 = Gamepad('GetAxis', gamepadIndex, 2);

                                if (sqrt(xCoordinate1^2+yCoordinate1^2) >= 0.99*AxisMax)
                                    thetaLine = atan2(-yCoordinate1, xCoordinate1);
                                    mglBltTexture(boundaryTexture,[0 0],0,0,angleReference) 
                                    mglBltTexture(lineTexture,[0 0],0,0,rad2deg(thetaLine));
                                    mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
                                    mglFlush
                                elseif (sqrt(xCoordinate2^2+yCoordinate2^2) >= 0.99*AxisMax)
                                    thetaLine = atan2(-yCoordinate2, xCoordinate2);
                                    mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)  
                                    mglBltTexture(lineTexture,[0 0],0,0,rad2deg(thetaLine));
                                    mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
                                    mglFlush
                                end
                            % Finish estimation
                            case 1
                                if ~isempty(thetaLine)
                                    thetaEstimate = wrapTo180(rad2deg(thetaLine));
                                    if (thetaEstimate > 180) || (thetaEstimate < 0) 
                                        dataTotal(trialOrder(trialIndex),7) = thetaEstimate-sign(thetaEstimate)*180;
                                    else
                                        dataTotal(trialOrder(trialIndex),7) = thetaEstimate;
                                    end                                
                                    keepLooping = false;  
                                end
                        end
                        key = mglGetKeyEvent;
                        if ~isempty(key) && (key.keyCode == 13)
                            % Quit
                            error('abort');                                
                        end 

                    end
                end
                
                % Record estimation time
                dataTotal(trialOrder(trialIndex),11) = mglGetSecs - startTimeEstimation;
                
                % Increment the trial counter and the sample counter
                if responseCorrect
                    params.score(1) = params.score(1) + 1;
                    params.score(floor((trialIndex-1)/params.nTrialDisplayScore)+2) = ...
                        params.score(floor((trialIndex-1)/params.nTrialDisplayScore)+2) + 1;
                end                                
                trialIndex = trialIndex + 1;
                params.currentTrialIndex = trialIndex;
            end
            mglClearScreen(0.5);
            mglFlush;

            % Wait for subject's signal for next trial
            if repeatFlag
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
            params.currentTrialIndex = trialIndex;
        end
    end  
    
    % Save the data
    if isempty(fileLoad)
        dataFolder = fullfile('Data', params.subject, 'MainExperiment', [params.experimentName num2str(params.session)]);
        if ~exist(dataFolder,'dir')
            mkdir(dataFolder)
        end
        dataFile = sprintf('%s/%s-%d.mat', dataFolder, params.experimentName, GetNextDataFileNumber(dataFolder, '.mat'));    
    end
    save(dataFile,'dataTotal','params')
    
    % Close the window
    ListenChar(0)
    mglClose
    
    % Return the screen to dark mode
    Datapixx('Open')
    Datapixx('EnableVideoScanningBacklight');
    Datapixx('RegWrRd')
catch e
    % Save the data
    if isempty(fileLoad)
        dataFolder = fullfile('Data', params.subject, 'MainExperiment', [params.experimentName num2str(params.session)]);
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
    
    % Return the screen to dark mode
    Datapixx('Open')
    Datapixx('EnableVideoScanningBacklight');
    Datapixx('RegWrRd')

end
