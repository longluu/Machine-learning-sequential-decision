function Experiment_DriverPTB_trainArray_Line(params, fileLoad)
% Set some parameters we'll use.
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
ringRadius1 = params.apertureSize/2 - params.barLength;
ringRadius = [ringRadius1; ringRadius1 - 1];
numBarTotal = params.numBars;
numBarOuter = 16;
numBar = [numBarOuter numBarTotal-numBarOuter];
priorWedgeRGB = params.priorWedgeRGB;
sweepAngle = params.sweepAngle;
gapAngle = params.gapAngle(1);

% Create the data array.  
%  Column 1: angle differences
%  Column 2: bar presentation time
%  Column 3: SOA
%  Column 4: bar aspect ratio
%  Column 5: bar reference angle
%  Column 6 to 3+5: data (categorical decision and estimation)
nTrials = nAngleDiff*nBarDuration*nSOA*nBarStd*nAngleReference;
if isempty(fileLoad)
    dataTotal = NaN(nTrials, 8);
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
else
    dataFile = fullfile('Data', params.subject, 'MainExperiment', [params.experimentName num2str(params.session)], fileLoad);
    load(dataFile);
end
nCategoryGivenTrial = round(params.proportionCategoryGiven * nTrials);
nEstimationTrial = round(params.proportionEstimation * (nTrials - nCategoryGivenTrial));
params.timeWait = timeWait;

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
    durationSound = 0.3;
    t = linspace(0,durationSound,fs*durationSound);
    yLow = sin(2*pi*fSound1*t);
    yHigh = sin(2*pi*fSound2*t);
    
    % Start the experiment
    if isempty(fileLoad)
        % New session
        trialOrder = randperm(nTrials);  
        params.trialOrder = trialOrder;
        trialIndex = 0;
        while trialIndex <= nTrials 
            % Dertermine which condition to use
            decisionGiven = (trialOrder(trialIndex) <= nCategoryGivenTrial);

            % boundary line  
            angleReference = dataTotal(trialOrder(trialIndex),5);
            
            % Prior wedge
            startAngle1 = 90-angleReference+gapAngle;
            startAngle2 = 90-angleReference-gapAngle;
            startAngle3 = 270-angleReference+gapAngle;
            startAngle4 = 270-angleReference-gapAngle;

            % Pick the angle difference
            if ~params.staircase
                angleDiff = dataTotal(trialOrder(trialIndex),1);
            else
                % Update the value using Quest
            end 
                        
            % Fixation point
            mglClearScreen(0.5); 
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)
            mglPolygon(xVertex, yVertex, params.fixation1(5:7))             
            startTimeFixation = mglGetSecs;
            mglFlush
            if params.instruction
                pause
            end
            
            % Present bars
            angleBarStd = dataTotal(trialOrder(trialIndex),4);
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
            barAngles = angleBarMean + angleBarStd * randn(1, numBarTotal);
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
            startTimeEndStim = mglGetSecs;

            % Get the response for estimation task
            if (trialOrder(trialIndex) <= nCategoryGivenTrial) ...
                    || ((trialOrder(trialIndex) >= nCategoryGivenTrial) ...
                            && (trialOrder(trialIndex) <= nCategoryGivenTrial + nEstimationTrial))
                mglClearScreen(0.5);
                mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)
                mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
                if params.usePriorWedge == 1;
                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                        startAngle1, sweepAngle, priorWedgeRGB, 100, 2);
                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                        startAngle3, sweepAngle, priorWedgeRGB, 100, 2);                         
                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                        startAngle2, -sweepAngle, priorWedgeRGB, 100, 2);
                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                        startAngle4, -sweepAngle, priorWedgeRGB, 100, 2);
                end
                while (mglGetSecs - startTimeEndStim) < params.lineToWedgeInterval
                end                
                mglFlush
                mglGetKeyEvent
                keepLooping = true;
                thetaLine = [];
                if ~decisionGiven
                    sound(yLow,fs)
                elseif dataTotal(trialOrder(trialIndex),6) == 1
                    sound(yLow,fs)
                elseif dataTotal(trialOrder(trialIndex),6) == 0
                    sound(yHigh,fs)
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
                                if params.usePriorWedge == 1;
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle1, sweepAngle, priorWedgeRGB, 100, 2);
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle3, sweepAngle, priorWedgeRGB, 100, 2);                         
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle2, -sweepAngle, priorWedgeRGB, 100, 2);
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle4, -sweepAngle, priorWedgeRGB, 100, 2);
                                end                                
                                mglFlush
                            elseif (sqrt(xCoordinate2^2+yCoordinate2^2) >= 0.99*AxisMax)
                                thetaLine = atan2(-yCoordinate2, xCoordinate2);
                                mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)  
                                mglBltTexture(lineTexture,[0 0],0,0,rad2deg(thetaLine));
                                mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
                                if params.usePriorWedge == 1;
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle1, sweepAngle, priorWedgeRGB, 100, 2);
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle3, sweepAngle, priorWedgeRGB, 100, 2);                         
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle2, -sweepAngle, priorWedgeRGB, 100, 2);
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle4, -sweepAngle, priorWedgeRGB, 100, 2);
                                end                                                                
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

            % Increment the trial counter
            trialIndex = trialIndex + 1;
            params.currentTrialIndex = trialIndex;
            mglClearScreen(0.5);
            mglFlush;

            % Wait for subject's signal for next trial
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
            mglClearScreen(0.5);  
            mglFlush
            params.currentTrialIndex = trialIndex;
        end
    else
        % Continued session
        trialOrder = params.trialOrder;
        trialIndex = params.currentTrialIndex; 
        while trialIndex <= nTrials 
            % Dertermine which condition to use
            decisionGiven = (trialOrder(trialIndex) <= nCategoryGivenTrial);

            % boundary line  
            angleReference = dataTotal(trialOrder(trialIndex),5);
            
            % Prior wedge
            startAngle1 = 90-angleReference+gapAngle;
            startAngle2 = 90-angleReference-gapAngle;
            startAngle3 = 270-angleReference+gapAngle;
            startAngle4 = 270-angleReference-gapAngle;

            % Pick the angle difference
            if ~params.staircase
                angleDiff = dataTotal(trialOrder(trialIndex),1);
            else
                % Update the value using Quest
            end 
                        
            % Fixation point
            mglClearScreen(0.5); 
            mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)
            mglPolygon(xVertex, yVertex, params.fixation1(5:7))             
            startTimeFixation = mglGetSecs;
            mglFlush
            if params.instruction
                pause
            end
            
            % Present bars
            angleBarStd = dataTotal(trialOrder(trialIndex),4);
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
            barAngles = angleBarMean + angleBarStd * randn(1, numBarTotal);
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
            startTimeEndStim = mglGetSecs;

            % Get the response for estimation task
            if (trialOrder(trialIndex) <= nCategoryGivenTrial) ...
                    || ((trialOrder(trialIndex) >= nCategoryGivenTrial) ...
                            && (trialOrder(trialIndex) <= nCategoryGivenTrial + nEstimationTrial))
                mglClearScreen(0.5);
                mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)
                mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
                if params.usePriorWedge == 1;
                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                        startAngle1, sweepAngle, priorWedgeRGB, 100, 2);
                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                        startAngle3, sweepAngle, priorWedgeRGB, 100, 2);                         
                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                        startAngle2, -sweepAngle, priorWedgeRGB, 100, 2);
                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                        startAngle4, -sweepAngle, priorWedgeRGB, 100, 2);
                end
                while (mglGetSecs - startTimeEndStim) < params.lineToWedgeInterval
                end                
                mglFlush
                mglGetKeyEvent
                keepLooping = true;
                thetaLine = [];
                if ~decisionGiven
                    sound(yLow,fs)
                elseif dataTotal(trialOrder(trialIndex),6) == 1
                    sound(yLow,fs)
                elseif dataTotal(trialOrder(trialIndex),6) == 0
                    sound(yHigh,fs)
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
                                if params.usePriorWedge == 1;
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle1, sweepAngle, priorWedgeRGB, 100, 2);
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle3, sweepAngle, priorWedgeRGB, 100, 2);                         
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle2, -sweepAngle, priorWedgeRGB, 100, 2);
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle4, -sweepAngle, priorWedgeRGB, 100, 2);
                                end                                
                                mglFlush
                            elseif (sqrt(xCoordinate2^2+yCoordinate2^2) >= 0.99*AxisMax)
                                thetaLine = atan2(-yCoordinate2, xCoordinate2);
                                mglBltTexture(boundaryTexture,[0 0],0,0,angleReference)  
                                mglBltTexture(lineTexture,[0 0],0,0,rad2deg(thetaLine));
                                mglPolygon(xVertex, yVertex, params.fixation1(5:7)) 
                                if params.usePriorWedge == 1;
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle1, sweepAngle, priorWedgeRGB, 100, 2);
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle3, sweepAngle, priorWedgeRGB, 100, 2);                         
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle2, -sweepAngle, priorWedgeRGB, 100, 2);
                                    mglGluPartialDisk( 0, 0, params.apertureSize/2+params.wedgeRadius(1), params.apertureSize/2+params.wedgeRadius(2),...
                                        startAngle4, -sweepAngle, priorWedgeRGB, 100, 2);
                                end                                                                
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

            % Increment the trial counter
            trialIndex = trialIndex + 1;
            params.currentTrialIndex = trialIndex;
            mglClearScreen(0.5);
            mglFlush;
            
            % Wait for subject's signal for next trial
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
            mglClearScreen(0.5);  
            mglFlush
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
