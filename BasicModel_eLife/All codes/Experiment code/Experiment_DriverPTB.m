function Experiment_DriverPTB(params)
% Set some parameters we'll use.
screenParams = mglDescribeDisplays;   
nBlocks = ceil(params.nTrialPerCondition/length(params.ellipseAngleReference));
nAngleDiff = size(params.ellipseAngleDiff,2);
nEllipseDuration = size(params.ellipseDuration,2);
nSOA = size(params.SOA,2);
nAngleReference = size(params.ellipseAngleReference,2);
if ~params.staircase
    % Don't use staircase
    nAspectRatio = size(params.ellipseAspectRatio,2);
else
    % Use staircase
end
radiusAxisLong = visangle2stimsize(params.ellipseLongAxis,0,params.viewDistance,...
    screenParams(2).screenSizeMM(1),screenParams(2).screenSizePixel(1));
fixationRadius = visangle2stimsize(params.fixation1(3),0,params.viewDistance,...
    screenParams(2).screenSizeMM(1),screenParams(2).screenSizePixel(1));
lineEstimateSize = [2 2*radiusAxisLong];

% Create the data array.  
%  Column 1: angle differences
%  Column 2: ellipse presentation time
%  Column 3: SOA
%  Column 4: ellipse aspect ratio
%  Column 5: ellipse reference angle
%  Column 6 to 3*nBlocks+5: data (categorical decision and estimation)
nTrials = nAngleDiff*nEllipseDuration*nSOA*nAspectRatio*nAngleReference;
dataTotal = NaN(nTrials, 3*nBlocks+5);
index = 1;
for ii = 1:nAngleDiff
    for jj = 1:nEllipseDuration
        for ll = 1:nSOA
            for kk = 1:nAspectRatio	
                for mm = 1:nAngleReference
                    dataTotal(index,1) = params.ellipseAngleDiff(ii);
                    dataTotal(index,2) = params.ellipseDuration(jj);
                    dataTotal(index,3) = params.SOA(ll);
                    if ~params.staircase
                        dataTotal(index,4) = params.ellipseAspectRatio(kk);
                    end
                    dataTotal(index,5) = params.ellipseAngleReference(mm);                
                    index = index + 1;
                end
            end
        end
    end
end
nCategoryGivenTrial = round(params.proportionCategoryGiven * nTrials);
nEstimationTrial = round(params.proportionEstimation * (nTrials - nCategoryGivenTrial));

% Prevent any key to Matlab terminal
ListenChar(2)

% Open mgl window
mglOpen;

try
    % clear both buffers to gray
    mglClearScreen(0.5);mglFlush;
    mglClearScreen(0.5);mglFlush;

    % Clear the keyboard buffer.
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
    fixation1 = ellipseTextureCreate(fixationRadius,fixationRadius,0,params.fixation1(4:end),params.backgroundRGB, [4 4 1.5]);
    fixationTexture1 = mglCreateTexture(fixation1);
    fixation2 = ellipseTextureCreate(fixationRadius,fixationRadius,0,params.fixation2(4:end),params.backgroundRGB, [4 4 1.5]);
    fixationTexture2 = mglCreateTexture(fixation2);
    lineMat = rectTextureCreate(lineEstimateSize(1), lineEstimateSize(2), params.ellipseRGB, params.backgroundRGB, fixation1, [],[]);
    lineTexture = mglCreateTexture(lineMat);

    for blockIndex = 1 : nBlocks
        trialOrder = randperm(nTrials);         
        for trialIndex = 1 : nTrials
            % Dertermine which condition to use
            decisionGiven = (trialOrder(trialIndex) <= nCategoryGivenTrial);
            
            % Fixation point
            mglClearScreen(0.5); 
            if decisionGiven
                mglBltTexture(fixationTexture2,[0 0],0,0,0);
            else
                mglBltTexture(fixationTexture1,[0 0],0,0,0);
            end
            mglFlush;    
            startTimeFixation = mglGetSecs;
            while (mglGetSecs - startTimeFixation) < params.fixationDuration
            end
            
            % Pick the angle difference
            if ~params.staircase
                angleDiff = dataTotal(trialOrder(trialIndex),1);
            else
                % Update the value using Quest
            end 
            
            
            % First ellipse 
            referenceDuration = params.refDuration;
            if ~params.referenceLine
                radiusAxisShort1 = radiusAxisLong / params.referenceAspectRatio;
                angleReference = dataTotal(trialOrder(trialIndex),5);
                ellipseMat1 = ellipseTextureCreate(radiusAxisLong,radiusAxisShort1,0,params.ellipseRGB,params.backgroundRGB);
                ellipseTexture1 = mglCreateTexture(ellipseMat1);
                mglClearScreen(0.5);  
                mglBltTexture(ellipseTexture1,[0 0],0,0,angleReference);  
            else
                angleReference = dataTotal(trialOrder(trialIndex),5);
                mglClearScreen(0.5);  
                mglBltTexture(lineTexture,[0 0],0,0,angleReference);                  
            end
            mglFlush
            startTimeEllipse1 = mglGetSecs;
            if trialOrder(trialIndex) <= nCategoryGivenTrial
                % Give the subject categorical decision
                startAngle = 90-angleReference+sign(angleDiff)*10;
                sweepAngle = 90;
                if angleDiff > 0
                    directionCue = 'clockwise';
                    mglGluPartialDisk( 0, 0, params.ellipseLongAxis+0.2, params.ellipseLongAxis+0.4,...
                        startAngle, sweepAngle, [0 1 0], 100, 2);
                elseif angleDiff < 0
                    directionCue = 'counterclockwise';                    
                    mglGluPartialDisk( 0, 0, params.ellipseLongAxis+0.2, params.ellipseLongAxis+0.4,...
                        startAngle, -sweepAngle, [0 1 0], 100, 2);
                else
                    if rand > 0.5
                        startAngle = 90-angleReference+10;
                        directionCue = 'clockwise';                        
                        mglGluPartialDisk( 0, 0, params.ellipseLongAxis+0.2, params.ellipseLongAxis+0.4,...
                            startAngle, sweepAngle, [0 1 0], 100, 2);
                    else
                        startAngle = 90-angleReference-10;
                        directionCue = 'counterclockwise';                        
                        mglGluPartialDisk( 0, 0, params.ellipseLongAxis+0.2, params.ellipseLongAxis+0.4,...
                            startAngle, -sweepAngle, [0 1 0], 100, 2);
                    end
                end
                if strcmp(directionCue,'clockwise')
                    dataTotal(trialOrder(trialIndex),blockIndex+2*nBlocks+5) = 1;
                elseif strcmp(directionCue,'counterclockwise')
                    dataTotal(trialOrder(trialIndex),blockIndex+2*nBlocks+5) = -1;
                end
                if ~params.referenceLine
                    mglBltTexture(ellipseTexture1,[0 0],0,0,angleReference);  
                else
                    mglBltTexture(lineTexture,[0 0],0,0,angleReference);                  
                end
                while (mglGetSecs - startTimeEllipse1) < params.cueDuration
                end
                mglFlush
            end             
            while (mglGetSecs - startTimeEllipse1) < referenceDuration
            end
            mglClearScreen(0.5);
            if decisionGiven
                mglBltTexture(fixationTexture2,[0 0],0,0,0);
            else
                mglBltTexture(fixationTexture1,[0 0],0,0,0);
            end
            mglFlush;
            if ~params.referenceLine
                mglDeleteTexture(ellipseTexture1)
            end
            
            % Second ellipse
            radiusAxisShort2 = radiusAxisLong / dataTotal(trialOrder(trialIndex),4);
            angleEllipse2 = angleReference - angleDiff;
            ellipseMat2 = ellipseTextureCreate(radiusAxisLong,radiusAxisShort2,0,params.ellipseRGB,params.backgroundRGB);
            ellipseTexture2 = mglCreateTexture(ellipseMat2);
            mglClearScreen(0.5);              
            mglBltTexture(ellipseTexture2,[0 0],0,0,angleEllipse2);
            while (mglGetSecs - startTimeEllipse1) < dataTotal(trialOrder(trialIndex),3)
            end
            mglFlush
            startTimeellipse2 = mglGetSecs;
            while (mglGetSecs - startTimeellipse2) < dataTotal(trialOrder(trialIndex),2)
            end
            mglClearScreen(0.5);
            if decisionGiven
                mglBltTexture(fixationTexture2,[0 0],0,0,0);
            else
                mglBltTexture(fixationTexture1,[0 0],0,0,0);
            end   
            mglFlush;
            mglDeleteTexture(ellipseTexture2)  
            
            if trialOrder(trialIndex) > nCategoryGivenTrial
                % Get the categorical response
                mglGetKeyEvent;
                keepLooping = true;
                while keepLooping
                    % Process any gamepad input.
                    buttonStates = NaN(1,numButtons);
                    for ii = 1 : numButtons
                        buttonStates(ii) = Gamepad('GetButton', gamepadIndex, ii);
                    end
                    if sum(buttonStates) == 1
                        keyPosition = find(buttonStates,1);
                        switch keyPosition
                            % Left
                            case 5
                                dataTotal(trialOrder(trialIndex),blockIndex+5) = -1;
                                keepLooping = false;

                            % Right
                            case 6
                                dataTotal(trialOrder(trialIndex),blockIndex+5) = 1;
                                keepLooping = false;
                        end
                    end
                    key = mglGetKeyEvent;
                    if ~isempty(key) && (key.keyCode == 13)
                        % Abort.
                        error('abort');
                    end                    
                end

                % Present the feedback if no estimation needed
                if (trialOrder(trialIndex) > nCategoryGivenTrial + nEstimationTrial)                   
                    responseCorrect = (sign(angleDiff) == dataTotal(trialOrder(trialIndex),blockIndex+5));
                    if responseCorrect
                        textName = 'correct';
                        textRGB = [0 1 0];
                    else
                        textName = 'incorrect';
                        textRGB = [1 0 0];
                    end

                    % Enable the appropriate feedback text.
                    mglTextSet('Helvetica',32,textRGB,0,0,0,0,0,0,0);
                    mglClearScreen(0.5);  
                    mglTextDraw(textName,[0 0]);
                    mglFlush
                    timeStart = mglGetSecs;
                    while (mglGetSecs - timeStart) < params.feedbackDuration
                    end

                    % Turn off the feedback text.
                    mglClearScreen(0.5);
                    mglFlush;
                end  
            end                
                
            % Get the response for estimation task
            if (trialOrder(trialIndex) <= nCategoryGivenTrial) ...
                    || ((trialOrder(trialIndex) >= nCategoryGivenTrial) ...
                            && (trialOrder(trialIndex) <= nCategoryGivenTrial + nEstimationTrial))
                mglClearScreen(0.5);
                mglBltTexture(fixationTexture2,[0 0],0,0,0);
                mglFlush
                mglGetKeyEvent
                keepLooping = true;
                thetaLine = [];
                while keepLooping
                    buttonState = Gamepad('GetButton', gamepadIndex, 5);
                    switch buttonState
                        % Update line position
                        case 0
                            mglClearScreen(0.5);
                            xCoordinate = Gamepad('GetAxis', gamepadIndex, 3);
                            yCoordinate = Gamepad('GetAxis', gamepadIndex, 4);
                            if (sqrt(xCoordinate^2+yCoordinate^2) >= 0.95*AxisMax)
                                thetaLine = atan2(-yCoordinate, xCoordinate);
                                mglBltTexture(lineTexture,[0 0],0,0,rad2deg(thetaLine));
                                mglFlush
                            end
                        % Finish estimation
                        case 1
                            if ~isempty(thetaLine)
                                thetaEstimate = wrapTo180(rad2deg(thetaLine));
                                if (thetaEstimate > 180) || (thetaEstimate < 0) 
                                    dataTotal(trialOrder(trialIndex),blockIndex+nBlocks+5) = thetaEstimate-sign(thetaEstimate)*180;
                                else
                                    dataTotal(trialOrder(trialIndex),blockIndex+nBlocks+5) = thetaEstimate;
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
            mglClearScreen(0.5);
            mglFlush;
            
            % Wait some time before next trial
            timeToNextTrial = rand*(params.intertrialInterval(2)-params.intertrialInterval(1)) + params.intertrialInterval(1);
            timeLastTrial = mglGetSecs;
            while (mglGetSecs - timeLastTrial) < timeToNextTrial
            end
        end
    end
    
    % Save the data
    dataFolder = fullfile('Data', params.subject, 'MainExperiment', ['Session' num2str(params.session)]);
    if ~exist(dataFolder,'dir')
        mkdir(dataFolder)
    end
    dataFile = sprintf('%s/%s-%d.mat', dataFolder, params.experimentName, GetNextDataFileNumber(dataFolder, '.mat'));    
    save(dataFile,'dataTotal','params')
    
    % Close the window
    ListenChar(0)
    mglClose
catch e
    % Save the data
    dataFolder = fullfile('Data', params.subject, 'MainExperiment', ['Session' num2str(params.session)]);
    if ~exist(dataFolder,'dir')
        mkdir(dataFolder)
    end
    dataFile = sprintf('%s/%s-%d.mat', dataFolder, params.experimentName, GetNextDataFileNumber(dataFolder, '.mat'));    
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
