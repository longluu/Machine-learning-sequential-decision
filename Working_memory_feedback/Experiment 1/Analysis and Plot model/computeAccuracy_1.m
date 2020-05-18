%%%%%%%%%%%%%%%%%%%%%%% Compute the error of the model and subjects' estimates %%%%%%%%%%%%%%%%%%%%%%

%% Compute the subject error
subjectID = 'average';
experimentNumber = 1:5;
experimentType = 'MainExperiment';
experiment = 'Original';
session = 1;
dataAll = [];
fontSize = 20;
lineWidth = 2;

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
biasMotor = zeros(1,15);  
if strcmp(experimentType, 'MotorNoise')
    a = dataAll(:,7);
    dataAll(isnan(a),:) = [];    
end

% Extract data
%  Column 1: angle differences (theta1 - theta2)
%  Column 2: bar presentation time
%  Column 3: SOA
%  Column 4: bar Stimulus noise
%  Column 5: bar reference angle
%  Column 6: CW/CCW OR Red/Green
%  Column 7: estimated angle
%  Column 8: given CW/CCW
stdNoiseLevel = unique(dataAll(:,4));
angleDiff = (unique(dataAll(:,1)))';
dataZeroDiff = dataAll(dataAll(:,1) == 0,:);
angleEstimate = dataAll(:, 7);
if strcmp(experimentType, 'MotorNoise')
    stdNoiseLevel = 0;
end

% Average across conditions
angleTrue = dataAll(:,5) - dataAll(:,1);
angleTrue(angleTrue>180 | angleTrue<0) = angleTrue(angleTrue>180 | angleTrue<0)...
                                            - 180 * sign(angleTrue(angleTrue>180 | angleTrue<0));
indexCorrect = find(abs(angleTrue-angleEstimate) > 70);
angleEstimate(indexCorrect) = angleEstimate(indexCorrect) + sign(angleTrue(indexCorrect)-angleEstimate(indexCorrect))*180;    
angleDiffEst = dataAll(:, 5) - angleEstimate;
indexAdjust1 = abs(angleDiffEst)>90;
indexAdjust2 = sign(angleDiffEst) ~= sign(dataAll(:, 1));
indexAdjust = indexAdjust1 & indexAdjust2;
angleDiffEst(indexAdjust) = angleDiffEst(indexAdjust) - 180 * sign(angleDiffEst(indexAdjust));
indexCorrect = dataAll(:, 6) == dataAll(:, 8);
errorCorrect_Data = sum((angleDiffEst(indexCorrect) - dataAll(indexCorrect,1)).^2) / sum(~isnan(dataAll(indexCorrect,1)));
errorIncorrect_Data = sum((angleDiffEst(~indexCorrect) - dataAll(~indexCorrect,1)).^2) / sum(~isnan(dataAll(~indexCorrect,1)));

%% Compute the mean errors of models
flagSC = 1; % 1: self-conditioned model
           % 0: standard Bayes
includeIncongruentTrials = 0;
incorrectType = 2; % 1: flip the decision bit
                   % 2: flip the estimates
                   % 3: resample the measurement mm until getting a consistent sample

dstep = 0.1;
paramsAll = [5.3271    9.8390           0.0000     45.8895     0.0085   0.9968    3.3234];
lapseRate = paramsAll(3);

% stimulus orientation
thetaStim = -21:0.1:21; % 
thetaStim = round(thetaStim, -log10(dstep));

% sensory noise
stdSensory = paramsAll(1:2);

% memory recall noise
stdMemory = paramsAll(5);
stdMemoryIncorrect = sqrt(stdMemory^2 + 0^2);

% motor noise;
stdMotor = paramsAll(7);

% priors
smoothFactor = paramsAll(6);

% LOOP - noise levels
pC = [0.5, 0.5]'; % [cw ccw]
pthcw = paramsAll(4);
pthccw = -paramsAll(4);

rangeth = [-60 60];
th = rangeth(1):dstep:rangeth(2);
th = round(th, -log10(dstep));
nth = length(th);

pthGC = zeros(2,nth);

if flagSC
    pthGC(1,:) = TukeyWindow([0 pthcw], 0, smoothFactor, th);
    pthGC(2,:) = TukeyWindow([pthccw 0], 1, smoothFactor, th);
else
    pth = (TukeyWindow([0 pthcw], 0, smoothFactor, th) + TukeyWindow([pthccw 0], 1, smoothFactor, th))/2;
    pth(th==0) = 0;
    pth(th==0) = max(pth);
    pthGC(1,:) = pth;
    pthGC(2,:) = pth;
end

errorIncorrect_Model = NaN(2, 3);
errorCorrect_Model = NaN(2, 1);
for kk=1:length(stdSensory)  
    rangeM = [min(thetaStim)-5*stdSensory(kk) max(thetaStim)+5*stdSensory(kk)];
    if rangeM(2) < rangeth(2)
        rangeM = rangeth;
    end
    nm = 1000;
    m = linspace(rangeM(1), rangeM(2), nm);

    nmm = 1200;
    rangeMM = [min(rangeM)-6*stdMemory max(rangeM)+6*stdMemory];
    if rangeMM(2) < rangeth(2)
        rangeMM = rangeth;
    end        
    mm = linspace(rangeMM(1), rangeMM(2), nmm);

    nmr = nm;
    mr = m;

    M = repmat(m',1,nth);
    MM_m = repmat(mm',1,nm);
    MM_th = repmat(mm',1,nth); 
    MM_ths = repmat(mm',1,length(thetaStim));
    MR_mm = repmat(mr', 1, nmm);
    MR_th = repmat(mr', 1, nth);
    THm = repmat(th, nm, 1); 
    THmm = repmat(th, nmm, 1);
    THmr = repmat(th, nmr, 1);
    THSmm = repmat(thetaStim, nmm, 1);
    
    %% Correct trials
    % Generative (forward)
    % orientation noise
    pmGth = exp(-((M-THm).^2)./(2*stdSensory(kk)^2));
    pmGth = pmGth./(repmat(sum(pmGth,1),nm,1)); 

    % Inference
    % 1: categorical judgment
    PCGm = (pthGC * pmGth') .* repmat(pC,1,nm);
    % fix the issue when sensory noise is too low
    indFirstNonZero = find(PCGm(2,:), 1);
    PCGm(2, 1: indFirstNonZero-1) = PCGm(2, indFirstNonZero);
    indLastNonZero = find(PCGm(1,:), 1, 'last');
    PCGm(1, indLastNonZero+1:end) = PCGm(1, indLastNonZero);
    PCGm = PCGm./(repmat(sum(PCGm,1),2,1));
    % max posterior decision
    PChGm = round(PCGm);
    % marginalization
    PChGtheta = PChGm * pmGth(:, ismember(th, thetaStim));
    PChGtheta_lapse = lapseRate + (1 - 2*lapseRate) * PChGtheta;
    PChGtheta_lapse = PChGtheta_lapse ./ repmat(sum(PChGtheta_lapse, 1), 2, 1);
    
    % 2: estimation
    pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); % p(mm|th) = N(th, sm^2 + smm^2)
    pmmGth = pmmGth./(repmat(sum(pmmGth,1),nmm,1)); 
    pthGmmChcw = (pmmGth.*repmat(pthGC(1,:),nmm,1))';
    pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
    pthGmmChcw(isnan(pthGmmChcw)) = 0;

    pthGmmChccw = (pmmGth.*repmat(pthGC(2,:),nmm,1))';
    pthGmmChccw = pthGmmChccw./repmat(sum(pthGmmChccw,1),nth,1);
    pthGmmChccw(isnan(pthGmmChccw)) = 0;

    EthChcw = th * pthGmmChcw;
    EthChccw = th * pthGmmChccw;
    % discard repeating/decreasing values (required for interpolation) 
    indKeepCw = 1:length(EthChcw);
    while sum(diff(EthChcw)<=0) >0
        indDiscardCw = [false diff(EthChcw)<=0];
        EthChcw(indDiscardCw) = [];
        indKeepCw(indDiscardCw) = [];
    end
    indKeepCcw = 1:length(EthChccw);
    while sum(diff(EthChccw)<=0) >0
        indDiscardCcw = [diff(EthChccw)<=0 false];
        EthChccw(indDiscardCcw) = [];
        indKeepCcw(indDiscardCcw) = [];
    end
    
    a = 1./gradient(EthChcw,dstep);
    % memory noise
    pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemory^2)); 
    pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1));   

    % attention marginalization: compute distribution only over those ms that lead to cw decision!
    pmmGthChcw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim)));
    b = repmat(a',1,length(thetaStim)) .* pmmGthChcw(indKeepCw, :);        
    pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
    % add motor noise
    pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
    pthhGthChcw(pthhGthChcw < 0) = 0; 
    
    a = 1./gradient(EthChccw,dstep);
    % attention marginalization: compute distribution only over those ms that lead to cw decision!
    pmmGthChccw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim)));        
    b = repmat(a',1,length(thetaStim)) .* pmmGthChccw(indKeepCcw, :);        
    pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
    % add motor noise
    pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
    pthhGthChccw(pthhGthChccw < 0) = 0; 
    
    pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
    pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            
    
    if includeIncongruentTrials == 0
        % modify psychometric curve p(Chat|theta, Congruent) ~ p(Congruent| Chat, theta) * p(Chat|Theta)
        pCongruentGcwTh = sum(pthhGthChcw(th' >= 0, :));
        pCongruentGccwTh = sum(pthhGthChccw(th' <= 0, :));
        PChGtheta_lapse_new = PChGtheta_lapse .* [pCongruentGcwTh; pCongruentGccwTh];
        PChGtheta_lapse_new = PChGtheta_lapse_new ./ repmat(sum(PChGtheta_lapse_new, 1), 2, 1);

        % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
        pthhGthChccw(th'>= 0, :) = 0;
        pthhGthChcw(th'< 0, :) = 0;
    else
        PChGtheta_lapse_new = PChGtheta_lapse;
    end
    
    if incorrectType == 2
        pthhGthChcw_Incorrect = pthhGthChcw;
        pthhGthChccw_Incorrect = pthhGthChccw;
        
        % remove correct trials
        pthhGthChcw_Incorrect(:, thetaStim > 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          hcw(:, thetaStim < 0) = 0;
        pthhGthChccw_Incorrect(:, thetaStim < 0) = 0;
        
        % flip the estimate
        pthhGthChcw_Incorrect = flipud(pthhGthChcw_Incorrect);
        pthhGthChccw_Incorrect = flipud(pthhGthChccw_Incorrect);
    end
    
    % remove incorrect trials
    pthhGthChcw(:, thetaStim < 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          hcw(:, thetaStim < 0) = 0;
    pthhGthChccw(:, thetaStim > 0) = 0;
    
    
    pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
    pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            
    mthhGthChcw_correct = th * pthhGthChcw_norm;
    mthhGthChccw_correct = th * pthhGthChccw_norm;
    mthhGthChcw_correct(thetaStim < 0) = NaN;
    mthhGthChccw_correct(thetaStim > 0) = NaN;
    
    pthhANDth_correct = pthhGthChcw.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw.*repmat(PChGtheta_lapse(2,:),nth,1);
    pthhANDth_correct(:, thetaStim == 0) = pthhANDth_correct(:, thetaStim == 0) /2;
    pthhANDth_correct = pthhANDth_correct / sum(pthhANDth_correct(:));
    
    errorCorrect = ((repmat(th', 1, length(thetaStim)) - repmat(thetaStim, length(th), 1)) .^ 2) .* pthhANDth_correct;
    errorCorrect_Model(kk) = sum(errorCorrect(:));
    
    %% Incorrect type 1
        pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemoryIncorrect^2))); % p(mm|th) = N(th, sm^2 + smm^2)
        pmmGth = pmmGth./(repmat(sum(pmmGth,1),nmm,1)); 
        
        pthGmmChcw = (pmmGth.*repmat(pthGC(2,:),nmm,1))';
        pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
        pthGmmChcw(isnan(pthGmmChcw)) = 0;

        pthGmmChccw = (pmmGth.*repmat(pthGC(1,:),nmm,1))';
        pthGmmChccw = pthGmmChccw./repmat(sum(pthGmmChccw,1),nth,1);
        pthGmmChccw(isnan(pthGmmChccw)) = 0;

        EthChcw = th * pthGmmChcw;
        EthChccw = th * pthGmmChccw;
        % discard repeating/decreasing values (required for interpolation) 
        indKeepCw = 1:length(EthChcw);
        while sum(diff(EthChcw)<=0) >0
            indDiscardCw = [false diff(EthChcw)<=0];
            EthChcw(indDiscardCw) = [];
            indKeepCw(indDiscardCw) = [];
        end
        indKeepCcw = 1:length(EthChccw);
        while sum(diff(EthChccw)<=0) >0
            indDiscardCcw = [diff(EthChccw)<=0 false];
            EthChccw(indDiscardCcw) = [];
            indKeepCcw(indDiscardCcw) = [];
        end

        a = 1./gradient(EthChcw,dstep);
        % memory noise
        pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemoryIncorrect^2)); 
        pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1));   

        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChcw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim)));
        b = repmat(a',1,length(thetaStim)) .* pmmGthChcw(indKeepCw, :);        

        pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChcw(pthhGthChcw < 0) = 0; 

        a = 1./gradient(EthChccw,dstep);
        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChccw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim)));        
        b = repmat(a',1,length(thetaStim)) .* pmmGthChccw(indKeepCcw, :);        
        pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChccw(pthhGthChccw < 0) = 0; 
        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            

        if includeIncongruentTrials == 0
            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'<= 0, :) = 0;
            pthhGthChcw(th'> 0, :) = 0;
        end

        % remove 'correct' trials
        pthhGthChcw(:, thetaStim > 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          hcw(:, thetaStim < 0) = 0;
        pthhGthChccw(:, thetaStim < 0) = 0;
        pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
        pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
        pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
        pthhANDth_incorrect = pthhGthChcw_norm.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw_norm.*repmat(PChGtheta_lapse(2,:),nth,1);
        pthhANDth_incorrect(:, thetaStim == 0) = pthhANDth_incorrect(:, thetaStim == 0)/2;
        pthhANDth_incorrect = pthhANDth_incorrect / sum(pthhANDth_incorrect(:));

        errorIncorrect = ((repmat(th', 1, length(thetaStim)) - repmat(thetaStim, length(th), 1)) .^ 2) .* pthhANDth_incorrect;
        errorIncorrect_Model(kk, 1) = sum(errorIncorrect(:));
    
    %% Incorrect type 2
        pthhGthChcw = pthhGthChcw_Incorrect;
        pthhGthChccw = pthhGthChccw_Incorrect;
        pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
        pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
        pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
        pthhANDth_incorrect = pthhGthChcw_norm.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw_norm.*repmat(PChGtheta_lapse(2,:),nth,1);
        pthhANDth_incorrect(:, thetaStim == 0) = pthhANDth_incorrect(:, thetaStim == 0)/2;
        pthhANDth_incorrect = pthhANDth_incorrect / sum(pthhANDth_incorrect(:));

        errorIncorrect = ((repmat(th', 1, length(thetaStim)) - repmat(thetaStim, length(th), 1)) .^ 2) .* pthhANDth_incorrect;
        errorIncorrect_Model(kk, 2) = sum(errorIncorrect(:));
        
    %% Incorrect type 3
        % Measurement m is degraded by memory noise p(mm|th, Chat) = p(mm|th) = N(th, sm^2 + smm^2)           
        pmrGth = exp(-((MR_th-THmr).^2)./(2*(stdSensory(kk)^2 + stdMemory^2)));
        pmrGth = pmrGth./(repmat(sum(pmrGth,1),nmr,1)); 

        pthGmrChcw = (pmrGth.*repmat(pthGC(2,:),nmr,1))';
        pthGmrChcw = pthGmrChcw./repmat(sum(pthGmrChcw,1),nth,1);
        pthGmrChcw(isnan(pthGmrChcw)) = 0;

        pthGmrChccw = (pmrGth.*repmat(pthGC(1,:),nmr,1))';
        pthGmrChccw = pthGmrChccw./repmat(sum(pthGmrChccw,1),nth,1);
        pthGmrChccw(isnan(pthGmrChccw)) = 0;

        EthChcw = th * pthGmrChcw;
        EthChccw = th * pthGmrChccw;
        % discard repeating/decreasing values (required for interpolation) 
        indKeepCw = 1:length(EthChcw);
        while sum(diff(EthChcw)<=0) >0
            indDiscardCw = [false diff(EthChcw)<=0];
            EthChcw(indDiscardCw) = [];
            indKeepCw(indDiscardCw) = [];
        end
        indKeepCcw = 1:length(EthChccw);
        while sum(diff(EthChccw)<=0) >0
            indDiscardCcw = [diff(EthChccw)<=0 false];
            EthChccw(indDiscardCcw) = [];
            indKeepCcw(indDiscardCcw) = [];
        end

        % Resample mm until we have a sample that is consistent with feedback
        % p(mr|mm, theta, Chat)
        pmrGmmth = exp(-((MR_mm-repmat(mm, nmr, 1)).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); 

        pmrGmmthChcw = pmrGmmth;
        pmrGmmthChcw(mr > 0, :) = 0;
        pmrGmmthChcw = pmrGmmthChcw./(repmat(sum(pmrGmmthChcw,1),nmr,1));

        pmrGmmthChccw = pmrGmmth;
        pmrGmmthChccw(mr < 0, :) = 0;
        pmrGmmthChccw = pmrGmmthChccw./(repmat(sum(pmrGmmthChccw,1),nmr,1));

        % Marginalize over mm that lead to cw decision to compute likelihood p(mr|theta, Chat)
        pmrGthChcw = pmrGmmthChcw * pmmGthChcw;   
        pmrGthChcw = pmrGthChcw ./ (repmat(sum(pmrGthChcw,1),nmr,1)); 
        pmrGthChcw(isnan(pmrGthChcw)) = 0;

        pmrGthChccw = pmrGmmthChccw * pmmGthChccw;
        pmrGthChccw = pmrGthChccw ./ (repmat(sum(pmrGthChccw,1),nmr,1)); 
        pmrGthChccw(isnan(pmrGthChccw)) = 0;

        a = 1./gradient(EthChcw,dstep);
        b = repmat(a',1,length(thetaStim)) .* pmrGthChcw(indKeepCw, :);        

        pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChcw(pthhGthChcw < 0) = 0; 

        a = 1./gradient(EthChccw,dstep);
        b = repmat(a',1,length(thetaStim)) .* pmrGthChccw(indKeepCcw, :);        
        pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChccw(pthhGthChccw < 0) = 0; 

        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            

        if includeIncongruentTrials == 0
            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'<= 0, :) = 0;
            pthhGthChcw(th'> 0, :) = 0;
        end

        % remove 'correct' trials
        pthhGthChcw(:, thetaStim > 0) = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          hcw(:, thetaStim < 0) = 0;
        pthhGthChccw(:, thetaStim < 0) = 0; 
        
        pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
        pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
        pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
        pthhANDth_incorrect = pthhGthChcw_norm.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw_norm.*repmat(PChGtheta_lapse(2,:),nth,1);
        pthhANDth_incorrect(:, thetaStim == 0) = pthhANDth_incorrect(:, thetaStim == 0)/2;
        pthhANDth_incorrect = pthhANDth_incorrect / sum(pthhANDth_incorrect(:));

        errorIncorrect = ((repmat(th', 1, length(thetaStim)) - repmat(thetaStim, length(th), 1)) .^ 2) .* pthhANDth_incorrect;
        errorIncorrect_Model(kk, 3) = sum(errorIncorrect(:));        
end

%% Plot the error
maxY = max(errorIncorrect_Data)+5;
figure
subplot(1, 2, 1)
bar(1:2, [errorCorrect_Data mean(errorCorrect_Model)]);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Data', 'Model'}, 'FontSize', 20)
ylabel('Mean squared error (deg^2)')
xlim([0 3])
ylim([0 maxY]);
title('Correct trials')

subplot(1, 2, 2)
bar(1:4, [errorIncorrect_Data mean(errorIncorrect_Model, 1)]);
set(gca, 'XTick', 1:4, 'XTickLabel', {'Data', 'Flip decision', 'Flip estimate', 'Resample'}, 'FontSize', 20)
ylim([0 maxY]);
ylabel('Mean squared error (deg^2)')
title('Incorrect trials')

