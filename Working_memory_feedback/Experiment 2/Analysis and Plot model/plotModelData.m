%%%% Plot the model along with subjects' data
%% Set the condition
subjectID = {'average'};
experimentNumber = 1:7;
experimentType = 'MainExperiment';
experiment = 'Original';
session = 1;
dataAll = [];
fontSize = 20;
lineWidth = 3;
if strcmp(subjectID{1}, 'average')
    subjectIDD = '';
else
    subjectIDD = subjectID{1};
end

includeIncongruentTrials = 0; 

correctType = 1; % 1: no resampling
                 % 2: resampling (center m, variance: memory)
incorrectType = 4; % 1: flip the decision bit
                   % 2: flip the estimates
                   % 3: resample (center: memory, variance: sensory + memory)
                   % 4: resample (center: sensory, variance: sensory + memory)
flagSC = 1; % 1: Self-conditioned 
            % 0: Full Bayes
            
if correctType == 1
    paramsFit = [7.3559    9.9743           0.0000     25.7467   -15.1199   0.1192    4.9975    0.9760    0.5276]; % no resample    
elseif correctType == 2
    paramsFit = [6.9011    8.7172           0.0000     33.0973   -19.3675   0.8412    2.6349    0.9895    0.5503]; % resample    
end
stdSensory = paramsFit(1:2);
priorRange = paramsFit(4);
smoothFactor = paramsFit(8);
stdMemory = paramsFit(6);
if stdMemory == 0
    stdMemory = 0.0001;
end
lapseRate = paramsFit(3);
stdMotor = paramsFit(7);

%% Load data
for ii = 1 : length(experimentNumber)
    if strcmp(experiment, 'PilotData')
        dataFullPath = fullfile('Data', subjectID{1}, experimentType, experiment, ['Session' num2str(session)], ...
                                    ['ConditionDecisionMaking-' num2str(experimentNumber(ii))]);
    elseif strcmp(experimentType, 'MotorNoise')
        dataFullPath = fullfile('Data', subjectID{1}, experimentType, ['Session' num2str(session)], ...
                                    ['MotorNoise-' num2str(experimentNumber(ii))]);
    elseif strcmp(experimentType, 'PerceptNoise')
        dataFullPath = fullfile('Data', subjectID{1}, experimentType, ['Session' num2str(session)], ...
                                    ['PerceptNoise-' num2str(experimentNumber(ii))]);                                
    else
        dataFullPath = fullfile('Data', subjectID{1}, experimentType,[experiment num2str(session)], ...
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

stdNoiseExp = unique(dataAll(:,4));
angleDiff = (unique(dataAll(:,1)))';
dataZeroDiff = dataAll(dataAll(:,1) == 0,:);
angleEstimate = dataAll(:, 7);

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
biasAll = angleTrue - angleEstimate;
binaryFeedbackAll = dataAll(:,8);  
binaryDecisionAll = dataAll(:,6); 
indicatorConsistent = ones(size(binaryDecisionAll));
if includeIncongruentTrials == 0
    indicatorConsistent = sign(angleDiffEst) == binaryFeedbackAll;
    angleDiffEst(~indicatorConsistent) = NaN;
end
angleDiffEstimate = cell(length(stdNoiseExp), length(angleDiff));
biasEstimate = cell(length(stdNoiseExp), length(angleDiff));
binaryDecisionSelf = cell(length(stdNoiseExp), length(angleDiff));
binaryFeedback = binaryDecisionSelf;
for ii = 1 : length(stdNoiseExp)
    for jj = 1 : length(angleDiff)
        tempBias = biasAll(dataAll(:,4)==stdNoiseExp(ii) ...
                           & dataAll(:,1)==angleDiff(jj), :);
        biasEstimate{ii,jj} = tempBias(:);
        tempangleDiffEstimate = angleDiffEst(dataAll(:,4)==stdNoiseExp(ii) ...
                                            & dataAll(:,1)==angleDiff(jj), :);
        binaryDecisionSelf{ii,jj} = binaryDecisionAll(dataAll(:,4)==stdNoiseExp(ii) ...
                                        & dataAll(:,1)==angleDiff(jj), :);  
        binaryFeedback{ii,jj} = binaryFeedbackAll(dataAll(:,4)==stdNoiseExp(ii) ...
                                        & dataAll(:,1)==angleDiff(jj), :);  
        angleDiffEstimate{ii,jj} = tempangleDiffEstimate;                            
    end
end
binaryDecision = binaryDecisionSelf;  

fractionCW = NaN(length(stdNoiseExp),length(params.barAngleDiff));
nTrialsPerCondition = NaN(length(stdNoiseExp),length(params.barAngleDiff));
if strcmp(experiment, 'DecisionGiven')
    binaryDecisionAll = dataAll(:, 8);
else
    binaryDecisionAll = dataAll(:, 6);
end    
for ii = 1 : length(stdNoiseExp)
    for jj = 1 : length(fractionCW)
        indexAll = dataAll(:,1) == params.barAngleDiff(jj) ...
                    & dataAll(:,4) == stdNoiseExp(ii)...
                    & indicatorConsistent;
        indexCW = dataAll(:,1) == params.barAngleDiff(jj) ...
                    & dataAll(:,4) == stdNoiseExp(ii)...
                    & indicatorConsistent ...
                    & binaryDecisionAll == 1;
        nTrialsPerCondition(ii,jj) = nansum((indexAll));
        fractionCW(ii,jj) = nansum(indexCW)/nansum(indexAll);
    end
end

%% Modeling
dstep = 0.1;
thetaStim = -12:0.1:30; % 
thetaStim = round(thetaStim, -log10(dstep));
pCw = paramsFit(9);
pC = [pCw, 1-pCw]'; % [cw ccw]
pthcw = paramsFit(4);
pthccw = paramsFit(5); % paramsFit(4)

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

pCw = NaN(length(stdSensory), length(thetaStim));
estimateModel_Correct.Xval = cell(length(stdSensory));
estimateModel_Correct.Yval = cell(length(stdSensory));  
estimateTheoryCCW_Correct = NaN(length(stdSensory), length(thetaStim));
estimateTheoryCW_Correct = NaN(length(stdSensory), length(thetaStim));
estimateModel_Incorrect.Xval = cell(length(stdSensory));
estimateModel_Incorrect.Yval = cell(length(stdSensory));  
estimateTheoryCCW_Incorrect = NaN(length(stdSensory), length(thetaStim));
estimateTheoryCW_Incorrect = NaN(length(stdSensory), length(thetaStim));

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
    
    nmr = nmm;
    mr = mm;

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
    if correctType == 1
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
    elseif correctType == 2
        pmrGth = exp(-((MR_th-THmr).^2)./(2*(stdSensory(kk)^2 + stdMemory^2)));
        pmrGth = pmrGth./(repmat(sum(pmrGth,1),nmr,1)); 
        pthGmrChcw = (pmrGth.*repmat(pthGC(1,:),nmr,1))';
        pthGmrChcw = pthGmrChcw./repmat(sum(pthGmrChcw,1),nth,1);
        pthGmrChcw(isnan(pthGmrChcw)) = 0;

        pthGmrChccw = (pmrGth.*repmat(pthGC(2,:),nmr,1))';
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

        % Resample m until we have a sample that is consistent with feedback
        % p(mr|m, theta, Chat)
        MR_m = repmat(mr', 1, nm);
        pmrGmth = exp(-((MR_m-repmat(m, nmr, 1)).^2)./(2*(stdMemory^2))); 

        pmrGmthChcw = pmrGmth;
        pmrGmthChcw(mr < 0, :) = 0;
        % put the tail with all 0 to 1 (deal with small memory noise)
        indZero = sum(pmrGmthChcw, 1) == 0;
        pmrGmthChcw(mr > 0, indZero) = 1;
        pmrGmthChcw = pmrGmthChcw./(repmat(sum(pmrGmthChcw,1),nmr,1));
        pmrGmthChcw(mr > 0, indZero) = 1e-50;

        pmrGmthChccw = pmrGmth;
        pmrGmthChccw(mr > 0, :) = 0;
        % put the tail with all 0 to 1 (deal with small memory noise)
        indZero = sum(pmrGmthChccw, 1) == 0;
        pmrGmthChccw(mr < 0, indZero) = 1;        
        pmrGmthChccw = pmrGmthChccw./(repmat(sum(pmrGmthChccw,1),nmr,1));
        pmrGmthChccw(mr < 0, indZero) = 1e-50;

        % Marginalize over m that lead to cw/ccw decision to compute likelihood p(mr|theta, Chat)
        pmGthChcw = pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim));
        pmrGthChcw = pmrGmthChcw * pmGthChcw;   
        pmrGthChcw = pmrGthChcw ./ (repmat(sum(pmrGthChcw,1),nmr,1)); 
        pmrGthChcw(isnan(pmrGthChcw)) = 0;

        pmGthChccw = pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim));
        pmrGthChccw = pmrGmthChccw * pmGthChccw;
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
    end
    
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
    
    pthhGth_correct = pthhGthChcw.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw.*repmat(PChGtheta_lapse(2,:),nth,1);
    pthhAndth_correct = pthhGth_correct;
    
    pCw(kk, :) = PChGtheta_lapse_new(1,:);
    estimateModel_Correct.Xval{kk} = th;
    estimateModel_Correct.Yval{kk} = pthhAndth_correct;
    estimateTheoryCCW_Correct(kk,:) = mthhGthChccw_correct;
    estimateTheoryCW_Correct(kk,:) = mthhGthChcw_correct;
    dummyVar = estimateTheoryCW_Correct;
    
    %% Incorrect trials  
    if incorrectType == 1
        pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemory(kk)^2))); % p(mm|th) = N(th, sm^2 + smm^2)
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
        pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemory(kk)^2)); 
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
    elseif incorrectType == 2
        pthhGthChcw = pthhGthChcw_Incorrect;
        pthhGthChccw = pthhGthChccw_Incorrect;
    elseif incorrectType == 3
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
    elseif incorrectType == 4        
        % Likelihood function the same as correct decision p(mm|th) = N(th, sm^2 + smm^2)           
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

        % Resample m until we have a sample that is consistent with feedback
        % p(mr|m, theta, Chat)
        MR_m = repmat(mr', 1, nm);
        pmrGmth = exp(-((MR_m-repmat(m, nmr, 1)).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); 

        pmrGmthChcw = pmrGmth;
        pmrGmthChcw(mr > 0, :) = 0;
        pmrGmthChcw = pmrGmthChcw./(repmat(sum(pmrGmthChcw,1),nmr,1));

        pmrGmthChccw = pmrGmth;
        pmrGmthChccw(mr < 0, :) = 0;
        pmrGmthChccw = pmrGmthChccw./(repmat(sum(pmrGmthChccw,1),nmr,1));

        % Marginalize over m that lead to cw/ccw decision to compute likelihood p(mr|theta, Chat)
        pmGthChcw = pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim));
        pmrGthChcw = pmrGmthChcw * pmGthChcw;   
        pmrGthChcw = pmrGthChcw ./ (repmat(sum(pmrGthChcw,1),nmr,1)); 
        pmrGthChcw(isnan(pmrGthChcw)) = 0;

        pmGthChccw = pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim));
        pmrGthChccw = pmrGmthChccw * pmGthChccw;
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
    end
    pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
    pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            
    mthhGthChcw_incorrect = th * pthhGthChcw_norm;
    mthhGthChccw_incorrect = th * pthhGthChccw_norm;
    mthhGthChcw_incorrect(thetaStim > 0) = NaN;
    mthhGthChccw_incorrect(thetaStim < 0) = NaN;

    pthhGth_incorrect = pthhGthChcw.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw.*repmat(PChGtheta_lapse(2,:),nth,1);
    pthhAndth_incorrect = pthhGth_incorrect;
    
    estimateModel_Incorrect.Xval{kk} = th;
    estimateModel_Incorrect.Yval{kk} = pthhAndth_incorrect;
    estimateTheoryCCW_Incorrect(kk,:) = mthhGthChccw_incorrect;
    estimateTheoryCW_Incorrect(kk,:) = mthhGthChcw_incorrect;

end    
    
    
% Analyze
rangeCollapseTheory = round(length(thetaStim)/2);
biasTheory_Correct = NaN(length(stdSensory), rangeCollapseTheory);
biasTheory_Incorrect = NaN(length(stdSensory), rangeCollapseTheory);

for ii = 1 : length(stdSensory)
    tempBiasTheory1 = estimateTheoryCW_Correct(ii, rangeCollapseTheory:end) - thetaStim(rangeCollapseTheory:end);
    tempBiasTheory2 = estimateTheoryCCW_Correct(ii, 1:rangeCollapseTheory) - thetaStim(1:rangeCollapseTheory);
    tempBiasTheory = [tempBiasTheory1;-tempBiasTheory2(:,sort(1:size(tempBiasTheory2,2),'descend'))];  
    biasTheory_Correct(ii,:) = squeeze(nanmean(tempBiasTheory,1));
    
    tempBiasTheory1 = estimateTheoryCW_Incorrect(ii, 1:rangeCollapseTheory) - thetaStim(1:rangeCollapseTheory);
    tempBiasTheory2 = estimateTheoryCCW_Incorrect(ii, rangeCollapseTheory:end) - thetaStim(rangeCollapseTheory:end);
    tempBiasTheory = [tempBiasTheory2;-tempBiasTheory1(:,sort(1:size(tempBiasTheory2,2),'descend'))];  
    biasTheory_Incorrect(ii,:) = squeeze(nanmean(tempBiasTheory,1));    
end

%% Plot the psychometric curves
colorName = {'g', 'r', 'b', 'magenta', 'y'};
figure
set(gcf, 'Position',  [300, 200, 500, 500])
hold on
set(gca, 'FontSize', fontSize)
for ii = 1 : size(pCw, 1)
    plot(thetaStim, pCw(ii, :), 'Color', colorName{ii},  'LineWidth', 5)
    plot(angleDiff, fractionCW(ii, :), 'o', 'MarkerSize', 15, 'Color', colorName{ii}, 'MarkerFaceColor', colorName{ii})
end
xlabel('Stimulus orientation (deg)')
ylabel('Fraction "cw"')
xlim([min(angleDiff) max(angleDiff)])

%% Plot the average
% hAverage = figure;
% figPos = [0.3, 0.2, 0.5, 0.7];
% set(hAverage,'Units','normalized','Position',figPos)
% hold on
% rangeCollapse = round(length(angleDiff)/2);
% legendName = cell(1,length(stdSensory));
% maxYplot = max(params.barAngleDiff)+1;
% maxXplot = max(params.barAngleDiff);
% indexColor = 1;
% hLegend = NaN(1, length(stdSensory));
% boxcarLength = 1;
% 
% % Plot true angle
% biasMean_Correct = NaN(length(stdSensory), rangeCollapse);
% biasStd_Correct = NaN(length(stdSensory), rangeCollapse);
% biasMean_Incorrect = NaN(length(stdSensory), rangeCollapse);
% biasStd_Incorrect = NaN(length(stdSensory), rangeCollapse);
% 
% for ii = 1 : length(stdSensory)
%     
%     tempAngleDiffEstimate = [angleDiffEstimate{ii,:}];
%     tempBinaryDecision = [binaryDecision{ii,:}];
%     tempFeedBack = [binaryFeedback{ii,:}];
%     indexCorrect = tempBinaryDecision == tempFeedBack;    
% 
%     % Correct trials
%     subplot(1, length(stdSensory), 1)
%     hold on
%     set(gca,'FontSize',fontSize)
%     tempAngleDiffEstimateCW = NaN(size(tempAngleDiffEstimate));
%     tempAngleDiffEstimateCW(tempBinaryDecision == 1 & indexCorrect) = tempAngleDiffEstimate(tempBinaryDecision == 1 & indexCorrect);
%     tempAngleDiffEstimateCWAve = squeeze(nanmean(tempAngleDiffEstimateCW,1));
%     tempAngleDiffEstimateCWStd = squeeze(nanstd(tempAngleDiffEstimateCW,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCW),1));
%     tempBias1 = tempAngleDiffEstimateCW(:,rangeCollapse:end) - repmat(angleDiff(rangeCollapse:end),size(tempAngleDiffEstimateCW,1),1);
%     tempBias1 = tempBias1 - repmat(biasMotor(rangeCollapse:end), size(tempBias1,1),1);
%     hShade = shadedErrorBar(angleDiff, tempAngleDiffEstimateCWAve, tempAngleDiffEstimateCWStd,... 
%                          {'Color', colorName{ii}, 'LineWidth', lineWidth}, 1, 0, 0);        
% 
%     tempAngleDiffEstimateCCW = NaN(size(tempAngleDiffEstimate));
%     tempAngleDiffEstimateCCW(tempBinaryDecision == -1 & indexCorrect) = tempAngleDiffEstimate(tempBinaryDecision == -1 & indexCorrect);
%     tempAngleDiffEstimateCCWAve = squeeze(nanmean(tempAngleDiffEstimateCCW,1));
%     tempAngleDiffEstimateCCWStd = squeeze(nanstd(tempAngleDiffEstimateCCW,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCCW),1));
%     tempBias2 = tempAngleDiffEstimateCCW(:,1:rangeCollapse) - repmat(angleDiff(1:rangeCollapse),size(tempAngleDiffEstimateCCW,1),1);
%     tempBias2 = tempBias2 - repmat(biasMotor(1:rangeCollapse), size(tempBias2,1),1);
%     tempBias = [tempBias1;-tempBias2(:,sort(1:size(tempBias2,2),'descend'))];
%     tempBiasMean = nanmean(tempBias,1);
%     bias = [ones(1,boxcarLength-1)*tempBiasMean(1) tempBiasMean ones(1,boxcarLength-1)*tempBiasMean(end)];
%     smoothBias = conv(bias, ones(1,boxcarLength)/boxcarLength, 'full');
%     biasMean_Correct(ii,:) = smoothBias(boxcarLength:boxcarLength+rangeCollapse-1);
%     biasStd_Correct(ii,:) = squeeze(nanstd(tempBias,1))./sqrt(nansum(~isnan(tempBias),1));        
%     shadedErrorBar(angleDiff, tempAngleDiffEstimateCCWAve, tempAngleDiffEstimateCCWStd,... 
%                          {'Color', colorName{ii}, 'LineWidth', lineWidth}, 1, 0, 0); 
%     hLegend(ii) = plot(thetaStim, estimateTheoryCCW_Correct(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);
%     plot(thetaStim, estimateTheoryCW_Correct(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);
% 
%     legendName{ii} = ['Noise: ' num2str(stdSensory(ii))];
%     xlabel('Stimulus orientation (deg)')
%     ylabel('Estimated orientation (deg)')  
%     title('Correct trials');
%     plot([angleDiff(1) angleDiff(end)],  [angleDiff(1) angleDiff(end)], '--k')
%     xlim([-maxXplot maxXplot])
%     ylim([-maxYplot maxYplot])
%     box on
%     text(6, -15, '+/-1SEM', 'FontSize', 15)
% 
%     % Incorrect trials
%     subplot(1, length(stdSensory), 2)
%     hold on
%     set(gca,'FontSize',fontSize)
%     tempAngleDiffEstimateCW = NaN(size(tempAngleDiffEstimate));
%     tempAngleDiffEstimateCW(tempBinaryDecision == 1 & ~indexCorrect) = tempAngleDiffEstimate(tempBinaryDecision == 1 & ~indexCorrect);
%     tempAngleDiffEstimateCWAve = squeeze(nanmean(tempAngleDiffEstimateCW,1));
%     tempAngleDiffEstimateCWStd = squeeze(nanstd(tempAngleDiffEstimateCW,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCW),1));
%     tempBias1 = tempAngleDiffEstimateCW(:, 1:rangeCollapse) - repmat(angleDiff(1:rangeCollapse),size(tempAngleDiffEstimateCCW,1),1);
%     tempBias1 = tempBias1 - repmat(biasMotor(rangeCollapse:end), size(tempBias1,1),1);
%     hShade = shadedErrorBar(angleDiff, tempAngleDiffEstimateCWAve, tempAngleDiffEstimateCWStd,... 
%                          {'Color', colorName{ii}, 'LineWidth', lineWidth}, 1, 0, 0);        
% 
%     tempAngleDiffEstimateCCW = NaN(size(tempAngleDiffEstimate));
%     tempAngleDiffEstimateCCW(tempBinaryDecision == -1 & ~indexCorrect) = tempAngleDiffEstimate(tempBinaryDecision == -1 & ~indexCorrect);
%     tempAngleDiffEstimateCCWAve = squeeze(nanmean(tempAngleDiffEstimateCCW,1));
%     tempAngleDiffEstimateCCWStd = squeeze(nanstd(tempAngleDiffEstimateCCW,1))./sqrt(nansum(~isnan(tempAngleDiffEstimateCCW),1));
%     tempBias2 = tempAngleDiffEstimateCCW(:, rangeCollapse:end) - repmat(angleDiff(rangeCollapse:end),size(tempAngleDiffEstimateCW,1),1);
%     tempBias2 = tempBias2 - repmat(biasMotor(1:rangeCollapse), size(tempBias2,1),1);
%     tempBias = [tempBias2;-tempBias1(:,sort(1:size(tempBias2,2),'descend'))];
%     tempBiasMean = nanmean(tempBias,1);
%     bias = [ones(1,boxcarLength-1)*tempBiasMean(1) tempBiasMean ones(1,boxcarLength-1)*tempBiasMean(end)];
%     smoothBias = conv(bias, ones(1,boxcarLength)/boxcarLength, 'full');
%     biasMean_Incorrect(ii,:) = smoothBias(boxcarLength:boxcarLength+rangeCollapse-1);
%     biasStd_Incorrect(ii,:) = squeeze(nanstd(tempBias,1))./sqrt(nansum(~isnan(tempBias),1));        
%     shadedErrorBar(angleDiff, tempAngleDiffEstimateCCWAve, tempAngleDiffEstimateCCWStd,... 
%                          {'Color', colorName{ii}, 'LineWidth', lineWidth}, 1, 0, 0); 
%     hLegend(ii) = plot(thetaStim, estimateTheoryCCW_Incorrect(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);
%     plot(thetaStim, estimateTheoryCW_Incorrect(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);
% 
%     legendName{ii} = ['Noise: ' num2str(stdSensory(ii))];
%     xlabel('Stimulus orientation (deg)')
%     title('Incorrect trials');
%     plot([angleDiff(1) angleDiff(end)],  [angleDiff(1) angleDiff(end)], '--k')
%     xlim([-maxXplot maxXplot])
%     ylim([-maxYplot maxYplot])
%     box on
%     text(6, -15, '+/-1SEM', 'FontSize', 15)    
% end
% legend(hLegend, legendName, 'Location', 'NorthWest')
% 
% % Plot the bias
% hBias = figure;
% hold on
% set(gca,'FontSize',fontSize)
% for ii = 1 : length(stdSensory)
%     % Correct trials
%     subplot(1, length(stdSensory), 1)
%     hold on
%     hShade = shadedErrorBar(angleDiff(rangeCollapse:end), biasMean_Correct(ii,:), biasStd_Correct(ii,:),... 
%                          {'Color', colorName{ii}, 'LineWidth', lineWidth},1,0,0); 
%     hLegend(ii) = plot(thetaStim(rangeCollapseTheory:end), biasTheory_Correct(ii,:), ... 
%                     'Color', colorName{ii}, 'LineWidth', lineWidth);
%     legendName{ii} = ['noise = ' num2str(stdSensory(ii))];
%     xlabel('Absolute true angle (degree)')
%     ylabel('Bias (degree)')  
%     title('Correct trials');
%     xlim([0 maxXplot])
%     ylim([-10 15])
%     set(gca, 'YTick', [-10:5:15]) 
%     plot([0 21], [0 0], 'k--', 'LineWidth', lineWidth)
%     
%     % Incorrect trials
%     subplot(1, length(stdSensory), 2)
%     hold on
%     hShade = shadedErrorBar(angleDiff(rangeCollapse:end), biasMean_Incorrect(ii,:), biasStd_Correct(ii,:),... 
%                          {'Color', colorName{ii}, 'LineWidth', lineWidth},1,0,0); 
%     hLegend(ii) = plot(thetaStim(rangeCollapseTheory:end), biasTheory_Incorrect(ii,:), ... 
%                     'Color', colorName{ii}, 'LineWidth', lineWidth);
%     legendName{ii} = ['noise = ' num2str(stdSensory(ii))];
%     xlabel('Absolute true angle (degree)')
%     ylabel('Bias (degree)')  
%     title('Incorrect trials');
%     xlim([0 maxXplot])
%     ylim([-10 15])
%     set(gca, 'YTick', [-10:5:15]) 
%     plot([0 21], [0 0], 'k--', 'LineWidth', lineWidth)
%    
% end

%% Extract the collapsed data
correctIndex = 1;
[~, ~, ~, estimateDataCorrect, angleDiff, ~] = dataForFitting(subjectID, includeIncongruentTrials, correctIndex);
correctIndex = 0;
[~, ~, ~, estimateDataIncorrect, ~, ~] = dataForFitting(subjectID, includeIncongruentTrials, correctIndex);


%% Plot the estimate distribution
% Correct trials
figure
counter = 1;
binCenter = linspace(-45, 45, 33);
tempEstimateModelX = estimateModel_Correct.Xval{1};
for kk = 1 : length(estimateModel_Correct.Yval)
    estimateModelY = estimateModel_Correct.Yval{kk};
    estimateModelY = estimateModelY(:, ismember(thetaStim, angleDiff));
    for jj = 1 : size(estimateModelY, 2)            
        tempEstimateModelY = estimateModelY(:, jj); 
        tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
        
        tempestimateDataX = estimateDataCorrect{kk,jj};
        
        subplot(length(estimateModel_Correct.Yval), size(estimateModelY, 2), counter)
        hold on
        binCount = hist(tempestimateDataX, binCenter);
        binCount = binCount / (sum(binCount) * diff(binCenter(1:2)));
        bar(binCenter, binCount, 1);
 
        plot(tempEstimateModelX, tempEstimateModelY,'r','LineWidth',3)
        plot([0 0], [0 0.11], 'k--')
        xlim([-47 47])
        ylim([0 0.13])
        set(gca,'YTickLabel',[],'YTick',[])
        if kk == 1
            title(num2str(angleDiff(jj)))
        end
        box off
        
        counter = counter + 1;
    end
end

% Incorrect trials
figure
counter = 1;
binCenter = linspace(-45, 45, 33);
tempEstimateModelX = estimateModel_Incorrect.Xval{1};
for kk = 1 : length(estimateModel_Incorrect.Yval)
    estimateModelY = estimateModel_Incorrect.Yval{kk};
    estimateModelY = estimateModelY(:, ismember(thetaStim, angleDiff));
    for jj = 1 : size(estimateModelY, 2)            
        tempEstimateModelY = estimateModelY(:, jj); 
        tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
        
        tempestimateDataX = estimateDataIncorrect{kk,jj};
        
        subplot(length(estimateModel_Correct.Yval), size(estimateModelY, 2), counter)
        hold on
        binCount = hist(tempestimateDataX, binCenter);
        binCount = binCount / (sum(binCount) * diff(binCenter(1:2)));
        bar(binCenter, binCount, 1);
 
        plot(tempEstimateModelX, tempEstimateModelY,'r','LineWidth',3)
        plot([0 0], [0 0.12], 'k--')
        xlim([-47 47])
        ylim([0 0.15])
        set(gca,'YTickLabel',[],'YTick',[])
        if kk == 1
            title(num2str(angleDiff(jj)))
        end
        box off
        
        counter = counter + 1;
    end
end
