%%%%%%%%%%%%%%%%%%%%%%% Compute the LLH of the models %%%%%%%%%%%%%%%%%%%%%%

%% Compute the subject LLH
subjectIDAll = {'ll', 'an', 'ep', 'jp', 'kc'};
logLH_AllModel = NaN(length(subjectIDAll), 4);
logLH_LowNoiseModel = NaN(length(subjectIDAll), 4);
logLH_Data = NaN(length(subjectIDAll), 1);

paramsAllSubject = [2.5019    5.9805           0.0000     27.5778   -24.6865   1.9712    2.0976;
                    3.0195    9.4765           0.0000     34.6523   -38.0195   0.0456    3.1069;
                    4.3879    9.9848           0.0000     31.4439   -34.8956   0.1682    3.8106;
                    7.8494   11.5576           0.0000     53.6605   -59.1985   0.0058    3.8551
                    5.1500   10.0831           0.0000     44.0107   -56.0947   4.3459    3.3313];

for nn = 1 : length(subjectIDAll)
    subjectID = subjectIDAll{nn};
    if strcmp(subjectID, 'average')
        experimentNumber = 1:5;
    else
        experimentNumber = 1;    
    end
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
    angleEstimate = dataAll(:, 7);

    [~, ~, ~, estimateData, ~, ~] = dataForFittingNew(subjectID, 0, 0);

    % Compute LLH
    logLH = 0;
    binCenter = linspace(-45, 45, 33);

    for ii = 1 : size(estimateData, 1)
        for jj = 1 : size(estimateData, 2)
            binCount = hist(estimateData{ii, jj}, binCenter);
            pEmpirical = binCount / sum(binCount);
            logLH = binCount .* log(pEmpirical);
            logLH = nansum(logLH);
            if ~isnan(logLH)
                logLH = logLH -  logLH;
            end
        end
    end
    logLH_Data(nn) = logLH;
    
    %% Compute the mean LLH of models
    flagSC = 1; % 1: self-conditioned model
               % 0: standard Bayes
    includeIncongruentTrials = 0;
    incorrectType = 2; % 1: flip the decision bit
                       % 2: flip the estimates
                       % 3: resample the measurement mm until getting a consistent sample

    dstep = 0.1;
    paramsAll = paramsAllSubject(nn, :);
    lapseRate = paramsAll(3);

    % stimulus orientation
    thetaStim = angleDiff; % 
    thetaStim = round(thetaStim, -log10(dstep));

    % sensory noise
    stdSensory = paramsAll(1:2);

    % memory recall noise
    stdMemory = paramsAll(6);
    stdMemoryIncorrect = sqrt(stdMemory^2 + 0^2);

    % motor noise;
    stdMotor = paramsAll(7);

    % priors
    smoothFactor = 0;



    % LOOP - noise levels
    pC = [0.5, 0.5]'; % [cw ccw]
    pthcw = paramsAll(4);
    pthccw = paramsAll(5); % paramsAll(4)

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
    logLH_Model = NaN(2, size(logLH_AllModel, 2));

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
        pthhANDth_incorrect = 2 * pthhANDth_incorrect / sum(pthhANDth_incorrect(:));

        % compute LLH
        logLikelihoodEstimate = 0;
        tempEstimateModelX = th;
        epsZero = 10^(-10);

        for jj = 1 : length(thetaStim)
            tempEstimateModelY = pthhANDth_incorrect(:, jj); 
            tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
            tempestimateDataX = estimateData{kk,jj};
            tempestimateDataX(tempestimateDataX < min(tempEstimateModelX)) = min(tempEstimateModelX);
            tempestimateDataX(tempestimateDataX > max(tempEstimateModelX)) = max(tempEstimateModelX);               
            tempLikelihood = interp1(tempEstimateModelX, tempEstimateModelY, tempestimateDataX, 'pchip');
            tempLikelihood(tempLikelihood == 0) = epsZero; 
            tempLikelihood(isnan(tempLikelihood)) = [];            
            logLikelihoodEstimate = logLikelihoodEstimate + nansum(log(tempLikelihood));             
        end
        logLH_Model(kk, 1) = logLikelihoodEstimate;

        %% Incorrect type 2
        pthhGthChcw = pthhGthChcw_Incorrect;
        pthhGthChccw = pthhGthChccw_Incorrect;
        pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
        pthhGthChcw_norm(isnan(pthhGthChcw_norm)) = 0;    
        pthhGthChccw_norm(isnan(pthhGthChccw_norm)) = 0;
        pthhANDth_incorrect = pthhGthChcw_norm.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw_norm.*repmat(PChGtheta_lapse(2,:),nth,1);
        pthhANDth_incorrect(:, thetaStim == 0) = pthhANDth_incorrect(:, thetaStim == 0)/2;
        pthhANDth_incorrect = 2 * pthhANDth_incorrect / sum(pthhANDth_incorrect(:));

        % compute LLH
        logLikelihoodEstimate = 0;
        tempEstimateModelX = th;
        epsZero = 10^(-10);

        for jj = 1 : length(thetaStim)
            tempEstimateModelY = pthhANDth_incorrect(:, jj); 
            tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
            tempestimateDataX = estimateData{kk,jj};
            tempestimateDataX(tempestimateDataX < min(tempEstimateModelX)) = min(tempEstimateModelX);
            tempestimateDataX(tempestimateDataX > max(tempEstimateModelX)) = max(tempEstimateModelX);               
            tempLikelihood = interp1(tempEstimateModelX, tempEstimateModelY, tempestimateDataX, 'pchip');
            tempLikelihood(tempLikelihood == 0) = epsZero; 
            tempLikelihood(isnan(tempLikelihood)) = [];            
            logLikelihoodEstimate = logLikelihoodEstimate + nansum(log(tempLikelihood));             
        end
        logLH_Model(kk, 2) = logLikelihoodEstimate;

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
        pthhANDth_incorrect = 2 * pthhANDth_incorrect / sum(pthhANDth_incorrect(:));

        % compute LLH
        logLikelihoodEstimate = 0;
        tempEstimateModelX = th;
        epsZero = 10^(-10);

        for jj = 1 : length(thetaStim)
            tempEstimateModelY = pthhANDth_incorrect(:, jj); 
            tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
            tempestimateDataX = estimateData{kk,jj};
            tempestimateDataX(tempestimateDataX < min(tempEstimateModelX)) = min(tempEstimateModelX);
            tempestimateDataX(tempestimateDataX > max(tempEstimateModelX)) = max(tempEstimateModelX);               
            tempLikelihood = interp1(tempEstimateModelX, tempEstimateModelY, tempestimateDataX, 'pchip');
            tempLikelihood(tempLikelihood == 0) = epsZero; 
            tempLikelihood(isnan(tempLikelihood)) = [];            
            logLikelihoodEstimate = logLikelihoodEstimate + nansum(log(tempLikelihood));             
        end
        logLH_Model(kk, 3) = logLikelihoodEstimate;    
        
        %% Incorrect type 4
        pthhGthChcw = repmat(normpdf(th', pthccw/2, stdMotor), 1, length(thetaStim));
        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);   
        pthhGthChcw = pthhGthChcw  .* repmat(PChGtheta_lapse(1,:),nth,1);

        pthhGthChccw = repmat(normpdf(th', pthcw/2, stdMotor), 1, length(thetaStim)) .* repmat(PChGtheta_lapse(2,:),nth,1); 
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1); 
        pthhGthChccw =  pthhGthChccw .* repmat(PChGtheta_lapse(2,:),nth,1); 
           

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
        pthhANDth_incorrect = 2 * pthhANDth_incorrect / sum(pthhANDth_incorrect(:));

        % compute LLH
        logLikelihoodEstimate = 0;
        tempEstimateModelX = th;
        epsZero = 10^(-10);

        for jj = 1 : length(thetaStim)
            tempEstimateModelY = pthhANDth_incorrect(:, jj); 
            tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
            tempestimateDataX = estimateData{kk,jj};
            tempestimateDataX(tempestimateDataX < min(tempEstimateModelX)) = min(tempEstimateModelX);
            tempestimateDataX(tempestimateDataX > max(tempEstimateModelX)) = max(tempEstimateModelX);               
            tempLikelihood = interp1(tempEstimateModelX, tempEstimateModelY, tempestimateDataX, 'pchip');
            tempLikelihood(tempLikelihood == 0) = epsZero; 
            tempLikelihood(isnan(tempLikelihood)) = [];            
            logLikelihoodEstimate = logLikelihoodEstimate + nansum(log(tempLikelihood));             
        end
        logLH_Model(kk, 4) = logLikelihoodEstimate;            
    end
    logLH_AllModel(nn, :) = nansum(logLH_Model, 1);
    logLH_LowNoiseModel(nn, :) = logLH_Model(1,:);
end

%% Plot the LLH
% Normalize the -LLH of models such that the oracle model is 0 and the
% prior only model is 1
upperBound = repmat(logLH_Data, 1, size(logLH_AllModel, 2));
lowerBound = repmat(logLH_AllModel(:, end), 1, size(logLH_AllModel, 2));
logLH_AllModel = (logLH_AllModel - lowerBound) ./ (upperBound - lowerBound);
logLH_AllModel = logLH_AllModel(:, 1:3);
maxY = 1;
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

hFig = figure;
[~, hPanel] = errorbar_groups(logLH_AllModel', zeros(size(logLH_AllModel')), zeros(size(logLH_AllModel')), ...
                'bar_width', 0.6, 'errorbar_width', 0, 'bar_colors', colorIndex, 'FigID', hFig);
ylabel('Log Likelihood (normalized)')
xlabel('Subject') 
set(gca, 'FontSize', 20, 'XTickLabel', {'1'; '2'; '3'; '4'; '5'; 'Average'})
ylim([0 maxY]);
title('Incorrect trials')
legend('Flip decision', 'Flip estimate', 'Resample', 'Location', 'NorthWest')
