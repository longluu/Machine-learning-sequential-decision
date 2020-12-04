%%%%%%%%%%%%%%%%%%%%%%% Compute the mean estimates of the models %%%%%%%%%%%%%%%%%%%%%%

%% Compute the subject LLH
subjectIDAll = {'ll', 'an', 'ep', 'jp', 'kc'};
estimate_FlipEst = NaN(length(subjectIDAll), 2, 8);
estimate_FlipDecision_2a1 = NaN(length(subjectIDAll), 2, 8);
estimate_LHboundary_1c = NaN(length(subjectIDAll), 2, 8);
estimate_LHestimate_1d = NaN(length(subjectIDAll), 2, 8);
estimate_Resample_2b2 = NaN(length(subjectIDAll), 2, 8);
estimate_Prior = NaN(length(subjectIDAll), 2, 8);
estimate_Data = NaN(length(subjectIDAll), 2, 8);

% % No resample for correct trials
% paramsAllSubject = [2.6500    6.0895           0.0000     22.2852     1.6506   0.9414    2.0976;
%                     3.0023    9.7384           0.0000     34.4053     0.0615   0.9480    3.1069;
%                     4.6136   10.4165           0.0000     29.8375     0.1325   0.9940    3.8106;
%                     7.7094   11.9114           0.0000     55.7419     0.0083   0.2850    3.8551;
%                     5.1033   10.3703           0.0000     46.6421     4.7921   0.8187    3.3313];

% Resample for correct trials
paramsAllSubject = [2.7354    6.1605           0.0000     21.7403     1.5041   0.9972    2.0976;
                    2.7818    9.6676           0.0000     34.4606     0.1772   0.9475    3.1069;
                    3.9181   10.3898           0.0000     29.4658     1.1074   0.9985    3.8106;
                    7.8828   12.2344           0.0000     54.5629     0.1530   0.6860    3.8551;
                    4.2004    9.4129           0.0000     46.4438     6.0270   0.8464    3.3313];

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

    [~, ~, ~, estimateData, ~, ~] = dataForFitting(subjectID, 0, 0);

    % Compute meanEstimate
    meanEstimate = NaN(size(estimateData));
    for ii = 1 : size(estimateData, 1)
        for jj = 1 : size(estimateData, 2)
            tempEst = abs(estimateData{ii, jj});
            if sum(~isnan(tempEst)) > 5
                meanEstimate(ii, jj) = nanmean(tempEst);
            end
        end
    end
    estimate_Data(nn, :, :) = meanEstimate;
    
    %% Compute the mean estimate of models
    flagSC = 1; % 1: self-conditioned model
               % 0: standard Bayes
    includeIncongruentTrials = 0;
    dstep = 0.1;
    paramsAll = paramsAllSubject(nn, :);
    lapseRate = paramsAll(3);

    % stimulus orientation
    thetaStim = angleDiff; % 
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
    pthccw = -paramsAll(4); % paramsAll(4)

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
        
        %% Prior model
        % Compute the estimate
        pthhGthChcw = repmat(normpdf(th', pthccw/2, stdMotor), 1, length(thetaStim));
        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);   
        pthhGthChcw = pthhGthChcw  .* repmat(PChGtheta_lapse(1,:),nth,1);

        pthhGthChccw = repmat(normpdf(th', pthcw/2, stdMotor), 1, length(thetaStim)) .* repmat(PChGtheta_lapse(2,:),nth,1); 
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1); 
        pthhGthChccw =  pthhGthChccw .* repmat(PChGtheta_lapse(2,:),nth,1); 
           

        if includeIncongruentTrials == 0
            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'<= 0, :) = 0;
        end

        % remove 'correct' trials
        pthhGthChccw(:, thetaStim < 0) = 0;
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
        mthhGthChccw= th * pthhGthChccw_norm;
        estimate_Prior(nn, kk, :) = mthhGthChccw(thetaStim >= 0);       

        %% Model 2a1 (Flip decision, Resample: memory, No rejection)
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
        pthhGthChccw(:, thetaStim < 0) = 0;
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
        mthhGthChccw= th * pthhGthChccw_norm;
        estimate_FlipDecision_2a1(nn, kk, :) = mthhGthChccw(thetaStim >= 0);        

        %% Incorrect type 3 (Resample)
        % Likelihood centers on mr, variance: sum of sensory and memory           
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
        end

        % remove 'correct' trials
        pthhGthChccw(:, thetaStim < 0) = 0;
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);  
        mthhGthChccw= th * pthhGthChccw_norm;
        estimate_Resample_2b2(nn, kk, :) = mthhGthChccw(thetaStim >= 0);        
                
        %% Likelihood center at boundary
        % Compute the estimate
        pthGmm = normpdf(th, 0, sqrt(stdSensory(kk)^2 + stdMemory^2));
        pthGmmChcw = pthGmm;
        pthGmmChcw(th<0) = 0;
        pthGmmChcw = pthGmmChcw / sum(pthGmmChcw);
        thhChcw = th * pthGmmChcw';
        
        pthhGthChccw = repmat(normpdf(th', thhChcw, stdMotor), 1, length(thetaStim)) .* repmat(PChGtheta_lapse(2,:),nth,1); 
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1); 
        pthhGthChccw =  pthhGthChccw .* repmat(PChGtheta_lapse(2,:),nth,1); 
           

        if includeIncongruentTrials == 0
            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th'<= 0, :) = 0;
        end

        % remove 'correct' trials
        pthhGthChccw(:, thetaStim < 0) = 0; 
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1); 
        mthhGthChccw= th * pthhGthChccw_norm;
        estimate_LHboundary_1c(nn, kk, :) = mthhGthChccw(thetaStim >= 0); 
        
        %% Likelihood center at estimate
        % Inference: p(thetaHat|mr, cHat) = N(th, sm^2 + smm^2)           
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


        % Get the distribution of LH center p(mr| theta, Chat) = p(thetaHat_temp|theta, Chat)
        % memory noise
        pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemory^2)); 
        pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1)); 

        a = 1./gradient(EthChcw,dstep);
        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChcw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim)));
        b = repmat(a',1,length(thetaStim)) .* pmmGthChcw(indKeepCw, :);        
        pmrGthChcw = interp1(EthChcw,b,th,'linear','extrap');
        pmrGthChcw(pmrGthChcw < 0) = 0; 

        a = 1./gradient(EthChccw,dstep);
        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChccw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim)));        
        b = repmat(a',1,length(thetaStim)) .* pmmGthChccw(indKeepCcw, :);        
        pmrGthChccw = interp1(EthChccw,b,th,'linear','extrap');
        pmrGthChccw(pmrGthChccw < 0) = 0; 

        % Marginalize over the LH center to get the predictive distribution p(thetaHat|theta, Chat)
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
        end

        % remove 'correct' trials
        pthhGthChccw(:, thetaStim < 0) = 0; 
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1); 
        mthhGthChccw= th * pthhGthChccw_norm;
        estimate_LHestimate_1d(nn, kk, :) = mthhGthChccw(thetaStim >= 0); 
            
    end
end

%% Plot the estimates
estimate_FlipDecisionAll = estimate_FlipDecision_2a1;
estimate_ResampleAll = estimate_Resample_2b2;
estimate_PriorAll = estimate_Prior;
estimate_LHboundaryAll = estimate_LHboundary_1c;
estimate_LHestimateAll = estimate_LHestimate_1d;
estimate_DataAll = estimate_Data;

estimate_FlipDecision_2a1 = estimate_FlipDecision_2a1(:);
estimate_Resample_2b2 = estimate_Resample_2b2(:);
estimate_LHboundary_1c = estimate_LHboundary_1c(:);
estimate_LHestimate_1d = estimate_LHestimate_1d(:);
estimate_Prior = estimate_Prior(:);
estimate_Data = estimate_Data(:);

indExclude = isnan(estimate_Data);
estimate_FlipDecision_2a1(indExclude) = [];
estimate_Resample_2b2(indExclude) = [];
estimate_LHboundary_1c(indExclude) = [];
estimate_LHestimate_1d(indExclude) = [];
estimate_Prior(indExclude) = [];
estimate_Data(indExclude) = [];

figure
maxPlot = max([estimate_FlipDecision_2a1; estimate_Prior; estimate_Resample_2b2; estimate_LHboundary_1c; estimate_LHestimate_1d; estimate_Data]) + 1;
minPlot = min([estimate_FlipDecision_2a1; estimate_Prior; estimate_Resample_2b2; estimate_LHboundary_1c; estimate_LHestimate_1d; estimate_Data]) - 1;
colorName = {'SlateGray', 'DarkMagenta', 'DarkOrange', 'DarkGoldenRod', 'SpringGreen',...
            'Teal', 'DodgerBlue', 'Navy'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

subplot(1, 5, 1)
hold on
set(gca, 'FontSize', 12)
for nn = 1 : length(subjectIDAll)
    plot(squeeze(estimate_DataAll(nn, 1, :)), squeeze(estimate_PriorAll(nn, 1, :)), 'o', 'Color', colorIndex(nn, :))
    plot(squeeze(estimate_DataAll(nn, 2, :)), squeeze(estimate_PriorAll(nn, 2, :)), 'x', 'Color', colorIndex(nn, :))
end
plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
xlabel('Mean estimate - data (deg)')
ylabel('Mean estimate - model (deg)')
r = round(corr(estimate_Data, estimate_Prior, 'type', 'Pearson'), 2);
MSE = round(sum((estimate_Data - estimate_Prior).^2) / length(estimate_Data), 1);
title (['Prior, r: ' num2str(r) ', MSE: ' num2str(MSE)])

subplot(1, 5, 2)
hold on
set(gca, 'FontSize', 12)
for nn = 1 : length(subjectIDAll)
    plot(squeeze(estimate_DataAll(nn, 1, :)), squeeze(estimate_FlipDecisionAll(nn, 1, :)), 'o', 'Color', colorIndex(nn, :))
    plot(squeeze(estimate_DataAll(nn, 2, :)), squeeze(estimate_FlipDecisionAll(nn, 2, :)), 'x', 'Color', colorIndex(nn, :))
end
plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
xlabel('Mean estimate - data (deg)')
ylabel('Mean estimate - model (deg)')
r = round(corr(estimate_Data, estimate_FlipDecision_2a1, 'type', 'Pearson'), 2);
MSE = round(sum((estimate_Data - estimate_FlipDecision_2a1).^2) / length(estimate_Data), 1);
title (['Flip Decision, r: ' num2str(r) ', MSE: ' num2str(MSE)])

subplot(1, 5, 3)
hold on
set(gca, 'FontSize', 12)
for nn = 1 : length(subjectIDAll)
    plot(squeeze(estimate_DataAll(nn, 1, :)), squeeze(estimate_ResampleAll(nn, 1, :)), 'o', 'Color', colorIndex(nn, :))
    plot(squeeze(estimate_DataAll(nn, 2, :)), squeeze(estimate_ResampleAll(nn, 2, :)), 'x', 'Color', colorIndex(nn, :))
end
plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
xlabel('Mean estimate - data (deg)')
ylabel('Mean estimate - model (deg)')
r = round(corr(estimate_Data, estimate_Resample_2b2, 'type', 'Pearson'), 2);
MSE = round(sum((estimate_Data - estimate_Resample_2b2).^2) / length(estimate_Data), 1);
title (['Resample, r: ' num2str(r) ', MSE: ' num2str(MSE)])

subplot(1, 5, 4)
hold on
set(gca, 'FontSize', 12)
for nn = 1 : length(subjectIDAll)
    plot(squeeze(estimate_DataAll(nn, 1, :)), squeeze(estimate_LHboundaryAll(nn, 1, :)), 'o', 'Color', colorIndex(nn, :))
    plot(squeeze(estimate_DataAll(nn, 2, :)), squeeze(estimate_LHboundaryAll(nn, 2, :)), 'x', 'Color', colorIndex(nn, :))
end
plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
xlabel('Mean estimate - data (deg)')
ylabel('Mean estimate - model (deg)')
r = round(corr(estimate_Data, estimate_LHboundary_1c, 'type', 'Pearson'), 2);
MSE = round(sum((estimate_Data - estimate_LHboundary_1c).^2) / length(estimate_Data), 1);
title (['LH at boundary, r: ' num2str(r) ', MSE: ' num2str(MSE)])

subplot(1, 5, 5)
hold on
set(gca, 'FontSize', 12)
for nn = 1 : length(subjectIDAll)
    plot(squeeze(estimate_DataAll(nn, 1, :)), squeeze(estimate_LHestimateAll(nn, 1, :)), 'o', 'Color', colorIndex(nn, :))
    plot(squeeze(estimate_DataAll(nn, 2, :)), squeeze(estimate_LHestimateAll(nn, 2, :)), 'x', 'Color', colorIndex(nn, :))
end
plot([minPlot maxPlot], [minPlot maxPlot], 'k--')
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
xlabel('Mean estimate - data (deg)')
ylabel('Mean estimate - model (deg)')
r = round(corr(estimate_Data, estimate_LHestimate_1d, 'type', 'Pearson'), 2);
MSE = round(sum((estimate_Data - estimate_LHestimate_1d).^2) / length(estimate_Data), 1);
title (['LH at estimate, r: ' num2str(r) ', MSE: ' num2str(MSE)])