%%%%%%%% The sampling Bayesian model of conditioned perception %%%%%%%%%%
% ********** Old version - condition the prior only ********** 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
% 


%% Set the parameters of the model
subjectID = 'average';
flagSC = 1; % 1: self-conditioned model
            % 0: standard Bayes
incorrectType = 3; % 1: flip the decision bit
                   % 2: flip the estimates thetaHat
                   % 3: resample from the measurement, use old noise model for the likelihood
                   % 4: resample from the measurement, use new noise model for the likelihood
            
params = [3.2048    6.4420           0.0000    30   -12     0.0080       2.0976];
stdSensory = params(1:2);
lapseRate = params(3);
stdMemory = params(6);
stdMotor = params(7);                     
smoothFactor = 0;
boundaryCutoff = 0;
thetaStim = [-12:2:0 5:5:30];
pC = [0.5, 0.5]'; % [cw ccw]
pthcw = params(4);
pthccw = params(5);
nTrialPerCondition = 1000;
fontSize = 20;
fixedConsistencyIndex = 1;
marginalizedVersion = 1;
if fixedConsistencyIndex
    consistencyIndex = 1;
end
includeIncongruentTrials = 0;
rangeCollapse = round(length(thetaStim)/2);

%% Modeling
dstep = 0.1;
rangeth = [-60 60];
th = rangeth(1):dstep:rangeth(2);
nth = length(th);

pthGC = zeros(2,nth);
if flagSC
    pthGC(1,:) = TukeyWindowNew([0 pthcw], smoothFactor, th, boundaryCutoff);
    pthGC(2,:) = TukeyWindowNew([pthccw 0], smoothFactor, th, boundaryCutoff);
else
    pth = (TukeyWindow([0 pthcw], 0, smoothFactor, th) + TukeyWindow([pthccw 0], 1, smoothFactor, th))/2;
    pth(th==0) = 0;
    pth(th==0) = max(pth);
    pthGC(1,:) = pth;
    pthGC(2,:) = pth;    
end

pTheta_CCW = pthGC(2,:);
pTheta_CW = pthGC(1,:);
pTheta_CCW = repmat(pTheta_CCW, nTrialPerCondition, 1);
pTheta_CW = repmat(pTheta_CW, nTrialPerCondition, 1);
pTheta_C_new = zeros(nTrialPerCondition, nth);

%%% SAMPLING version
discriminateSim = NaN(nTrialPerCondition, length(stdSensory), length(thetaStim));
estimateSim = NaN(nTrialPerCondition, length(stdSensory), length(thetaStim));
biasEstimateSim = NaN(nTrialPerCondition, length(stdSensory), length(thetaStim));
indexCorrectSim = NaN(nTrialPerCondition, length(stdSensory), length(thetaStim));
if incorrectType == 1
    for jj = 1 : length(stdSensory)
        for kk = 1 : length(thetaStim)   
            % Draw a sample and measure
            m = normrnd(thetaStim(kk), stdSensory(jj),[nTrialPerCondition 1]);

            % Discrimination (Inference from p(C|m))
            pM_Theta = normpdf(repmat(m,1,nth), repmat(th,nTrialPerCondition,1), stdSensory(jj));
            pCCW_m = pC(1) .* squeeze(sum((pTheta_CCW .* pM_Theta),2)); 
            pCW_m = pC(2) .* squeeze(sum((pTheta_CW .* pM_Theta),2)); 
            discriminateSim(pCCW_m >= pCW_m, jj, kk) = -1;
            discriminateSim(pCCW_m < pCW_m, jj, kk) = 1;
            tempDiscriminate = squeeze(discriminateSim(:,jj,kk));

            % Correct feedback
            correctFeedback = sign(thetaStim(kk)) * ones(size(m));
            if thetaStim(kk) < 0
                pTheta_C_new = pTheta_CCW;
            elseif thetaStim(kk) > 0
                pTheta_C_new = pTheta_CW;
            else
                inCCW = round(rand(1,nTrialPerCondition));
                pTheta_C_new(inCCW==1, :) = squeeze(pTheta_CCW(inCCW==1,:));
                pTheta_C_new(inCCW==0, :) = squeeze(pTheta_CW(inCCW==0,:));
                correctFeedback(inCCW==1) = -1;
                correctFeedback(inCCW==0) = 1;
            end

            % Estimation
            Mm = normrnd(m, stdMemory);                    
            pMm_Theta = normpdf(repmat(Mm,1,nth), repmat(th,nTrialPerCondition,1), sqrt(stdSensory(jj)^2 + stdMemory^2));
            pTheta_Mm_c = pTheta_C_new .* pMm_Theta;
            pMm_c = trapz(th, pTheta_Mm_c, 2);
            pTheta_m_Norm = pTheta_Mm_c ./ repmat(pMm_c, 1, nth);
            thetaEstimate = trapz(th, repmat(th, nTrialPerCondition, 1) .* pTheta_m_Norm, 2);
            thetaResponse = normrnd(thetaEstimate, stdMotor);
            if includeIncongruentTrials == 0
                indexExcludeIncongruent = sign(thetaResponse) ~= sign(thetaEstimate);
                thetaResponse(indexExcludeIncongruent) = NaN;
                discriminateSim(indexExcludeIncongruent, jj, kk) = NaN;
            end
            estimateSim(:, jj, kk) = thetaResponse;
            biasEstimateSim(:, jj, kk) = thetaStim(kk) - thetaResponse;
            indexCorrectSim(:, jj, kk) = (tempDiscriminate == correctFeedback);
        end
    end
elseif incorrectType == 2
    for jj = 1 : length(stdSensory)
        for kk = 1 : length(thetaStim)   
            % Draw a sample and measure
            m = normrnd(thetaStim(kk), stdSensory(jj),[nTrialPerCondition 1]);

            % Discrimination (Inference from p(C|m))
            pM_Theta = normpdf(repmat(m,1,nth), repmat(th,nTrialPerCondition,1), stdSensory(jj));
            pCCW_m = pC(1) .* squeeze(sum((pTheta_CCW .* pM_Theta),2)); 
            pCW_m = pC(2) .* squeeze(sum((pTheta_CW .* pM_Theta),2)); 
            discriminateSim(pCCW_m >= pCW_m, jj, kk) = -1;
            discriminateSim(pCCW_m < pCW_m, jj, kk) = 1;
            tempDiscriminate = squeeze(discriminateSim(:,jj,kk));

            % Correct feedback
            correctFeedback = sign(thetaStim(kk)) * ones(size(m));
            if thetaStim(kk) == 0
                inCCW = round(rand(1,nTrialPerCondition));
                correctFeedback(inCCW==1) = -1;
                correctFeedback(inCCW==0) = 1;
            end
            pTheta_C_new(tempDiscriminate == -1, :) = pTheta_CCW(tempDiscriminate == -1, :);
            pTheta_C_new(tempDiscriminate == 1, :) = pTheta_CW(tempDiscriminate == 1, :);

            % Estimation
            Mm = normrnd(m, stdMemory);                    
            pMm_Theta = normpdf(repmat(Mm,1,nth), repmat(th,nTrialPerCondition,1), sqrt(stdSensory(jj)^2 + stdMemory^2));
            pTheta_Mm_c = pTheta_C_new .* pMm_Theta;
            pMm_c = trapz(th, pTheta_Mm_c, 2);
            pTheta_m_Norm = pTheta_Mm_c ./ repmat(pMm_c, 1, nth);
            thetaEstimate = trapz(th, repmat(th, nTrialPerCondition, 1) .* pTheta_m_Norm, 2);
            
            % Flip the estimate
            thetaEstimate(sign(thetaEstimate) ~= correctFeedback) = -thetaEstimate(sign(thetaEstimate) ~= correctFeedback);
            thetaResponse = normrnd(thetaEstimate, stdMotor);
            if includeIncongruentTrials == 0
                indexExcludeIncongruent = sign(thetaResponse) ~= sign(thetaEstimate);
                thetaResponse(indexExcludeIncongruent) = NaN;
                discriminateSim(indexExcludeIncongruent, jj, kk) = NaN;
            end
            estimateSim(:, jj, kk) = thetaResponse;
            biasEstimateSim(:, jj, kk) = thetaStim(kk) - thetaResponse;
            indexCorrectSim(:, jj, kk) = (tempDiscriminate == correctFeedback);
        end
    end 
elseif incorrectType == 3
    for jj = 1 : length(stdSensory)
        for kk = 1 : length(thetaStim)   
            % Draw a sample and measure
            m = normrnd(thetaStim(kk), stdSensory(jj),[nTrialPerCondition 1]);

            % Discrimination (Inference from p(C|m))
            pM_Theta = normpdf(repmat(m,1,nth), repmat(th,nTrialPerCondition,1), stdSensory(jj));
            pCCW_m = pC(1) .* squeeze(sum((pTheta_CCW .* pM_Theta),2)); 
            pCW_m = pC(2) .* squeeze(sum((pTheta_CW .* pM_Theta),2)); 
            discriminateSim(pCCW_m >= pCW_m, jj, kk) = -1;
            discriminateSim(pCCW_m < pCW_m, jj, kk) = 1;
            tempDiscriminate = squeeze(discriminateSim(:,jj,kk));

            % Correct feedback
            correctFeedback = sign(thetaStim(kk)) * ones(size(m));
            if thetaStim(kk) < 0
                pTheta_C_new = pTheta_CCW;
            elseif thetaStim(kk) > 0
                pTheta_C_new = pTheta_CW;
            else
                inCCW = round(rand(1,nTrialPerCondition));
                pTheta_C_new(inCCW==1, :) = squeeze(pTheta_CCW(inCCW==1,:));
                pTheta_C_new(inCCW==0, :) = squeeze(pTheta_CW(inCCW==0,:));
                correctFeedback(inCCW==1) = -1;
                correctFeedback(inCCW==0) = 1;
            end

            % Memory noise degrades the measurement
            Mm = normrnd(m, stdMemory);                    

            % Resample the incorrect measurement
            Mr = Mm;
            indIncorrect = tempDiscriminate ~= correctFeedback;
            mIncorrect = Mm(indIncorrect);
            if ~isempty(mIncorrect)
                pResample = normpdf(repmat(th, length(mIncorrect), 1), ...
                                    repmat(mIncorrect, 1, nth), sqrt(stdSensory(jj)^2 + stdMemory^2));            
                if thetaStim(kk) < 0
                    pResample(:, th > 0) = 0;
                elseif thetaStim(kk) > 0
                    pResample(:,th < 0) = 0;
                else
                    inCCW = inCCW(indIncorrect);
                    pResample(inCCW==1, th > 0) = 0;
                    pResample(inCCW==0, th < 0) = 0;
                end

                mResampled = NaN(size(mIncorrect));
                for tt = 1 : length(mIncorrect)
                    mResampled(tt) = randpdf(pResample(tt, :), th, [1 1]);
                end
                Mr(indIncorrect) = mResampled;
            end
            
            % Estimation
            pMr_Theta = normpdf(repmat(Mr,1,nth), repmat(th,nTrialPerCondition,1), sqrt(stdSensory(jj)^2 + stdMemory^2));
            pTheta_Mm_c = pTheta_C_new .* pMr_Theta;
            pMm_c = trapz(th, pTheta_Mm_c, 2);
            pTheta_m_Norm = pTheta_Mm_c ./ repmat(pMm_c, 1, nth);
            thetaEstimate = trapz(th, repmat(th, nTrialPerCondition, 1) .* pTheta_m_Norm, 2);
            thetaResponse = normrnd(thetaEstimate, stdMotor);
            if includeIncongruentTrials == 0
                indexExcludeIncongruent = sign(thetaResponse) ~= sign(thetaEstimate);
                thetaResponse(indexExcludeIncongruent) = NaN;
                discriminateSim(indexExcludeIncongruent, jj, kk) = NaN;
            end
            estimateSim(:, jj, kk) = thetaResponse;
            biasEstimateSim(:, jj, kk) = thetaStim(kk) - thetaResponse;
            indexCorrectSim(:, jj, kk) = (tempDiscriminate == correctFeedback);
        end
    end    
end
a = estimateSim(:);
sum(~isnan(a)) / length(a)

%% MARGINALIZED version
if marginalizedVersion == 1
    estimateModelCorrect.Xval = cell(length(stdSensory));
    estimateModelCorrect.Yval = cell(length(stdSensory));  
    estimateTheoryCCW_Correct = NaN(length(stdSensory), length(thetaStim));
    estimateTheoryCW_Correct = NaN(length(stdSensory), length(thetaStim));
    estimateModelIncorrect.Xval = cell(length(stdSensory));
    estimateModelIncorrect.Yval = cell(length(stdSensory));  
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
        % orientation noise p(m|th)
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

        pthhGth_correct = pthhGthChcw.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw.*repmat(PChGtheta_lapse(2,:),nth,1);
        pthhAndth_correct = pthhGth_correct;

        estimateModelCorrect.Xval{kk} = th;
        estimateModelCorrect.Yval{kk} = pthhAndth_correct;
        estimateTheoryCCW_Correct(kk,:) = mthhGthChccw_correct;
        estimateTheoryCW_Correct(kk,:) = mthhGthChcw_correct;

        %% Incorrect trials
        if incorrectType == 1 
            % Generative (forward)
            % orientation noise
            pmGth = exp(-((M-THm).^2)./(2*stdSensory(kk)^2));
            pmGth = pmGth./(repmat(sum(pmGth,1),nm,1)); 

            % Inference
            % 2: estimation
            pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); % p(mm|th) = N(th, sm^2 + smm^2)
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
         
        end
        
        pthhGthChcw_norm = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1);    
        pthhGthChccw_norm = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1); 
        
        mthhGthChcw_incorrect = th * pthhGthChcw_norm;
        mthhGthChccw_incorrect = th * pthhGthChccw_norm;
        mthhGthChcw_incorrect(thetaStim > 0) = NaN;
        mthhGthChccw_incorrect(thetaStim < 0) = NaN;

        pthhGth_incorrect = pthhGthChcw.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw.*repmat(PChGtheta_lapse(2,:),nth,1);
        pthhAndth_incorrect = pthhGth_incorrect;   
        
        estimateModelIncorrect.Xval{kk} = th;
        estimateModelIncorrect.Yval{kk} = pthhAndth_incorrect;
        estimateTheoryCCW_Incorrect(kk,:) = mthhGthChccw_incorrect;
        estimateTheoryCW_Incorrect(kk,:) = mthhGthChcw_incorrect;        
    end
end

% Analyze and plot
estimateResponseSimCCW_correct = NaN(size(estimateSim));
estimateResponseSimCCW_correct(discriminateSim == -1 & indexCorrectSim == 1) = estimateSim(discriminateSim == -1 & indexCorrectSim == 1);
estimateResponseSimCCWAve_correct = squeeze(nanmean(estimateResponseSimCCW_correct, 1));
estimateResponseSimCCWStd_correct = squeeze(nanstd(estimateResponseSimCCW_correct, 1))./(sqrt(squeeze(nansum(~isnan(estimateResponseSimCCW_correct),1))));
estimateResponseSimCW_correct = NaN(size(estimateSim));
estimateResponseSimCW_correct(discriminateSim == 1 & indexCorrectSim == 1) = estimateSim(discriminateSim == 1 & indexCorrectSim == 1);
estimateResponseSimCWAve_correct = squeeze(nanmean(estimateResponseSimCW_correct, 1));
estimateResponseSimCWStd_correct = squeeze(nanstd(estimateResponseSimCW_correct, 1))./(sqrt(squeeze(nansum(~isnan(estimateResponseSimCW_correct),1))));

estimateResponseSimCCW_incorrect = NaN(size(estimateSim));
estimateResponseSimCCW_incorrect(discriminateSim == -1 & indexCorrectSim == 0) = estimateSim(discriminateSim == -1 & indexCorrectSim == 0);
estimateResponseSimCCWAve_incorrect = squeeze(nanmean(estimateResponseSimCCW_incorrect, 1));
estimateResponseSimCCWStd_incorrect = squeeze(nanstd(estimateResponseSimCCW_incorrect, 1))./(sqrt(squeeze(nansum(~isnan(estimateResponseSimCCW_incorrect),1))));
estimateResponseSimCW_incorrect = NaN(size(estimateSim));
estimateResponseSimCW_incorrect(discriminateSim == 1 & indexCorrectSim == 0) = estimateSim(discriminateSim == 1 & indexCorrectSim == 0);
estimateResponseSimCWAve_incorrect = squeeze(nanmean(estimateResponseSimCW_incorrect, 1));
estimateResponseSimCWStd_incorrect = squeeze(nanstd(estimateResponseSimCW_incorrect, 1))./(sqrt(squeeze(nansum(~isnan(estimateResponseSimCW_incorrect),1))));

hAverage = figure;
lineWidth = 4;
figPos = [0.3, 0.2, 0.5, 0.7];
set(hAverage,'Units','normalized','Position',figPos)
hold on
legendName = cell(1,length(stdSensory));
colorName = {'g', 'r', 'b', 'cyan', 'magenta', 'y'};
maxPlot = max(thetaStim);
hLegend = NaN(1, length(stdSensory));
biasMean_correct = NaN(length(stdSensory), rangeCollapse);
biasStd_correct = NaN(length(stdSensory), rangeCollapse);
biasTheory_correct = NaN(length(stdSensory), rangeCollapse);
biasMean_incorrect = NaN(length(stdSensory), rangeCollapse);
biasStd_incorrect = NaN(length(stdSensory), rangeCollapse);
biasTheory_incorrect = NaN(length(stdSensory), rangeCollapse);

for ii = 1 : length(stdSensory)
    % Correct trials
    set(gca,'FontSize',fontSize)
    hShade = shadedErrorBar(thetaStim, estimateResponseSimCCWAve_correct(ii,:), estimateResponseSimCCWStd_correct(ii,:),... 
                         {'Color', colorName{ii}, 'LineWidth', lineWidth},1,0,0);        
    hLegend(ii) = hShade.patch;
    hold on
    shadedErrorBar(thetaStim, estimateResponseSimCWAve_correct(ii, :), estimateResponseSimCWStd_correct(ii, :),... 
                         {'Color', colorName{ii}, 'LineWidth', lineWidth},1,0,0);        
    if  marginalizedVersion
        hLegend(ii) = plot(thetaStim, estimateTheoryCCW_Correct(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);
        plot(thetaStim, estimateTheoryCW_Correct(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);
    end
    legendName{ii} = ['Noise level = ' num2str(stdSensory(ii))];
    xlabel('True angle (degree)')
    ylabel('Angle estimate (degree)')  
    plot([-maxPlot maxPlot], [-maxPlot maxPlot], 'k--', 'LineWidth', 2)
    xlim([-maxPlot maxPlot])
    ylim([-maxPlot maxPlot])
    text(6, -15, '+/-1SEM', 'FontSize', 15)

    % Get the bias
    tempEstimateCW = squeeze(estimateResponseSimCW_correct(:,ii,:));
    tempEstimateCCW = squeeze(estimateResponseSimCCW_correct(:,ii,:));
    tempBias1 = tempEstimateCW(:,rangeCollapse:end) - repmat(thetaStim(rangeCollapse:end),size(tempEstimateCW,1),1);
    tempBias2 = tempEstimateCCW(:,1:rangeCollapse) - repmat(thetaStim(1:rangeCollapse),size(tempEstimateCCW,1),1);
    tempBias = [tempBias1;-tempBias2(:,sort(1:size(tempBias2,2),'descend'))];
    biasMean_correct(ii,:) = squeeze(nanmean(tempBias,1));
    biasStd_correct(ii,:) = squeeze(nanstd(tempBias,1))./sqrt(nansum(~isnan(tempBias),1)); 

    if marginalizedVersion
        tempBiasTheory1 = estimateTheoryCW_Correct(ii,rangeCollapse:end) - thetaStim(rangeCollapse:end);
        tempBiasTheory2 = estimateTheoryCCW_Correct(ii,1:rangeCollapse) - thetaStim(1:rangeCollapse);
        tempBiasTheory = [tempBiasTheory1;-tempBiasTheory2(:,sort(1:size(tempBiasTheory2,2),'descend'))];  
        biasTheory_correct(ii,:) = squeeze(nanmean(tempBiasTheory,1));
    end
    
    % Incorrect trials
    set(gca,'FontSize',fontSize)
    hShade = shadedErrorBar(thetaStim, estimateResponseSimCCWAve_incorrect(ii,:), estimateResponseSimCCWStd_incorrect(ii,:),... 
                         {'Color', colorName{ii}, 'LineWidth', lineWidth},1,0,0);        
    hLegend(ii) = hShade.patch;
    hold on
    shadedErrorBar(thetaStim, estimateResponseSimCWAve_incorrect(ii, :), estimateResponseSimCWStd_incorrect(ii, :),... 
                         {'Color', colorName{ii}, 'LineWidth', lineWidth},1,0,0);        
    if  marginalizedVersion
        hLegend(ii) = plot(thetaStim, estimateTheoryCCW_Incorrect(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);
        plot(thetaStim, estimateTheoryCW_Incorrect(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);
    end
    legendName{ii} = ['Noise level = ' num2str(stdSensory(ii))];
    xlabel('True angle (degree)')
    ylabel('Angle estimate (degree)')  
    plot([-maxPlot maxPlot], [-maxPlot maxPlot], 'k--', 'LineWidth', 2)
    xlim([-maxPlot maxPlot])
    ylim([-maxPlot maxPlot])
    text(6, -15, '+/-1SEM', 'FontSize', 15)

    % Get the bias
    tempEstimateCW = squeeze(estimateResponseSimCW_incorrect(:,ii,:));
    tempEstimateCCW = squeeze(estimateResponseSimCCW_incorrect(:,ii,:));
    tempBias1 = tempEstimateCW(:,rangeCollapse:end) - repmat(thetaStim(rangeCollapse:end),size(tempEstimateCW,1),1);
    tempBias2 = tempEstimateCCW(:,1:rangeCollapse) - repmat(thetaStim(1:rangeCollapse),size(tempEstimateCCW,1),1);
    tempBias = [tempBias1;-tempBias2(:,sort(1:size(tempBias2,2),'descend'))];
    biasMean_incorrect(ii,:) = squeeze(nanmean(tempBias,1));
    biasStd_incorrect(ii,:) = squeeze(nanstd(tempBias,1))./sqrt(nansum(~isnan(tempBias),1)); 

    if marginalizedVersion
        tempBiasTheory1 = estimateTheoryCW_Incorrect(ii,rangeCollapse:end) - thetaStim(rangeCollapse:end);
        tempBiasTheory2 = estimateTheoryCCW_Incorrect(ii,1:rangeCollapse) - thetaStim(1:rangeCollapse);
        tempBiasTheory = [tempBiasTheory1;-tempBiasTheory2(:,sort(1:size(tempBiasTheory2,2),'descend'))];  
        biasTheory_incorrect(ii,:) = squeeze(nanmean(tempBiasTheory,1));
    end    
end
legend(hLegend, legendName, 'Location', 'NorthWest')

% hBiasTheory = figure;
% lineWidth = 4;
% legendName = cell(1,length(stdSensory));
% colorName = {'b', 'g', 'r', 'cyan', 'magenta', 'y'};
% maxPlot = max(thetaStim);
% hLegend = NaN(1, length(stdSensory));
% rangeCollapse = round(length(thetaStim)/2);
% set(gca,'FontSize',fontSize)
% for ii = 1 : length(stdSensory)
%     hold on
%     hTheory = plot(thetaStim(rangeCollapse:end), biasTheory_correct(ii,:), ... 
%                          'Color', colorName{ii}, 'LineWidth', lineWidth);        
%     hLegend(ii) = hTheory;
%         shadedErrorBar(thetaStim(rangeCollapse:end), biasMean_correct(ii,:), biasStd_correct(ii,:),... 
%                              {'Color', colorName{ii}, 'LineWidth', lineWidth},1,0,0);        
%     legendName{ii} = ['Noise = ' num2str(stdSensory(ii))];
%     xlabel('True angle (degree)')
%     ylabel('Bias (degree)')  
% %     title(['Subject ' upper(subjectID) ]);
%     xlim([0 maxPlot])
%     ylim([-10 15])
%     box on
%     text(6, -15, '+/-1SEM', 'FontSize', 15)
% end
% plot([0 maxPlot], [0 0], 'k--', 'LineWidth', lineWidth)

%% Plot the smoothed raw data
hScatter = figure;
figPos = [0.1, 0.2, 0.8, 0.6];
set(hScatter,'Units','normalized','Position',figPos)
hold on

% Correct trials
maxYplot = 40;
yAxis = -maxYplot:1:maxYplot;
lengthXaxis = length(thetaStim(1):thetaStim(end)) + diff(thetaStim(end-1:end)) -1;
maxthetaStim = max(thetaStim);
for ii = 1 : length(stdSensory)
    imageData = zeros(length(yAxis), lengthXaxis);
    counter = 1;
    for jj = 1 :length(thetaStim)
        tempthetaStimEst = squeeze(estimateSim(:,ii,jj));
        tempIndexCorrectSim = squeeze(indexCorrectSim(:,ii,jj));
        tempthetaStimEst(tempIndexCorrectSim == 0) = []; 
        tempthetaStimEst(isnan(tempthetaStimEst)) = [];
        tempthetaStimEst(abs(tempthetaStimEst)>maxYplot) = [];
        
        if thetaStim(jj) < 0
            imageData(:, counter:counter+1) = repmat((hist(tempthetaStimEst, yAxis))', [1 2]); 
            counter = counter + 2;
        else
            imageData(:, counter:counter+4) = repmat((hist(tempthetaStimEst, yAxis))', [1 5]); 
            counter = counter + 5;            
        end
        
    end

    myfilter = fspecial('gaussian', [1 1], 2);
    smoothImage = imfilter(imageData, myfilter, 'replicate');
    imageData = round(smoothImage*255/max(smoothImage(:)));

    subplot(2, length(stdSensory), ii)
    hold on
    tempImage = uint8(imageData);
    [height, width] = size(tempImage);
    tempImage = max(tempImage(:)) - tempImage;
    imagesc(tempImage)
    hold on;
    axis xy;
    colormap('gray');
    
    indYStart = find(yAxis == thetaStim(1));
    indYEnd = find(yAxis == thetaStim(end));
    xDisplay = -10:10:30;
    xTick = round((xDisplay + 12) * lengthXaxis / (angleDiff(end)-angleDiff(1)));
    xTick(xDisplay<=0) = xTick(xDisplay<=0) + 2;

    plot([1 width],[round(height/2) round(height/2)],'k:', 'LineWidth', 1.5);
    plot([15 15],[1 height],'k:', 'LineWidth', 1.5);
    plot([2 lengthXaxis-2],[indYStart indYEnd],'w:', 'LineWidth', 2);
    set(gca, 'ylim', [1 height], 'xlim', [1 width], ...
        'XTick', xTick, 'XTickLabel', num2cell(xDisplay),...
        'YTick', round(linspace(1,height,5)), 'YTickLabel', num2cell(round(linspace(-maxYplot,maxYplot,5))), ...
        'FontSize', 20)
    xlabel('True orientation (degree)')
    ylabel('Estimated orientation (degree)')  
end

% Incorrect trials
maxYplot = 40;
yAxis = -maxYplot:1:maxYplot;
lengthXaxis = length(thetaStim(1):thetaStim(end)) + diff(thetaStim(end-1:end)) -1;
maxthetaStim = max(thetaStim);
for ii = 1 : length(stdSensory)
    imageData = zeros(length(yAxis), lengthXaxis);
    counter = 1;
    for jj = 1 :length(thetaStim)
        tempthetaStimEst = squeeze(estimateSim(:,ii,jj));
        tempIndexCorrectSim = squeeze(indexCorrectSim(:,ii,jj));
        tempthetaStimEst(tempIndexCorrectSim == 1) = []; 
        tempthetaStimEst(isnan(tempthetaStimEst)) = [];
        tempthetaStimEst(abs(tempthetaStimEst)>maxYplot) = [];
        
        if thetaStim(jj) < 0
            imageData(:, counter:counter+1) = repmat((hist(tempthetaStimEst, yAxis))', [1 2]); 
            counter = counter + 2;
        else
            imageData(:, counter:counter+4) = repmat((hist(tempthetaStimEst, yAxis))', [1 5]); 
            counter = counter + 5;            
        end
        
    end

    myfilter = fspecial('gaussian', [1 1], 2);
    smoothImage = imfilter(imageData, myfilter, 'replicate');
    imageData = round(smoothImage*255/max(smoothImage(:)));

    subplot(2, length(stdSensory), ii+2)
    hold on
    tempImage = uint8(imageData);
    [height, width] = size(tempImage);
    tempImage = max(tempImage(:)) - tempImage;
    imagesc(tempImage)
    hold on;
    axis xy;
    colormap('gray');
    
    indYStart = find(yAxis == thetaStim(1));
    indYEnd = find(yAxis == thetaStim(end));
    xDisplay = -10:10:30;
    xTick = round((xDisplay + 12) * lengthXaxis / (angleDiff(end)-angleDiff(1)));
    xTick(xDisplay<=0) = xTick(xDisplay<=0) + 2;

    plot([1 width],[round(height/2) round(height/2)],'k:', 'LineWidth', 1.5);
    plot([15 15],[1 height],'k:', 'LineWidth', 1.5);
    plot([2 lengthXaxis-2],[indYStart indYEnd],'w:', 'LineWidth', 2);
    set(gca, 'ylim', [1 height], 'xlim', [1 width], ...
        'XTick', xTick, 'XTickLabel', num2cell(xDisplay),...
        'YTick', round(linspace(1,height,5)), 'YTickLabel', num2cell(round(linspace(-maxYplot,maxYplot,5))), ...
        'FontSize', 20)
    xlabel('True orientation (degree)')
    ylabel('Estimated orientation (degree)')  
end


%% Plot the estimate distribution
% Correct trials
figure
counter = 1;
binCenter = linspace(-40, 40, 41);
tempEstimateModelX = estimateModelCorrect.Xval{1};

for kk = 1 : length(estimateModelCorrect.Yval)
    estimateModelY = estimateModelCorrect.Yval{kk};
    estimateModelY = estimateModelY(:, thetaStim >=0);
    estimateSimulation = estimateSim(:, kk, thetaStim >=0);
    indexCorrectSimulation = squeeze(indexCorrectSim(:, kk, thetaStim >=0));
    for jj = 1 : size(estimateModelY, 2)            
        tempEstimateModelY = estimateModelY(:, jj); 
        tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
        tempestimateDataX = estimateSimulation(indexCorrectSimulation(:, jj) == 1, 1, jj);
        
        subplot(length(estimateModelCorrect.Yval), size(estimateModelY, 2), counter)
        hold on
        binCount = hist(tempestimateDataX, binCenter);
        binCount = binCount / (sum(binCount) * diff(binCenter(1:2)));
        bar(binCenter, binCount, 1);
 
        plot(tempEstimateModelX, tempEstimateModelY,'r','LineWidth',3)
        plot([0 0], [0 0.1], 'k--')
        xlim([-40 40])
        set(gca,'YTickLabel',[],'YTick',[])
        if kk == 1
            title(num2str(thetaStim(jj+7)))
        end
        box off
        
        counter = counter + 1;
    end
end
tightfig

% Incorrect trials
figure
counter = 1;
binCenter = linspace(-45, 45, 33);
tempEstimateModelX = estimateModelIncorrect.Xval{1};

for kk = 1 : length(estimateModelIncorrect.Yval)
    estimateModelY = estimateModelIncorrect.Yval{kk};
    estimateModelY = estimateModelY(:, thetaStim >=0);
    estimateSimulation = estimateSim(:, kk, thetaStim >=0);
    indexCorrectSimulation = squeeze(indexCorrectSim(:, kk, thetaStim >=0));
    for jj = 1 : size(estimateModelY, 2)            
        tempEstimateModelY = estimateModelY(:, jj); 
        tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
        tempestimateDataX = estimateSimulation(indexCorrectSimulation(:, jj) == 0, 1, jj);
        
        subplot(length(estimateModelIncorrect.Yval), size(estimateModelY, 2), counter)
        hold on
        binCount = hist(tempestimateDataX, binCenter);
        binCount = binCount / (sum(binCount) * diff(binCenter(1:2)));
        bar(binCenter, binCount, 1);
 
        plot(tempEstimateModelX, tempEstimateModelY,'r','LineWidth',3)
        plot([0 0], [0 0.1], 'k--')
        xlim([-47 47])
        set(gca,'YTickLabel',[],'YTick',[])
        if kk == 1
            title(num2str(thetaStim(jj+7)))
        end
        box off
        
        counter = counter + 1;
    end
end
tightfig

%% Plot the psychometric curve
pCW_Simulation = NaN(length(stdSensory), length(thetaStim));
for ii = 1 : length(stdSensory)
    for jj = 1 : length(thetaStim)
        tempDecision = squeeze(discriminateSim(:, ii, jj));
        pCW_Simulation(ii, jj) = sum(tempDecision==1) / sum(~isnan(tempDecision));
    end
end

figure
plot(thetaStim, pCW_Simulation, 'o-')
xlabel('Stimulus orientation (deg)')
ylabel('Fraction CW')