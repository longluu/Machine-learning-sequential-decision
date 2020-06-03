%%%%%%%% The sampling Bayesian model of conditioned perception %%%%%%%%%%
%%%%%% Old version: Condition the prior only %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
% 


%% Set the parameters of the model
subjectID = 'average';
flagSC = 1; % 1: self-conditioned model
            % 0: standard Bayes
flagDecisionGiven = 0;
params = [3.66      5.50      9.24           0.0000       20       5.35     0.32      3.60];
stdSensory = params(1:3);
stdMemory = params(6);
priorRange = params(5);
stdMotor = params(8);                     
smoothFactor = params(7);
thetaStim = -21:1:21;
pC = [0.5, 0.5]'; % [cw ccw]
pthcw = priorRange;
pthccw = -priorRange;
nTrialPerCondition = 1000;
fontSize = 20;
fixedConsistencyIndex = 1;
marginalizedVersion = 1;
if fixedConsistencyIndex
    consistencyIndex = 1;
end
includeIncongruentTrials = '';
rangeCollapse = round(length(thetaStim)/2);
incongruentType = 0; % 0: no incongruent
                     % 1:laspe (wrong button)
                     % 2: motor noise
                     % 3: redo decision
                     % 4: wrong button and forget stimulus
removeBoundary = 1;
stdBoundary = 5;
collapseScatterPlot = 0;                              
lapseRate = 0;

if incongruentType == 1
    lapseRate = 1;
elseif incongruentType == 2
elseif incongruentType == 3
    lapseRateForget = 1;
elseif incongruentType == 4
    lapseRate = 1;    
end

%% Modeling
dstep = 0.1;
rangeth = [-60 60];
th = rangeth(1):dstep:rangeth(2);
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

pTheta_CCW = pthGC(2,:);
pTheta_CW = pthGC(1,:);
pTheta_CCW = repmat(pTheta_CCW, nTrialPerCondition, 1);
pTheta_CW = repmat(pTheta_CW, nTrialPerCondition, 1);
pTheta_C_new = zeros(nTrialPerCondition, nth);

%%% SAMPLING version
discriminateSim = NaN(nTrialPerCondition, length(stdSensory), length(thetaStim));
estimateSim = NaN(nTrialPerCondition, length(stdSensory), length(thetaStim));
biasEstimateSim = NaN(nTrialPerCondition, length(stdSensory), length(thetaStim));

if ~flagDecisionGiven
    if incongruentType == 0
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

                if ~fixedConsistencyIndex
                    pCCW_m = pCCW_m ./ (pCCW_m + pCW_m);
                    consistencyIndex = repmat(pCCW_m, 1, nth);
                    pTheta_C_new(tempDiscriminate == -1, :) = consistencyIndex(tempDiscriminate == -1,:) .* pTheta_CCW(tempDiscriminate == -1, :) + ...
                                                           (1-consistencyIndex(tempDiscriminate == -1,:)) .* pTheta_CW(tempDiscriminate == -1, :);
                    pTheta_C_new(tempDiscriminate == 1, :) = consistencyIndex(tempDiscriminate == 1,:) .* pTheta_CCW(tempDiscriminate == 1, :) + ...
                                                            (1-consistencyIndex(tempDiscriminate == 1,:)) .* pTheta_CW(tempDiscriminate == 1, :);
                else
                    pTheta_C_new(tempDiscriminate == -1, :) = consistencyIndex * pTheta_CCW(tempDiscriminate == -1, :) + ...
                                                           (1-consistencyIndex) * pTheta_CW(tempDiscriminate == -1, :);
                    pTheta_C_new(tempDiscriminate == 1, :) = consistencyIndex * pTheta_CW(tempDiscriminate == 1, :) + ...
                                                            (1-consistencyIndex) * pTheta_CCW(tempDiscriminate == 1, :);                    
                end

                % Jitter the boundary
                if removeBoundary
                    thetaShift = normrnd(0, stdBoundary, nTrialPerCondition, 1);
                    stepShift = round(thetaShift/dstep);
                    for ii = 1 : nTrialPerCondition
                        pTheta_C_new(ii, :) = circshift(pTheta_C_new(ii, :), stepShift(ii), 2);
                    end
                end
                
                % Estimation
                Mm = normrnd(m, stdMemory);                    
                pMm_Theta = normpdf(repmat(Mm,1,nth), repmat(th,nTrialPerCondition,1), sqrt(stdSensory(jj)^2 + stdMemory^2));
                pTheta_Mm_c = pTheta_C_new .* pMm_Theta;
                pMm_c = trapz(th, pTheta_Mm_c, 2);
                pTheta_m_Norm = pTheta_Mm_c ./ repmat(pMm_c, 1, nth);
                thetaEstimate = trapz(th, repmat(th, nTrialPerCondition, 1) .* pTheta_m_Norm, 2);
                thetaResponse = normrnd(thetaEstimate, stdMotor);
                if isempty(includeIncongruentTrials)
                    indexExcludeIncongruent = sign(thetaResponse) ~= sign(thetaEstimate);
                    thetaResponse(indexExcludeIncongruent) = NaN;
                    discriminateSim(indexExcludeIncongruent, jj, kk) = NaN;
                end
                estimateSim(:, jj, kk) = thetaResponse;
                biasEstimateSim(:, jj, kk) = thetaStim(kk) - thetaResponse;
            end
        end
    elseif incongruentType == 1
        nFlipDecisionPerCondition = round(lapseRate * nTrialPerCondition);
        for jj = 1 : length(stdSensory)
            for kk = 1 : length(thetaStim)   
                % Draw a sample and measure
                m = normrnd(thetaStim(kk), stdSensory(jj),[nTrialPerCondition 1]);

                % Discrimination (Inference from p(H|m))
                pM_Theta = normpdf(repmat(th,nTrialPerCondition,1), repmat(m,1,nth), stdSensory(jj));
                pCCW_m = pC(1) .* squeeze(sum((pTheta_CCW .* pM_Theta),2)); 
                pCW_m = pC(2) .* squeeze(sum((pTheta_CW .* pM_Theta),2)); 
                discriminateSim(pCCW_m >= pCW_m, jj, kk) = -1;
                discriminateSim(pCCW_m < pCW_m, jj, kk) = 1;
                tempDiscriminate = squeeze(discriminateSim(:,jj,kk));

                % Choose nTrialPerCondition-nFlipDecisionPerCondition decisions that will not flip
                indShuffle = Shuffle(1:nTrialPerCondition);
                indNotFlip = indShuffle(1:nTrialPerCondition-nFlipDecisionPerCondition);

                if ~fixedConsistencyIndex
                    pCCW_m = pCCW_m ./ (pCCW_m + pCW_m);
                    consistencyIndex = repmat(pCCW_m, 1, nth);
                    pTheta_C_new(tempDiscriminate == -1, :) = consistencyIndex(tempDiscriminate == -1,:) .* pTheta_CCW(tempDiscriminate == -1, :) + ...
                                                           (1-consistencyIndex(tempDiscriminate == -1,:)) .* pTheta_CW(tempDiscriminate == -1, :);
                    pTheta_C_new(tempDiscriminate == 1, :) = consistencyIndex(tempDiscriminate == 1,:) .* pTheta_CCW(tempDiscriminate == 1, :) + ...
                                                            (1-consistencyIndex(tempDiscriminate == 1,:)) .* pTheta_CW(tempDiscriminate == 1, :);
                else
                    pTheta_C_new(tempDiscriminate == -1, :) = consistencyIndex * pTheta_CCW(tempDiscriminate == -1, :) + ...
                                                           (1-consistencyIndex) * pTheta_CW(tempDiscriminate == -1, :);
                    pTheta_C_new(tempDiscriminate == 1, :) = consistencyIndex * pTheta_CW(tempDiscriminate == 1, :) + ...
                                                            (1-consistencyIndex) * pTheta_CCW(tempDiscriminate == 1, :);                    
                end

                % Estimation
                Mm = normrnd(m, stdMemory);                    
                pMm_Theta = normpdf(repmat(Mm,1,nth), repmat(th,nTrialPerCondition,1), sqrt(stdSensory(jj)^2 + stdMemory^2));
                pTheta_Mm_c = pTheta_C_new .* pMm_Theta;
                pMm_c = trapz(th, pTheta_Mm_c, 2);
                pTheta_m_Norm = pTheta_Mm_c ./ repmat(pMm_c, 1, nth);
                thetaEstimate = trapz(th, repmat(th, nTrialPerCondition, 1) .* pTheta_m_Norm, 2);
                thetaResponse = normrnd(thetaEstimate, stdMotor);
                thetaResponse(sign(thetaResponse) ~= sign(thetaEstimate)) = NaN;
                estimateSim(:, jj, kk) = thetaResponse;
                biasEstimateSim(:, jj, kk) = thetaStim(kk) - thetaResponse;
            end
        end
    elseif incongruentType == 2
        for jj = 1 : length(stdSensory)
            for kk = 1 : length(thetaStim)   
                % Draw a sample and measure
                m = normrnd(thetaStim(kk), stdSensory(jj),[nTrialPerCondition 1]);

                % Discrimination (Inference from p(H|m))
                pM_Theta = normpdf(repmat(th,nTrialPerCondition,1), repmat(m,1,nth), stdSensory(jj));
                pCCW_m = pC(1) .* squeeze(sum((pTheta_CCW .* pM_Theta),2)); 
                pCW_m = pC(2) .* squeeze(sum((pTheta_CW .* pM_Theta),2)); 
                discriminateSim(pCCW_m >= pCW_m, jj, kk) = -1;
                discriminateSim(pCCW_m < pCW_m, jj, kk) = 1;
                tempDiscriminate = squeeze(discriminateSim(:,jj,kk));

                if ~fixedConsistencyIndex
                    pCCW_m = pCCW_m ./ (pCCW_m + pCW_m);
                    consistencyIndex = repmat(pCCW_m, 1, nth);
                    pTheta_C_new(tempDiscriminate == -1, :) = consistencyIndex(tempDiscriminate == -1,:) .* pTheta_CCW(tempDiscriminate == -1, :) + ...
                                                           (1-consistencyIndex(tempDiscriminate == -1,:)) .* pTheta_CW(tempDiscriminate == -1, :);
                    pTheta_C_new(tempDiscriminate == 1, :) = consistencyIndex(tempDiscriminate == 1,:) .* pTheta_CCW(tempDiscriminate == 1, :) + ...
                                                            (1-consistencyIndex(tempDiscriminate == 1,:)) .* pTheta_CW(tempDiscriminate == 1, :);
                else
                    pTheta_C_new(tempDiscriminate == -1, :) = consistencyIndex * pTheta_CCW(tempDiscriminate == -1, :) + ...
                                                           (1-consistencyIndex) * pTheta_CW(tempDiscriminate == -1, :);
                    pTheta_C_new(tempDiscriminate == 1, :) = consistencyIndex * pTheta_CW(tempDiscriminate == 1, :) + ...
                                                            (1-consistencyIndex) * pTheta_CCW(tempDiscriminate == 1, :);                    
                end

                % Estimation
                Mm = normrnd(m, stdMemory);                    
                pMm_Theta = normpdf(repmat(Mm,1,nth), repmat(th,nTrialPerCondition,1), sqrt(stdSensory(jj)^2 + stdMemory^2));
                pTheta_Mm_c = pTheta_C_new .* pMm_Theta;
                pMm_c = trapz(th, pTheta_Mm_c, 2);
                pTheta_m_Norm = pTheta_Mm_c ./ repmat(pMm_c, 1, nth);
                thetaEstimate = trapz(th, repmat(th, nTrialPerCondition, 1) .* pTheta_m_Norm, 2);
                thetaResponse = normrnd(thetaEstimate, stdMotor);
                indRemove = sign(thetaResponse) == sign(thetaEstimate);
                thetaResponse(indRemove) = NaN;
                estimateSim(:, jj, kk) = thetaResponse;
                biasEstimateSim(:, jj, kk) = thetaStim(kk) - thetaResponse;
            end
        end
    elseif incongruentType == 3
        for jj = 1 : length(stdSensory)
            for kk = 1 : length(thetaStim)   
                % Draw a sample and measure
                m = normrnd(thetaStim(kk), stdSensory(jj),[nTrialPerCondition 1]);

                % Discrimination (Inference from p(H|m))
                pM_Theta = normpdf(repmat(th,nTrialPerCondition,1), repmat(m,1,nth), stdSensory(jj));
                pCCW_m = pC(1) .* squeeze(sum((pTheta_CCW .* pM_Theta),2)); 
                pCW_m = pC(2) .* squeeze(sum((pTheta_CW .* pM_Theta),2)); 
                discriminateSim(pCCW_m >= pCW_m, jj, kk) = -1;
                discriminateSim(pCCW_m < pCW_m, jj, kk) = 1;
                tempDiscriminate = squeeze(discriminateSim(:,jj,kk));

                if ~fixedConsistencyIndex
                    pCCW_m = pCCW_m ./ (pCCW_m + pCW_m);
                    consistencyIndex = repmat(pCCW_m, 1, nth);
                    pTheta_C_new(tempDiscriminate == -1, :) = consistencyIndex(tempDiscriminate == -1,:) .* pTheta_CCW(tempDiscriminate == -1, :) + ...
                                                           (1-consistencyIndex(tempDiscriminate == -1,:)) .* pTheta_CW(tempDiscriminate == -1, :);
                    pTheta_C_new(tempDiscriminate == 1, :) = consistencyIndex(tempDiscriminate == 1,:) .* pTheta_CCW(tempDiscriminate == 1, :) + ...
                                                            (1-consistencyIndex(tempDiscriminate == 1,:)) .* pTheta_CW(tempDiscriminate == 1, :);
                else
                    pTheta_C_new(tempDiscriminate == -1, :) = consistencyIndex * pTheta_CCW(tempDiscriminate == -1, :) + ...
                                                           (1-consistencyIndex) * pTheta_CW(tempDiscriminate == -1, :);
                    pTheta_C_new(tempDiscriminate == 1, :) = consistencyIndex * pTheta_CW(tempDiscriminate == 1, :) + ...
                                                            (1-consistencyIndex) * pTheta_CCW(tempDiscriminate == 1, :);                    
                end
                if thetaStim(kk) ~= 0
                    decisionTrue = sign(thetaStim(kk));
                else
                    decisionTrue = sign(round(rand(1,nTrialPerCondition)) -0.5);
                end

                % Estimation
                Mm = normrnd(m, stdMemory);                    
                pMm_Theta = normpdf(repmat(Mm,1,nth), repmat(th,nTrialPerCondition,1), sqrt(stdSensory(jj)^2 + stdMemory^2));
                pTheta_Mm_c = pTheta_C_new .* pMm_Theta;
                pMm_c = trapz(th, pTheta_Mm_c, 2);
                pTheta_m_Norm = pTheta_Mm_c ./ repmat(pMm_c, 1, nth);
                thetaEstimate = trapz(th, repmat(th, nTrialPerCondition, 1) .* pTheta_m_Norm, 2);
%                     thetaResponse = thetaEstimate;
                thetaResponse = normrnd(thetaEstimate, stdMotor);
                thetaResponse(sign(thetaResponse)==decisionTrue') = NaN;
                estimateSim(:, jj, kk) = thetaResponse;
                biasEstimateSim(:, jj, kk) = thetaStim(kk) - thetaResponse;
            end
        end  
    elseif incongruentType == 4
        nFlipDecisionPerCondition = round(lapseRate * nTrialPerCondition);
        for jj = 1 : length(stdSensory)
            for kk = 1 : length(thetaStim)   
                % Draw a sample and measure
                m = normrnd(thetaStim(kk), stdSensory(jj),[nTrialPerCondition 1]);

                % Discrimination (Inference from p(H|m))
                pM_Theta = normpdf(repmat(th,nTrialPerCondition,1), repmat(m,1,nth), stdSensory(jj));
                pCCW_m = pC(1) .* squeeze(sum((pTheta_CCW .* pM_Theta),2)); 
                pCW_m = pC(2) .* squeeze(sum((pTheta_CW .* pM_Theta),2)); 
                discriminateSim(pCCW_m >= pCW_m, jj, kk) = -1;
                discriminateSim(pCCW_m < pCW_m, jj, kk) = 1;
                tempDiscriminate = squeeze(discriminateSim(:,jj,kk));

                % Choose nTrialPerCondition-nFlipDecisionPerCondition decisions that will not flip
                indShuffle = Shuffle(1:nTrialPerCondition);
                indNotFlip = indShuffle(1:nTrialPerCondition-nFlipDecisionPerCondition);

                if ~fixedConsistencyIndex
                    pCCW_m = pCCW_m ./ (pCCW_m + pCW_m);
                    consistencyIndex = repmat(pCCW_m, 1, nth);
                    pTheta_C_new(tempDiscriminate == -1, :) = consistencyIndex(tempDiscriminate == -1,:) .* pTheta_CCW(tempDiscriminate == -1, :) + ...
                                                           (1-consistencyIndex(tempDiscriminate == -1,:)) .* pTheta_CW(tempDiscriminate == -1, :);
                    pTheta_C_new(tempDiscriminate == 1, :) = consistencyIndex(tempDiscriminate == 1,:) .* pTheta_CCW(tempDiscriminate == 1, :) + ...
                                                            (1-consistencyIndex(tempDiscriminate == 1,:)) .* pTheta_CW(tempDiscriminate == 1, :);
                else
                    pTheta_C_new(tempDiscriminate == -1, :) = consistencyIndex * pTheta_CCW(tempDiscriminate == -1, :) + ...
                                                           (1-consistencyIndex) * pTheta_CW(tempDiscriminate == -1, :);
                    pTheta_C_new(tempDiscriminate == 1, :) = consistencyIndex * pTheta_CW(tempDiscriminate == 1, :) + ...
                                                            (1-consistencyIndex) * pTheta_CCW(tempDiscriminate == 1, :);                    
                end

                % Estimation
                thetaEstimate = tempDiscriminate * priorRange/2;
                thetaEstimate(indNotFlip) = NaN;
                thetaResponse = normrnd(thetaEstimate, stdMotor);
                estimateSim(:, jj, kk) = thetaResponse;
                biasEstimateSim(:, jj, kk) = thetaStim(kk) - thetaResponse;
            end
        end          
    end
else
    for jj = 1 : length(stdSensory)
        for kk = 1 : length(thetaStim)   
            % Discrimination (Inference from p(H|m))
            if incongruentType == 1
                discriminateSim(:, jj, kk) = -sign(thetaStim(kk));
                if thetaStim(kk) < 0
                    pTheta_C_new = pTheta_CW;
                elseif thetaStim(kk) > 0
                    pTheta_C_new = pTheta_CCW;
                else
                    inCCW = round(rand(1,nTrialPerCondition));
                    pTheta_C_new(inCCW==1, :) = squeeze(pTheta_CW(inCCW==1,:));
                    pTheta_C_new(inCCW==0, :) = squeeze(pTheta_CCW(inCCW==0,:));
                    discriminateSim(inCCW==1, jj, kk) = -1;
                    discriminateSim(inCCW==0, jj, kk) = 1;
                end                
            else
                discriminateSim(:, jj, kk) = sign(thetaStim(kk));
                if thetaStim(kk) < 0
                    pTheta_C_new = pTheta_CCW;
                elseif thetaStim(kk) > 0
                    pTheta_C_new = pTheta_CW;
                else
                    inCCW = round(rand(1,nTrialPerCondition));
                    pTheta_C_new(inCCW==1, :) = squeeze(pTheta_CCW(inCCW==1,:));
                    pTheta_C_new(inCCW==0, :) = squeeze(pTheta_CW(inCCW==0,:));
                    discriminateSim(inCCW==1, jj, kk) = -1;
                    discriminateSim(inCCW==0, jj, kk) = 1;
                end
            end
            
            % Estimation            
            m = normrnd(thetaStim(kk), stdSensory(jj),[nTrialPerCondition 1]);
            Mm = normrnd(m, stdMemory);                    
            pMm_Theta = normpdf(repmat(Mm,1,nth), repmat(th,nTrialPerCondition,1), sqrt(stdSensory(jj)^2 + stdMemory^2));
            pTheta_Mm_c = pTheta_C_new .* pMm_Theta;
            pMm_c = trapz(th, pTheta_Mm_c, 2);
            pTheta_m_Norm = pTheta_Mm_c ./ repmat(pMm_c, 1, nth);
            thetaEstimate = trapz(th, repmat(th, nTrialPerCondition, 1) .* pTheta_m_Norm, 2);
            thetaResponse = normrnd(thetaEstimate, stdMotor);
            if incongruentType == 2 
                indexExcludeCongruent = sign(thetaResponse) == sign(thetaEstimate);
                thetaResponse(indexExcludeCongruent) = NaN;
                discriminateSim(indexExcludeCongruent, jj, kk) = NaN;                
            elseif incongruentType == 0 || incongruentType == 1
                indexExcludeIncongruent = sign(thetaResponse) ~= sign(thetaEstimate);
                thetaResponse(indexExcludeIncongruent) = NaN;
                discriminateSim(indexExcludeIncongruent, jj, kk) = NaN;
            end
            estimateSim(:, jj, kk) = thetaResponse;
            biasEstimateSim(:, jj, kk) = thetaStim(kk) - thetaResponse;
        end
    end
end
a = estimateSim(:);
sum(~isnan(a)) / length(a)

%%% MARGINALIZED version
if marginalizedVersion == 1
    estimateModel.Xval = cell(length(stdSensory));
    estimateModel.Yval = cell(length(stdSensory));  
    biasEstimateTheoryCCW = NaN(length(stdSensory), round(length(thetaStim)/2));
    biasEstimateTheoryCW = NaN(length(stdSensory), round(length(thetaStim)/2));
    estimateTheoryCCW = NaN(length(stdSensory), length(thetaStim));
    estimateTheoryCW = NaN(length(stdSensory), length(thetaStim));
    

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

        M = repmat(m',1,nth);
        MM_m = repmat(mm',1,nm);
        MM_th = repmat(mm',1,nth); 
        THm = repmat(th, nm, 1); 
        THmm = repmat(th, nmm, 1); 

        %% Generative (forward)
        % orientation noise
        pmGth = exp(-((M-THm).^2)./(2*stdSensory(kk)^2));
        pmGth = pmGth./(repmat(sum(pmGth,1),nm,1)); 

        %% Inference
        % 1: categorical judgment
        if ~flagDecisionGiven
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
        else
            PChGtheta_lapse = NaN(2, length(thetaStim));
            PChGtheta_lapse(1, thetaStim > 0) = 1;
            PChGtheta_lapse(1, thetaStim < 0) = 0;
            PChGtheta_lapse(1, thetaStim == 0) = 0.5;
            PChGtheta_lapse(2, :) = 1 - PChGtheta_lapse(1, :);
        end
        
        % 2: estimation
        if ~flagDecisionGiven
            % memory noise
            pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemory^2)); 
            pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1));   

            % attention marginalization: compute distribution only over those ms that lead to cw decision!
            pmmGthChcw = pmmGm * (pmGth .* repmat(PChGm(1,:)',1,nth));
            pmmGthChccw = pmmGm * (pmGth .* repmat(PChGm(2,:)',1,nth));
            
            pthGmmChcw = (pmmGthChcw.*repmat(pthGC(1,:),nmm,1))';
            pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
            pthGmmChcw(isnan(pthGmmChcw)) = 0;      
            
            pthGmmChccw = (pmmGthChccw.*repmat(pthGC(2,:),nmm,1))';
            pthGmmChccw = pthGmmChccw./repmat(sum(pthGmmChccw,1),nth,1);
            pthGmmChccw(isnan(pthGmmChccw)) = 0;                        
        else
            pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); % p(mm|th) = N(th, sm^2 + smm^2)
            pmmGth = pmmGth./(repmat(sum(pmmGth,1),nmm,1)); 
            
            pthGmmChcw = (pmmGth.*repmat(pthGC(1,:),nmm,1))';
            pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
            pthGmmChcw(isnan(pthGmmChcw)) = 0;
            
            pthGmmChccw = (pmmGth.*repmat(pthGC(2,:),nmm,1))';
            pthGmmChccw = pthGmmChccw./repmat(sum(pthGmmChccw,1),nth,1);
            pthGmmChccw(isnan(pthGmmChccw)) = 0;            
        end
        

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
        if ~flagDecisionGiven
            % attention marginalization: compute distribution only over those ms that lead to cw decision!
            pmmGthChcw = pmmGthChcw(:, ismember(th, thetaStim));   
            b = repmat(a',1,length(thetaStim)) .* pmmGthChcw(indKeepCw, :);        
        else
            b = repmat(a',1,length(thetaStim)) .* pmmGth(indKeepCw, ismember(th, thetaStim));   
        end
        pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChcw(pthhGthChcw < 0) = 0; 

        a = 1./gradient(EthChccw,dstep);
        if ~flagDecisionGiven
            % attention marginalization: compute distribution only over those ms that lead to cw decision!
            pmmGthChccw = pmmGthChccw(:, ismember(th, thetaStim));        
            b = repmat(a',1,length(thetaStim)) .* pmmGthChccw(indKeepCcw, :);        
        else
            b = repmat(a',1,length(thetaStim)) .* pmmGth(indKeepCcw, ismember(th, thetaStim));
        end  
        pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
        % add motor noise
        pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
        pthhGthChccw(pthhGthChccw < 0) = 0; 

        % remove incongruent trial
        if isempty(includeIncongruentTrials)
            pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
            pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            
            
            % modify psychometric curve p(Chat|theta, Congruent) ~ p(Congruent| Chat, theta) * p(Chat|Theta)
            pCongruentGcwTh = sum(pthhGthChcw(th' >= 0, :));
            pCongruentGccwTh = sum(pthhGthChccw(th' <= 0, :));
            PChGtheta_lapse = PChGtheta_lapse .* [pCongruentGcwTh; pCongruentGccwTh];
            PChGtheta_lapse = PChGtheta_lapse ./ repmat(sum(PChGtheta_lapse, 1), 2, 1);
            
            % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
            pthhGthChccw(th' >= 0, :) = 0;
            pthhGthChcw(th' < 0, :) = 0;
        end 
        pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
        pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            
        mthhGthChcw = th * pthhGthChcw;
        mthhGthChccw = th * pthhGthChccw;

        
        pthhGth = pthhGthChcw.*repmat(PChGtheta_lapse(1,:),nth,1) + pthhGthChccw.*repmat(PChGtheta_lapse(2,:),nth,1);

        estimateModel.Xval{kk} = th;
        estimateModel.Yval{kk} = pthhGth;
        estimateTheoryCCW(kk,:) = mthhGthChccw;
        estimateTheoryCW(kk,:) = mthhGthChcw;
    end                
end

% Analyze and plot
estimateCollapse = cell(length(stdSensory), rangeCollapse);
estimateResponseSimCCW = NaN(size(estimateSim));
estimateResponseSimCCW(discriminateSim == -1) = estimateSim(discriminateSim == -1);
estimateResponseSimCCWAve = squeeze(nanmean(estimateResponseSimCCW, 1));
estimateResponseSimCCWStd = squeeze(nanstd(estimateResponseSimCCW, 1))./(sqrt(squeeze(nansum(~isnan(estimateResponseSimCCW),1))));
estimateResponseSimCW = NaN(size(estimateSim));
estimateResponseSimCW(discriminateSim == 1) = estimateSim(discriminateSim == 1);
estimateResponseSimCWAve = squeeze(nanmean(estimateResponseSimCW, 1));
estimateResponseSimCWStd = squeeze(nanstd(estimateResponseSimCW, 1))./(sqrt(squeeze(nansum(~isnan(estimateResponseSimCW),1))));

hAverage = figure;
lineWidth = 4;
figPos = [0.3, 0.2, 0.5, 0.7];
set(hAverage,'Units','normalized','Position',figPos)
hold on
legendName = cell(1,length(stdSensory));
colorName = {'r', 'g', 'b', 'cyan', 'magenta', 'y'};
maxPlot = max(thetaStim);
indexColor = 1;
hLegend = NaN(1, length(stdSensory));
biasMean = NaN(length(stdSensory), rangeCollapse);
biasStd = NaN(length(stdSensory), rangeCollapse);
biasTheory = NaN(length(stdSensory), rangeCollapse);

for ii = 1 : length(stdSensory)
    set(gca,'FontSize',fontSize)
    hShade = shadedErrorBar(thetaStim, estimateResponseSimCCWAve(ii,:), estimateResponseSimCCWStd(ii,:),... 
                         {'Color', colorName{ii}, 'LineWidth', lineWidth},1,0,0);        
    hLegend(ii) = hShade.patch;
    hold on
    shadedErrorBar(thetaStim, estimateResponseSimCWAve(ii, :), estimateResponseSimCWStd(ii, :),... 
                         {'Color', colorName{ii}, 'LineWidth', lineWidth},1,0,0);        
    if  marginalizedVersion
        hLegend(ii) = plot(thetaStim, estimateTheoryCCW(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);
        plot(thetaStim, estimateTheoryCW(ii,:), 'Color', colorName{ii}, 'LineWidth', lineWidth);
    end
    legendName{ii} = ['Noise level = ' num2str(stdSensory(ii))];
    xlabel('True angle (degree)')
    ylabel('Angle estimate (degree)')  
%         title(['Subject ' upper(subjectID) ]);
    plot([-maxPlot maxPlot], [-maxPlot maxPlot], 'k--', 'LineWidth', 2)
    xlim([-maxPlot maxPlot])
    ylim([-maxPlot maxPlot])
    text(6, -15, '+/-1SEM', 'FontSize', 15)

    % Get the bias
    tempEstimateCW = squeeze(estimateResponseSimCW(:,ii,:));
    tempEstimateCCW = squeeze(estimateResponseSimCCW(:,ii,:));
    tempBias1 = tempEstimateCW(:,rangeCollapse:end) - repmat(thetaStim(rangeCollapse:end),size(tempEstimateCW,1),1);
    tempBias2 = tempEstimateCCW(:,1:rangeCollapse) - repmat(thetaStim(1:rangeCollapse),size(tempEstimateCCW,1),1);
    tempBias = [tempBias1;-tempBias2(:,sort(1:size(tempBias2,2),'descend'))];
    biasMean(ii,:) = squeeze(nanmean(tempBias,1));
    biasStd(ii,:) = squeeze(nanstd(tempBias,1))./sqrt(nansum(~isnan(tempBias),1)); 

    if marginalizedVersion
        tempBiasTheory1 = estimateTheoryCW(ii,rangeCollapse:end) - thetaStim(rangeCollapse:end);
        tempBiasTheory2 = estimateTheoryCCW(ii,1:rangeCollapse) - thetaStim(1:rangeCollapse);
        tempBiasTheory = [tempBiasTheory1;-tempBiasTheory2(:,sort(1:size(tempBiasTheory2,2),'descend'))];  
        biasTheory(ii,:) = squeeze(nanmean(tempBiasTheory,1));
    end

    tempEstimateCollapse = [-fliplr(tempEstimateCCW); tempEstimateCW];
    for jj = 1 : rangeCollapse
        estimateCollapse{ii,jj} = tempEstimateCollapse(:,jj+rangeCollapse-1);
    end
end
legend(hLegend, legendName, 'Location', 'NorthWest')

hBiasTheory = figure;
lineWidth = 4;
legendName = cell(1,length(stdSensory));
colorName = {'b', 'g', 'r', 'cyan', 'magenta', 'y'};
maxPlot = max(thetaStim);
indexColor = 1;
hLegend = NaN(1, length(stdSensory));
rangeCollapse = round(length(thetaStim)/2);
set(gca,'FontSize',fontSize)
for ii = 1 : length(stdSensory)
    hold on
    hTheory = plot(thetaStim(rangeCollapse:end), biasTheory(ii,:), ... 
                         'Color', colorName{ii}, 'LineWidth', lineWidth);        
    hLegend(ii) = hTheory;
        shadedErrorBar(thetaStim(rangeCollapse:end), biasMean(ii,:), biasStd(ii,:),... 
                             {'Color', colorName{ii}, 'LineWidth', lineWidth},1,0,0);        
    legendName{ii} = ['Noise = ' num2str(stdSensory(ii))];
    xlabel('True angle (degree)')
    ylabel('Bias (degree)')  
%     title(['Subject ' upper(subjectID) ]);
    xlim([0 maxPlot])
    ylim([-10 15])
    box on
    text(6, -15, '+/-1SEM', 'FontSize', 15)
end
plot([0 maxPlot], [0 0], 'k--', 'LineWidth', lineWidth)

%% Plot the smoothed raw data
if length(thetaStim) == 15
    hScatter = figure;
    figPos = [0.1, 0.2, 0.8, 0.6];
    set(hScatter,'Units','normalized','Position',figPos)
    hold on
    maxYplot = 40;
    imageData = zeros(length(-maxYplot:maxYplot), length(-22:22), length(stdSensory));
    for ii = 1 : length(stdSensory)
        for jj = 1 :length(thetaStim)
            tempthetaStimEst = squeeze(estimateSim(:,ii,jj));
            tempthetaStimEst(isnan(tempthetaStimEst)) = [];
            for kk = 1 : length(tempthetaStimEst)
                imageData(maxYplot-round(tempthetaStimEst(kk))+1, round(thetaStim(jj))+23, ii) = ...
                    imageData(maxYplot-round(tempthetaStimEst(kk))+1, round(thetaStim(jj))+23, ii) + 1;
            end
        end
        for jj = 1 :length(thetaStim)
            imageData(:, 3*jj-2, ii) = imageData(:, 3*jj-1, ii);
            imageData(:, 3*jj, ii) = imageData(:, 3*jj-1, ii);
        end
        myfilter = fspecial('gaussian', [2 1], 0.5);
        smoothImage = imfilter(imageData(:,:,ii), myfilter, 'replicate');
        imageData(:,:,ii) = round(smoothImage*255/max(smoothImage(:)));

        subplot(1,length(stdSensory),ii)
        hold on
        tempImage = uint8(imageData(:,:,ii));
        [height, width] = size(tempImage);
        widthPlot = round(width);
        tempImage = imresize(tempImage, [height widthPlot], 'bicubic');
        imshow(tempImage)
        axis on
        hold on
        set(gca,'FontSize',fontSize)
        set(gca, 'ylim', [1 height], 'xlim', [1 widthPlot], ...
            'XTick', round(linspace(1,widthPlot,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
            'YTick', round(linspace(1,height,5)), 'YTickLabel', num2cell(round(linspace(maxYplot,-maxYplot,5))))
        plot([1 widthPlot], [round(length(-maxYplot:maxYplot)/2)+22 round(length(-maxYplot:maxYplot)/2)-22], '--b', 'LineWidth', 1.1)
        title(['Stimulus noise = ' num2str(stdSensory(ii))])
        xlabel('True angle (degree)')
        ylabel('Estimated angle (degree)')
    end
else
    if collapseScatterPlot
        hScatter = figure;
        figPos = [0.1, 0.2, 0.8, 0.6];
        set(hScatter,'Units','normalized','Position',figPos)
        hold on
        maxYplot = 35;
        yAxis = -maxYplot:0.1:maxYplot;
        maxthetaStim = max(thetaStim);
        imageData = zeros(length(yAxis), length(thetaStim)+2);
        for ii = 1 : length(stdSensory)
            for jj = 1 :length(thetaStim)
                tempthetaStimEst = squeeze(estimateSim(:,ii,jj));
                tempthetaStimEst(isnan(tempthetaStimEst)) = [];
                tempthetaStimEst(abs(tempthetaStimEst)>maxYplot) = [];
                for kk = 1 : length(tempthetaStimEst)
                    imageData(10*round(tempthetaStimEst(kk)+maxYplot, 1), jj+1) = ...
                        imageData(10*round(tempthetaStimEst(kk)+maxYplot, 1), jj+1) + 1;
                end
            end
        end
        myfilter = fspecial('gaussian', [20 1], 2);
        smoothImage = imfilter(imageData, myfilter, 'replicate');
        imageData = round(smoothImage*255/max(smoothImage(:)));

        hold on
        tempImage = uint8(imageData);
        [height, width] = size(tempImage);
        widthPlot = round(width);
        imagesc(tempImage)
        hold on;
        axis xy;
        colormap('gray');
        plot([1 widthPlot],[round(height/2) round(height/2)],'w:', 'LineWidth', 1.5);
        plot([round(widthPlot/2) round(widthPlot/2)],[1 height],'w:', 'LineWidth', 1.5);
        plot([1 widthPlot],[1 height],'b:', 'LineWidth', 2);
        set(gca, 'ylim', [1 height], 'xlim', [1 widthPlot], ...
            'XTick', round(linspace(1,widthPlot,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
            'YTick', round(linspace(1,height,5)), 'YTickLabel', num2cell(round(linspace(-maxYplot,maxYplot,5))), ...
            'FontSize', 20)
        xlabel('True orientation (degree)')
        ylabel('Estimated orientation (degree)')        
        
        
    else
        hScatter = figure;
        figPos = [0.1, 0.2, 0.8, 0.6];
        set(hScatter,'Units','normalized','Position',figPos)
        hold on
        maxYplot = 55;
        maxthetaStim = max(thetaStim);
        imageData = zeros(length(-maxYplot:maxYplot), length(thetaStim)+2, length(stdSensory));
        for ii = 1 : length(stdSensory)
            for jj = 1 :length(thetaStim)
                tempthetaStimEst = squeeze(estimateSim(:,ii,jj));
                tempthetaStimEst(isnan(tempthetaStimEst)) = [];
                for kk = 1 : length(tempthetaStimEst)
                    imageData(maxYplot-round(tempthetaStimEst(kk))+1, jj+1, ii) = ...
                        imageData(maxYplot-round(tempthetaStimEst(kk))+1, jj+1, ii) + 1;
                end
            end
            myfilter = fspecial('gaussian', [1 1], 1);
            smoothImage = imfilter(imageData(:,:,ii), myfilter, 'replicate');
            imageData(:,:,ii) = round(smoothImage*255/max(smoothImage(:)));

            subplot(1,length(stdSensory),ii)
            hold on
            tempImage = uint8(imageData(:,:,ii));
            [height, width] = size(tempImage);
            widthPlot = round(width);
            tempImage = imresize(tempImage, [height widthPlot], 'bicubic');
            imshow(tempImage)
            axis on
            hold on
            set(gca,'FontSize',fontSize)
            set(gca, 'ylim', [1 height], 'xlim', [1 widthPlot], ...
                'XTick', round(linspace(1,widthPlot,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
                'YTick', round(linspace(1,height,5)), 'YTickLabel', num2cell(round(linspace(maxYplot,-maxYplot,5))))
            plot([1 widthPlot], [round(length(-maxYplot:maxYplot)/2)+maxthetaStim+1 round(length(-maxYplot:maxYplot)/2)-maxthetaStim-2], '--b', 'LineWidth', 1.1)
            title(['Stimulus noise = ' num2str(stdSensory(ii))])
            xlabel('True angle (degree)')
            ylabel('Estimated angle (degree)')
        end    
    end
end

%% Plot the estimate distribution
figure
counter = 1;
binCenter = linspace(-40, 40, 35);
tempEstimateModelX = estimateModel.Xval{1};

for kk = 1 : length(estimateModel.Yval)
    estimateModelY = estimateModel.Yval{kk};
    estimateModelY = estimateModelY(:, thetaStim >=0);
    estimateSimulation = estimateSim(:, kk, thetaStim >=0);
    for jj = 1 : size(estimateModelY, 2)            
        tempEstimateModelY = estimateModelY(:, jj); 
        tempEstimateModelY = tempEstimateModelY ./ trapz(tempEstimateModelX, tempEstimateModelY);
        tempestimateDataX = estimateSimulation(:, 1, jj);
        
        subplot(length(estimateModel.Yval), size(estimateModelY, 2), counter)
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
