% Density plot for combined incongruent error
% ********** Old version - condition the prior only ********** 
% astocker - lluu
% 06.2017
startTime = mglGetSecs;
experiment = 1; 
lapseType = 1;
dstep = 0.1;
paramsAllSubject = [2.5611    4.8570    7.1153      0.0000    28.6016      0.8805    0.8461    2.6327;
                    4.0300    5.3505    8.1590      0.0000    41.7649   2.0957    0.0257    1.8431;
                    3.8692    6.4042   12.3882      0.0000    39.5010  14.0727    0.2322    2.9146;
                    3.7041    4.6566   10.9916      0.0000    36.4739   0.2433    0.6481    5.4045;
                    3.9226    5.7843   10.9330      0.0000    39.0534  14.2463    0.1492    5.8917];

% paramsAllSubject = [2.4118    4.6421    6.6797           0          22.8689     0.7788    0.8160    2.6327;
%                     7.5309    9.5173   15.1030           0          21.3314    0.5988    0.8698    2.9350;
%                     4.4382    5.9486    7.8313           0          15.3890    1.3829    0.4445    2.9916;
%                     5.3096    5.9372    7.0886           0          16.6193    1.3993    0.1479    2.5348;
%                     4.3799    4.4121    8.9346           0          33.0360   14.8571    0.6678    2.0456];

                
incongruentErrorWeight = [0.1111 0.7778 5.7778 4.6111 3.8333];
incongruentErrorWeight = incongruentErrorWeight / sum(incongruentErrorWeight);
lapseAll = [0 0.76 4.75 1.88 2.86]/100;
incongruentTotal = NaN(size(lapseAll));

% stimulus orientation
thetaStim = -21:0.1:21; % 
thetaStim = round(thetaStim, -log10(dstep));

rangeth = [-42 42];
th = rangeth(1):dstep:rangeth(2);
th = round(th, -log10(dstep));
nth = length(th);

pthhGthAverage = zeros(nth, length(thetaStim));

if experiment == 1 || experiment == 2
    for ii = 1 : size(paramsAllSubject, 1)
        paramsAll = paramsAllSubject(ii, :);
        lapseRate = paramsAll(4);

        % sensory noise
        stdSensory = paramsAll(1:3);

        % memory recall noise
        stdMemory = paramsAll(6);

        % motor noise;
        stdMotor = paramsAll(8);

        % priors
        pC = [0.5, 0.5]'; % [cw ccw]
        pthcw = paramsAll(5);
        pthccw = -paramsAll(5);
        smoothFactor = paramsAll(7);

        pthGC = zeros(2,nth);

        pthGC(1,:) = TukeyWindow([0 pthcw], 0, smoothFactor, th);
        pthGC(2,:) = TukeyWindow([pthccw 0], 1, smoothFactor, th);

        % LOOP - noise levels
        pthhGth_Collapse = zeros(nth, length(thetaStim));
        pCongruent = zeros(length(stdSensory) , length(thetaStim));
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
            
            %% Incongruent trials by motor noise
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
            pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not                
            pthhGthChcw = pthhGthChcw .* repmat(PChGtheta_lapse(1,:),nth,1);
            
            a = 1./gradient(EthChccw,dstep);
            % attention marginalization: compute distribution only over those ms that lead to cw decision!
            pmmGthChccw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim)));        
            b = repmat(a',1,length(thetaStim)) .* pmmGthChccw(indKeepCcw, :);        
            pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
            % add motor noise
            pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
            pthhGthChccw(pthhGthChccw < 0) = 0; 
            pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);                 
            pthhGthChccw = pthhGthChccw .* repmat(PChGtheta_lapse(2,:),nth,1);
            
            % modify the estimate distribution p(thetaHat|theta, Chat, Incongruent)
            pthhGthChccw_motor = pthhGthChccw;
            pthhGthChcw_motor = pthhGthChcw;
            pthhGthChccw_motor(th'<= 0, :) = 0;
            pthhGthChcw_motor(th'> 0, :) = 0;
            
            %% Incongruent trials by lapse 
            if lapseType == 1
                % Wrong button, right estimate
                pthhGthChccw_lapse = pthhGthChccw;
                pthhGthChcw_lapse = pthhGthChcw;
                pthhGthChccw_lapse(th'>= 0, :) = 0;
                pthhGthChcw_lapse(th'< 0, :) = 0;                
            elseif lapseType == 2
                % Wrong button, forget stimulus
                pthhGthChcw_lapse = repmat(normpdf(th', pthcw/2, stdMotor), 1, length(thetaStim));
                pthhGthChcw_lapse = pthhGthChcw_lapse./repmat(sum(pthhGthChcw_lapse,1),nth,1);   
                pthhGthChcw_lapse = pthhGthChcw_lapse  .* repmat(PChGtheta_lapse(1,:),nth,1);
                
                pthhGthChccw_lapse = repmat(normpdf(th', pthccw/2, stdMotor), 1, length(thetaStim)) .* repmat(PChGtheta_lapse(2,:),nth,1); 
                pthhGthChccw_lapse = pthhGthChccw_lapse./repmat(sum(pthhGthChccw_lapse,1),nth,1); 
                pthhGthChccw_lapse =  pthhGthChccw_lapse .* repmat(PChGtheta_lapse(2,:),nth,1); 
            end
            
            %% Combine the two distributions
            pthhGthChcw = lapseAll(ii) * pthhGthChcw_lapse + (1-lapseAll(ii)) * pthhGthChcw_motor;            
            pthhGthChccw = lapseAll(ii) * pthhGthChccw_lapse + (1-lapseAll(ii)) * pthhGthChccw_motor;
            
            pthhGth = pthhGthChcw+ pthhGthChccw;
            pthhGth_Collapse = pthhGth_Collapse + pthhGth;
            pCongruent(kk, :) = sum(pthhGth, 1);
        end
        pthhGthAverage = pthhGthAverage + incongruentErrorWeight(ii) * pthhGth_Collapse; 
        incongruentTotal(ii) = mean(pCongruent(:));
        
        % Plot
        figure
        xRange = [-22 22];
        indX = find(thetaStim >= xRange(1) & thetaStim <= xRange(2));
        xMax = length(indX);
        yRange = [-35 35];
        indY = find(th >= yRange(1) & th <= yRange(2));
        yMax = length(indY);
        thNew = th(indY);
        indYStart = find(thNew == xRange(1));
        indYEnd = find(thNew == xRange(2));

        imagesc(pthhGth_Collapse(indY, indX));
        hold on;
        axis xy;
        colormap('gray');
        plot([1 xMax],[round(yMax/2) round(yMax/2)],'w:', 'LineWidth', 1.5);
        plot([round(xMax/2) round(xMax/2)],[1 yMax],'w:', 'LineWidth', 1.5);
        plot([1 xMax],[indYStart indYEnd],'b:', 'LineWidth', 2);
        set(gca, 'ylim', [1 yMax], 'xlim', [1 xMax], ...
            'XTick', round(linspace(1,xMax,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
            'YTick', round(linspace(1,yMax,5)), 'YTickLabel', num2cell(round(linspace(yRange(1),yRange(2),5))), ...
            'FontSize', 20)
        xlabel('True orientation (degree)')
        ylabel('Estimated orientation (degree)')        
        
    end    
elseif experiment == 3
    for ii = 1 : size(paramsAllSubject, 1)
        paramsAll = paramsAllSubject(ii, :);
        lapseRate = paramsAll(4);

        % sensory noise
        stdSensory = paramsAll(1:3);

        % memory recall noise
        stdMemory = paramsAll(6);

        % motor noise;
        stdMotor = paramsAll(8);

        % priors
        pC = [0.5, 0.5]'; % [cw ccw]
        pthcw = paramsAll(5);
        pthccw = -paramsAll(5);
        smoothFactor = paramsAll(7);

        pthGC = zeros(2,nth);

        pthGC(1,:) = TukeyWindow([0 pthcw], 0, smoothFactor, th);
        pthGC(2,:) = TukeyWindow([pthccw 0], 1, smoothFactor, th);

        % LOOP - noise levels
        pthhGth_Collapse = zeros(nth, length(thetaStim));
        pCongruent = zeros(length(stdSensory) , length(thetaStim));
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
            PChGtheta_lapse = NaN(2, length(thetaStim));
            PChGtheta_lapse(1, thetaStim > 0) = 1;
            PChGtheta_lapse(1, thetaStim < 0) = 0;
            PChGtheta_lapse(1, thetaStim == 0) = 0.5;
            PChGtheta_lapse(2, :) = 1 - PChGtheta_lapse(1, :);
            
            %% Incongruent trials by motor noise
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
            b = repmat(a',1,length(thetaStim)) .* pmmGth(indKeepCw, ismember(th, thetaStim));   
            pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
            pthhGthChcw(pthhGthChcw < 0) = 0; 
            pthhGthChcw(:, thetaStim < 0) = 0;
            
            a = 1./gradient(EthChccw,dstep);
            b = repmat(a',1,length(thetaStim)) .* pmmGth(indKeepCcw, ismember(th, thetaStim));
            pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
            pthhGthChccw(pthhGthChccw < 0) = 0; 
            pthhGthChccw(:, thetaStim > 0) = 0;
            
            %% Incongruent trials by lapse (misrember given decision)
            pthGC_new = flipud(pthGC);
            pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); % p(mm|th) = N(th, sm^2 + smm^2)
            pmmGth = pmmGth./(repmat(sum(pmmGth,1),nmm,1)); 
            pthGmmChcw = (pmmGth.*repmat(pthGC_new(1,:),nmm,1))';
            pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
            pthGmmChcw(isnan(pthGmmChcw)) = 0;

            pthGmmChccw = (pmmGth.*repmat(pthGC_new(2,:),nmm,1))';
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
            b = repmat(a',1,length(thetaStim)) .* pmmGth(indKeepCw, ismember(th, thetaStim));   
            pthhGthChcw_lapse = interp1(EthChcw,b,th,'linear','extrap');
            pthhGthChcw_lapse(pthhGthChcw_lapse < 0) = 0; 
            
            a = 1./gradient(EthChccw,dstep);
            b = repmat(a',1,length(thetaStim)) .* pmmGth(indKeepCcw, ismember(th, thetaStim));
            pthhGthChccw_lapse = interp1(EthChccw,b,th,'linear','extrap');
            pthhGthChccw_lapse(pthhGthChccw_lapse < 0) = 0; 

            pthhGthChcw_lapse(:, thetaStim < 0) = 0;
            pthhGthChccw_lapse(:, thetaStim > 0) = 0;
            
            %% Combine the two distributions
            pthhGthChcw = lapseAll(ii) * pthhGthChcw_lapse + (1-lapseAll(ii)) * pthhGthChcw;
            % add motor noise
            pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
            pthhGthChcw(pthhGthChcw < 0) = 0; 
            pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not  
            pthhGthChcw(:, thetaStim == 0) = pthhGthChcw(:, thetaStim == 0) / 2;            
            pthhGthChcw(isnan(pthhGthChcw)) = 0;
            pthhGthChcw(th > 0, :) = 0;
            
            pthhGthChccw = lapseAll(ii) * pthhGthChccw_lapse + (1-lapseAll(ii)) * pthhGthChccw;
            % add motor noise
            pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
            pthhGthChccw(pthhGthChccw < 0) = 0; 
            pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);
            pthhGthChccw(:, thetaStim == 0) = pthhGthChccw(:, thetaStim == 0) / 2;
            pthhGthChccw(isnan(pthhGthChccw)) = 0;            
            pthhGthChccw(th < 0, :) = 0;
            
            pthhGth = pthhGthChcw+ pthhGthChccw;
            pthhGth_Collapse = pthhGth_Collapse + pthhGth;
            pCongruent(kk, :) = sum(pthhGth, 1);
        end
        pthhGthAverage = pthhGthAverage + incongruentErrorWeight(ii) * pthhGth_Collapse; 
        incongruentTotal(ii) = mean(pCongruent(:));
        
        % Plot
        figure
        xRange = [-22 22];
        indX = find(thetaStim >= xRange(1) & thetaStim <= xRange(2));
        xMax = length(indX);
        yRange = [-35 35];
        indY = find(th >= yRange(1) & th <= yRange(2));
        yMax = length(indY);
        thNew = th(indY);
        indYStart = find(thNew == xRange(1));
        indYEnd = find(thNew == xRange(2));

        imagesc(pthhGth_Collapse(indY, indX));
        hold on;
        axis xy;
        colormap('gray');
        plot([1 xMax],[round(yMax/2) round(yMax/2)],'w:', 'LineWidth', 1.5);
        plot([round(xMax/2) round(xMax/2)],[1 yMax],'w:', 'LineWidth', 1.5);
        plot([1 xMax],[indYStart indYEnd],'b:', 'LineWidth', 2);
        set(gca, 'ylim', [1 yMax], 'xlim', [1 xMax], ...
            'XTick', round(linspace(1,xMax,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
            'YTick', round(linspace(1,yMax,5)), 'YTickLabel', num2cell(round(linspace(yRange(1),yRange(2),5))), ...
            'FontSize', 20)
        xlabel('True orientation (degree)')
        ylabel('Estimated orientation (degree)')        
        
    end
end

% Plot
figure
xRange = [-22 22];
indX = find(thetaStim >= xRange(1) & thetaStim <= xRange(2));
xMax = length(indX);
yRange = [-35 35];
indY = find(th >= yRange(1) & th <= yRange(2));
yMax = length(indY);
thNew = th(indY);
indYStart = find(thNew == xRange(1));
indYEnd = find(thNew == xRange(2));

imagesc(pthhGthAverage(indY, indX));
hold on;
axis xy;
colormap('gray');
plot([1 xMax],[round(yMax/2) round(yMax/2)],'w:', 'LineWidth', 1.5);
plot([round(xMax/2) round(xMax/2)],[1 yMax],'w:', 'LineWidth', 1.5);
plot([1 xMax],[indYStart indYEnd],'b:', 'LineWidth', 2);
set(gca, 'ylim', [1 yMax], 'xlim', [1 xMax], ...
    'XTick', round(linspace(1,xMax,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
    'YTick', round(linspace(1,yMax,5)), 'YTickLabel', num2cell(round(linspace(yRange(1),yRange(2),5))), ...
    'FontSize', 20)
xlabel('True orientation (degree)')
ylabel('Estimated orientation (degree)')        

mglGetSecs - startTime