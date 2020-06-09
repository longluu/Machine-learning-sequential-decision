%%%%%%%%%%% Plot log likelihood difference of resample vs. non-resample for correct trials %%%%%%%%%%%
% Negative log likelihood taken from the fit
neg_llh_noResample = [4539.91 5901.30 6051.89 6413.41 6429.30 30138.46 ...
                      4799.79 5045.18 4923.35 4403.22 4446.10 5528.20 5131.05 38739.74];
                  
neg_llh_Resample =   [4533.30 5887.88 6039.95 6409.87 6405.35 30131.41 ...
                      4820.28 5054.09 4921.72 4398.48 4400.77 5534.17 5134.66 38872.65];

% Normalize by the number of trials
neg_llh_noResample(1:5) = neg_llh_noResample(1:5) / 2100;
neg_llh_noResample(6) = neg_llh_noResample(6) / (5*2100);
neg_llh_noResample(7:end-1) = neg_llh_noResample(7:end-1) / 1820;
neg_llh_noResample(end) = neg_llh_noResample(end) / (7*1820);

neg_llh_Resample(1:5) = neg_llh_Resample(1:5) / 2100;
neg_llh_Resample(6) = neg_llh_Resample(6) / (5*2100);
neg_llh_Resample(7:end-1) = neg_llh_Resample(7:end-1) / 1820;
neg_llh_Resample(end) = neg_llh_Resample(end) / (7*1820);

llh_diff = -neg_llh_Resample + neg_llh_noResample;                 
figure
hold on
bar(1:length(neg_llh_Resample), llh_diff)
set(gca, 'XTick', 1:length(neg_llh_Resample), ...
            'XTickLabel', num2cell(1:length(neg_llh_Resample)))
        
xlabel('Subject')
ylabel('Log likelihood difference (resample - no resample)')