%%%%%%%%%%% Plot -LLH of resample vs. non-resample for correct trials %%%%%%%%%%%
neg_llh_noResample = [4539.91 5901.30 6051.89 6413.41 6429.30 30138.46 ...
                      4799.79 5045.18 4923.35 4403.22 4446.10 5528.20 5131.05 38739.74];
                  
neg_llh_Resample =   [4533.30 5887.88 6039.95 6409.87 6405.35 30131.41 ...
                      4820.28 5054.09 4921.72 4398.48 4400.77 5534.17 5134.66 38872.65];      
                  
figure
hold on
bar(1:length(neg_llh_Resample), neg_llh_noResample ./ neg_llh_Resample)
plot(1:length(neg_llh_Resample), ones(1,length(neg_llh_Resample)), 'k--')
set(gca, 'XTick', 1:length(neg_llh_Resample), ...
            'XTickLabel', num2cell(1:length(neg_llh_Resample)))
ylim([0.99 1.011])
xlabel('Subject')
ylabel('Ratio -LLH no-resample vs. resample')