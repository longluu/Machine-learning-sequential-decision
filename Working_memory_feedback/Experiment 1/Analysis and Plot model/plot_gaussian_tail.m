%%% Plot the tails of Gaussian and see the estimate if the likelihood
%%% center is farther away from the boundary
x = -20:0.001:50;
sigma = 1;
mean_test = 0:2:20;
y_test = NaN(length(mean_test), length(x));
for i = 1 : length(mean_test)
    full_pdf = normpdf(x, mean_test(i), sigma);
    chopped_pdf = full_pdf;
    chopped_pdf(x>0) = 0;
    chopped_pdf = chopped_pdf / sum(chopped_pdf);
    y_test(i, :) = chopped_pdf;
end

theta_est = x * y_test';

%% Plot
nLines = length(mean_test);
legend_str = cell(nLines,1);

for ii=1:nLines
    legend_str{ii} = num2str(mean_test(ii));
end

figure
subplot(1, 2, 1)
plot(x, y_test)
xlabel('theta')
ylabel('p')
xlim([-2 0])
legend(legend_str,'location','NorthWest')

subplot(1, 2, 2)
plot(mean_test, theta_est, 'o-')
xlabel('Likelihood center')
ylabel('Estimate')
