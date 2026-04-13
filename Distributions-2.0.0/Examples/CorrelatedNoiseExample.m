%% Correlated Gaussian noise
sigma = 1;
distribution = NormalDistribution(sigma=sigma);

n = 1000;
z = distribution.rand([n 1]);

figure
histogram(z, 50, 'Normalization', 'pdf')
zdist = linspace(min(z), max(z), 100)';
hold on
plot(zdist, distribution.pdf(zdist), 'LineWidth', 2)
title('NormalDistribution samples and pdf')
legend('samples', 'pdf')

%% Correlate the draws with a Gaussian autocorrelation model
tau = 10;
distribution.rho = @(z) exp(-(z/tau).^2);

t = (0:500).';
z = distribution.noise(t);

figure
plot(t, z)
xlabel('time')
ylabel('signal')
title('Correlated Gaussian noise')

%% Compare the empirical and theoretical autocorrelation
lagAxis = t(1:find(t < 3*tau, 1, 'last'));
zCentered = z - mean(z);
ac = zeros(size(lagAxis));
for iLag = 1:length(lagAxis)
    lag = iLag - 1;
    ac(iLag) = sum(zCentered(1:end-lag).*zCentered(1+lag:end))/sum(zCentered.^2);
end
dof = length(z) - (0:length(lagAxis)-1)';

SE_indep = lagAxis(2:end);
SE = sqrt((1 + 2*cumsum(ac.^2))./dof);
SE(1) = sqrt(1/dof(1));
SE(end) = [];

figure
tiledlayout(2, 2)

nexttile
plot(t, z)
title(sprintf('mean=%.2f, std=%.2f', mean(z), std(z)))
xlabel('time')

nexttile
histogram(z, 100, 'Normalization', 'pdf')
hold on
zdist = linspace(min(z), max(z), 100)';
plot(zdist, distribution.pdf(zdist))
title('pdf')

nexttile
plot(lagAxis, ac)
hold on
plot(lagAxis, distribution.rho(lagAxis))
plot(SE_indep, [3*SE, -3*SE], 'LineWidth', 1.5, 'Color', 0.4*[1 1 1])
legend('actual', 'theoretical', '3\sigma standard error')
xlabel('lag')
ylabel('autocorrelation')
ylim([0 1])

nexttile
histogram(z, 100, 'Normalization', 'cdf')
hold on
plot(zdist, distribution.cdf(zdist))
title('cdf')
