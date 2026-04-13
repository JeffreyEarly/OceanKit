%% Build a simple mixture distribution
pct = 0.05;
distribution = AddedDistribution(scalings=pct, distributions=[NormalDistribution(sigma=800); StudentTDistribution(sigma=8.5, nu=4.5)]);

z = linspace(-500, 500, 1000)';

figure
tiledlayout(2, 1)

nexttile
plot(z, distribution.pdf(z), 'LineWidth', 1.5)
xlabel('z')
ylabel('pdf')
title('Mixture distribution pdf')

nexttile
plot(z, distribution.cdf(z), 'LineWidth', 1.5)
xlabel('z')
ylabel('cdf')
title('Mixture distribution cdf')

%% Convert the 1D model into a 2D radial-distance distribution
distDistribution = TwoDimDistanceDistribution(distribution1d=distribution);

figure
plot(z, distDistribution.cdf(z), 'LineWidth', 1.5)
xlabel('radius')
ylabel('cdf')
title('TwoDimDistanceDistribution cdf')
