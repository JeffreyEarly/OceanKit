%% Compare the numerically constructed radial law to Rayleigh theory
sigma = 1;
normal = NormalDistribution(sigma=sigma);
normal2dDistance = TwoDimDistanceDistribution(distribution1d=normal);
rayleigh = RayleighDistribution(sigma=sigma);

z = linspace(0, 5, 1000)';

figure
plot(z, normal2dDistance.pdf(z), 'LineWidth', 1.5)
hold on
plot(z, rayleigh.pdf(z), '--', 'LineWidth', 1.5)
xlabel('radius')
ylabel('pdf')
legend('TwoDimDistanceDistribution', 'RayleighDistribution')
title('Radial pdf comparison')

figure
tiledlayout(2, 1)

nexttile
plot(z, 1 - normal2dDistance.cdf(z)./rayleigh.cdf(z))
xscale('log')
xlabel('radius')
ylabel('1 - F_{2D}/F_{Rayleigh}')
title('Relative cdf difference')

nexttile
plot(z, 1 - normal2dDistance.pdf(z)./rayleigh.pdf(z))
xscale('log')
xlabel('radius')
ylabel('1 - p_{2D}/p_{Rayleigh}')
title('Relative pdf difference')
