% Show the gridded-streamfunction estimator on synthetic strain-only trajectories.
sigma = 4e-6;
theta = pi/9;
zeta = 0;
u0 = 0.08;
v0 = -0.05;
kappa = 0;

velocityField = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta, u0=u0, v0=v0);
integrator = AdvectionDiffusionIntegrator(velocityField, kappa);

% Particles placed in a "plus" configuration.
x0 = [-500; -250; 0; 0; 0; 0; 0; 250; 500];
y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0];

% Use denser sampling than the direct linear-estimator example so the
% spline-based fit has enough support over the short synthetic window.
T = 12*3600;
dt = 900;
[t, x, y] = integrator.particleTrajectories(x0, y0, T, dt);

trajectories = TrajectorySpline.empty(0, 1);
for iDrifter = 1:size(x, 2)
    trajectories(end+1, 1) = TrajectorySpline.fromData(t, x(:, iDrifter), y(:, iDrifter), S=3);
end

fit = GriddedStreamfunction.fromTrajectories(trajectories);
tGrid = repmat(t, 1, size(x, 2));

uBackgroundFit = fit.uBackground(t);
vBackgroundFit = fit.vBackground(t);
uFit = fit.uMesoscale(tGrid, x, y) + repmat(uBackgroundFit, 1, size(x, 2));
vFit = fit.vMesoscale(tGrid, x, y) + repmat(vBackgroundFit, 1, size(x, 2));
sigmaNFit = fit.sigma_n(tGrid, x, y);
sigmaSFit = fit.sigma_s(tGrid, x, y);
zetaFit = fit.zeta(tGrid, x, y);

fprintf('Recovered background velocity samples\n');
fprintf('  rev3 note: background is the COM residual after the mesoscale solve,\n');
fprintf('             so constant-velocity truth need not match the fitted split.\n');
fprintf('  t (h)       fit u_bg       fit v_bg\n');
for iTime = 1:numel(t)
    fprintf(' %6.2f   %12.6f  %12.6f\n', t(iTime)/3600, uBackgroundFit(iTime), vBackgroundFit(iTime));
end

fprintf('\nRecovered linear diagnostics on observed trajectories\n');
fprintf('  quantity      true value      mean fitted      max abs error\n');
fprintf('  sigma_n    %12.6e  %12.6e  %12.6e\n', ...
    velocityField.sigma_n, mean(sigmaNFit, "all"), max(abs(sigmaNFit - velocityField.sigma_n), [], "all"));
fprintf('  sigma_s    %12.6e  %12.6e  %12.6e\n', ...
    velocityField.sigma_s, mean(sigmaSFit, "all"), max(abs(sigmaSFit - velocityField.sigma_s), [], "all"));
fprintf('  zeta       %12.6e  %12.6e  %12.6e\n', ...
    velocityField.zeta, mean(zetaFit, "all"), max(abs(zetaFit - velocityField.zeta), [], "all"));

fprintf('\nMax absolute trajectory-space errors\n');
fprintf('  u error        %12.6e m s^-1\n', max(abs(uFit - velocityField.u(tGrid, x, y)), [], "all"));
fprintf('  v error        %12.6e m s^-1\n', max(abs(vFit - velocityField.v(tGrid, x, y)), [], "all"));
fprintf('  sigma_n error  %12.6e s^-1\n', max(abs(sigmaNFit - velocityField.sigma_n), [], "all"));
fprintf('  sigma_s error  %12.6e s^-1\n', max(abs(sigmaSFit - velocityField.sigma_s), [], "all"));
fprintf('  zeta error     %12.6e s^-1\n', max(abs(zetaFit - velocityField.zeta), [], "all"));
