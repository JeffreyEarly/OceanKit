%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Show how the linearized ellipse fit recovers strain and diffusivity from
% synthetic strain-diffusive trajectories.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma = 8e-6;
theta = pi/4;
zeta = 0;
kappa = 1;
T = 3*86400;
dt = 4*3600;

velocityField = LinearVelocityField(sigma=sigma, theta=theta, zeta=zeta);
integrator = AdvectionDiffusionIntegrator(velocityField,kappa);

x0 = 4*[-500; -250; 0; 0; 0; 0; 0; 250; 500];
y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0];

totalIterations = 100;
kappaEst = zeros(totalIterations,1);
sigmaEst = zeros(totalIterations,1);
thetaEst = zeros(totalIterations,1);
zetaEst = zeros(totalIterations,1);

for iIteration = 1:totalIterations
    [t,x,y] = integrator.particleTrajectories(x0,y0,T,dt);

    parameters = FitTrajectoriesToLinearizedEllipseModel(x, y, t, 'strain-diffusive');

    kappaEst(iIteration) = parameters.kappa;
    sigmaEst(iIteration) = parameters.sigma;
    thetaEst(iIteration) = mod(parameters.theta + pi/2,pi) - pi/2;
    zetaEst(iIteration) = parameters.zeta;
end

fprintf('expected (kappa,sigma,theta,zeta) = (%.2f, %.2g, %.1f, %.2g), found (kappa,sigma,theta,zeta) = (%.2f±%.2f, %.2g±%.2g, %.1f±%.1f, %.2g±%.2g)\n',kappa,sigma,theta, zeta, mean(kappaEst), std(kappaEst), median(sigmaEst), std(sigmaEst), mean(thetaEst), std(thetaEst), median(zetaEst), std(zetaEst));

sigma_n = sigmaEst.*cos(2*thetaEst);
sigma_s = sigmaEst.*sin(2*thetaEst);

data = cat(2,sigma_n,sigma_s);
densityModel = KernelDensityEstimate.fromData(data);
[density, gridVectors] = densityModel.densityOnGrid();

figure
contourf(gridVectors{1}, gridVectors{2}, density.');
hold on
for iCircle = 2:2:10
    rectangle('Position',[-iCircle -iCircle 2*iCircle 2*iCircle]*1e-6, 'Curvature', [1 1]);
end
scatter(sigma_n, sigma_s, 5^2, 0*[1 1 1],'filled');
scatter(sigma*cos(2*theta), sigma*sin(2*theta), 20^2, 1*[1 1 1],'LineWidth',5);
scatter(sigma*cos(2*theta), sigma*sin(2*theta), 20^2, 0*[1 1 1],'LineWidth',3);
axis equal
xlabel('\sigma_n')
ylabel('\sigma_s')

cmap = colormap;
cmap(1,:) = 1;
colormap(cmap)

figure
subplot(1,3,1)
histogram(kappaEst,'Normalization','pdf')
xline(kappa,'k--')
xlabel('\kappa')

subplot(1,3,2)
histogram(sigmaEst,'Normalization','pdf')
xline(sigma,'k--')
xlabel('\sigma')

subplot(1,3,3)
histogram(thetaEst*180/pi,'Normalization','pdf')
xline(theta*180/pi,'k--')
xlabel('\theta (degrees)')
