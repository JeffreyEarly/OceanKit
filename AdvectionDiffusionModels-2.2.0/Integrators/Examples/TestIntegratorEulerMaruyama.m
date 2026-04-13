% Manual example for `IntegratorEulerMaruyama`.

reps = 500;
nDims = 2;
kappa = 1;

t = (0:5:900)';
deltaT = 1;

f = @(t, x) zeros(size(x));
g = @(t, x) sqrt(2*kappa)*ones(size(x));

p0 = 0;
x0 = p0*ones(reps, nDims);

integrator = IntegratorEulerMaruyama(f, g, x0, dt=deltaT);
pn = integrator.integrateToTime(t);
x = squeeze(pn(:,1,:)).';
y = squeeze(pn(:,2,:)).';

D2 = x(end,:).^2 + y(end,:).^2;
kappa_out = mean(D2)/(4*t(end))

figure
plot(t, x)

tau = 60;
sigma = sqrt(2*kappa)/tau;
f = @(t, x) cat(2, x(:,2), -x(:,2)/tau);
g = @(t, x) cat(2, zeros(size(x(:,1))), sigma*ones(size(x(:,2))));

x0 = zeros(reps, 2);

integrator = IntegratorEulerMaruyama(f, g, x0, dt=deltaT);
pn = integrator.integrateToTime(t);
x = squeeze(pn(:,1,:)).';
u = squeeze(pn(:,2,:)).';

D2 = x(end,:).^2;
kappa_out = mean(D2)/(2*t(end))

figure
plot(t, x)
