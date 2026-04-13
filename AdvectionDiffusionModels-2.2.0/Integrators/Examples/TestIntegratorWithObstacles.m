% Manual example for `IntegratorWithObstacles`.

obstacle = polyshape([0 1 1 0], [0 0 1 1]);
y0 = [-0.5 0.5; -0.5 1.5; -1.5 0.25];
t = (0:0.25:2).';

f = @(t, y) repmat([2 0], size(y, 1), 1);

integrator = IntegratorWithObstacles(f, y0, dt=0.25, obstacles=obstacle);

y = integrator.integrateToTime(t);
x = squeeze(y(:,1,:)).';
z = squeeze(y(:,2,:)).';

figure
plot(obstacle)
hold on
plot(x, z, 'LineWidth', 1.5)
scatter(y0(:,1), y0(:,2), 36, 'filled')
axis equal
xlabel('x')
ylabel('y')
title('IntegratorWithObstacles example')
