%% Exact geostrophic zero-APV modes
% Constant and exponential stratification have closed-form zero-APV
% solutions. The analytical API scales those solutions by their physical
% endpoint responses, so the surface column has response [1;0] and the
% bottom column has response [0;1].

clear
close all

D = 4000;
zDomain = [-D 0];
f0 = 1e-4;
g = 9.81;
k = 1e-4;
z = linspace(zDomain(1),zDomain(2),301).';

%% Constant stratification
% Requesting both endpoints returns two aligned F/G columns. For a free
% surface the surface response is G(0)-F(0); the bottom response is G(-D).

N0 = 5.2e-3;
constant = IMConstantStratificationSolution(N0=N0,zDomain=zDomain,f0=f0,g=g);
constantModes = constant.geostrophicZeroAPVModesAtWavenumber(k, ...
    endpoints=["surface" "bottom"],surfaceBoundary="freeSurface");

FConstant = constantModes.F(z);
GConstant = constantModes.G(z);

figure(Name="Exact constant-stratification zero-APV modes")
tiledlayout(1,2,TileSpacing="compact",Padding="compact")
nexttile
plot(FConstant(:,:,1),z,LineWidth=1.5)
xlabel("F")
ylabel("z (m)")
title("Pressure/velocity structure")
legend("surface response","bottom response",Location="best")
grid on
nexttile
plot(GConstant(:,:,1),z,LineWidth=1.5)
xlabel("G")
title("Displacement structure")
grid on

%% Endpoint choices and surface conventions
% A one-endpoint request still enforces zero response at the omitted
% endpoint. The surface convention changes only the definition of the
% surface response.

surfaceOnly = constant.geostrophicZeroAPVModesAtWavenumber(k, ...
    endpoints="surface",surfaceBoundary="rigidLid");
bottomOnly = constant.geostrophicZeroAPVModesAtWavenumber(k, ...
    endpoints="bottom",surfaceBoundary="rigidLid");

rigidSurfaceResponse = surfaceOnly.G(zDomain(2));
rigidSurfaceBottomResponse = surfaceOnly.G(zDomain(1));
rigidBottomSurfaceResponse = bottomOnly.G(zDomain(2));
rigidBottomResponse = bottomOnly.G(zDomain(1));

endpointCheck = table( ...
    [rigidSurfaceResponse; rigidBottomSurfaceResponse], ...
    [rigidSurfaceBottomResponse; rigidBottomResponse], ...
    RowNames=["surface mode" "bottom mode"], ...
    VariableNames=["surfaceResponse" "bottomResponse"]);
disp(endpointCheck)

%% Exponential stratification
% The exponential family evaluates scaled modified-Bessel functions and
% their exact derivatives. Its result has the same public contract.

b = 1300;
exponential = IMExponentialStratificationSolution(N0=N0,b=b,zDomain=zDomain,f0=f0,g=g);
exponentialModes = exponential.geostrophicZeroAPVModesAtWavenumber(k, ...
    endpoints=["surface" "bottom"],surfaceBoundary="freeSurface");

FExponential = exponentialModes.F(z);
GExponential = exponentialModes.G(z);

figure(Name="Exact exponential-stratification zero-APV modes")
tiledlayout(1,2,TileSpacing="compact",Padding="compact")
nexttile
plot(FExponential(:,:,1),z,LineWidth=1.5)
xlabel("F")
ylabel("z (m)")
title("Pressure/velocity structure")
legend("surface response","bottom response",Location="best")
grid on
nexttile
plot(GExponential(:,:,1),z,LineWidth=1.5)
xlabel("G")
title("Displacement structure")
grid on

%% Rotate the endpoint coordinates
% Canonical columns are convenient boundary coordinates. Applications may
% instead diagonalize generalized energy or surface buoyancy. The same
% rotation is applied to F and G.

g0 = -0.03;
gd = 0.01;
depthModes = exponentialModes.rotateBoundaryDepth(g0=g0,gd=gd);
surfaceModes = exponentialModes.rotateSurfaceBuoyancy(g0=g0,gd=gd);

rotationSummary = table( ...
    depthModes.rotationEigenvalues(:,1),depthModes.h0(:,1), ...
    surfaceModes.rotationEigenvalues(:,1), ...
    VariableNames=["boundaryDepthEigenvalue" "h0" "surfaceBuoyancyEigenvalue"]);
disp(rotationSummary)

%% Compare with a refined numerical solve
% The analytical object is an independent reference; it does not call a
% numerical mode solver.

N2 = @(zIn) N0*N0*exp(2*zIn/b);
problem = IMGeostrophicZeroAPVModes.atWavenumber( ...
    N2=N2,zDomain=zDomain,f0=f0,g=g,k=k, ...
    endpoints=["surface" "bottom"],surfaceBoundary="freeSurface");
numericalModes = IMSolverSpectral(nEVP=128).solveGeostrophicZeroAPVModes(problem);

relativeFError = norm(numericalModes.F(z)-FExponential)/norm(FExponential);
relativeGError = norm(numericalModes.G(z)-GExponential)/norm(GExponential);
fprintf("Refined numerical relative F error: %.3e\n",relativeFError)
fprintf("Refined numerical relative G error: %.3e\n",relativeGError)
