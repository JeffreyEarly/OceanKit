
%% Set stratification and determine quadrature points
Lz = 4000;
N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
N2 = @(z) N0*N0*exp(2*z/L_gm);
% N2 = @(z) N0*N0*ones(size(z));
Nz = 64;
latitude = 31;
g = 9.81;

z_q = WVStratification.quadraturePointsForStratifiedFlow(Lz,Nz,N2=N2,latitude=latitude);

%% Now init the internal modes at those points

im = InternalModesWKBSpectral(N2=N2,zIn=[-Lz 0],zOut=z_q,latitude=latitude,nEVP=max(256,floor(2.1*Nz)),nModes = Nz-1, g=g);
im.normalization = Normalization.kConstant;
im.upperBoundary = UpperBoundary.rigidLid;

%%
k = 0;
[Finv,Ginv,h] = im.ModesAtWavenumber(k);

Nz = size(Finv,1);
nModes = size(Finv,2);

% Make these matrices invertible by adding the barotropic mode
% to F, and removing the boundaries of G.
Finv = cat(2,ones(Nz,1),Finv);
h_F = cat(2,Lz,h);
Ginv = Ginv(2:end-1,1:end-1);
h_G = h(1:end-1);

% Compute the precondition matrices (really, diagonals)
P = max(abs(Finv),[],1); % ones(1,size(Finv,1)); %
Q = max(abs(Ginv),[],1); % ones(1,size(Ginv,1)); %

% Now create the actual transformation matrices
PFinv = Finv./P;
QGinv = Ginv./Q;
PF = inv(PFinv);
QG = inv(QGinv);

b = zeros(Nz,1);
b(1) = Lz;
dz = (PFinv.')\b;

%% Let's see how well modes at wavenumber k do, at these points
k = 2*pi/500;
[~,Ginv_zq_k] = im.ModesAtWavenumber(k);
Ginv_zq_k = Ginv_zq_k(2:end-1,1:end-1);
Q_zq_k = max(abs(Ginv_zq_k),[],1);
QGinv_zq_k = Ginv_zq_k./Q_zq_k;

[cond(QGinv), cond(QGinv_zq_k)]

%% Now find and init the internal modes at the k-quadrature points
z_k = im.GaussQuadraturePointsForModesAtWavenumber(Nz,k);
im_k = InternalModesWKBSpectral(N2=N2,zIn=[-Lz 0],zOut=z_k,latitude=latitude,nEVP=max(256,floor(2.1*Nz)),nModes = Nz-1, g=g);
im_k.normalization = Normalization.kConstant;
im_k.upperBoundary = UpperBoundary.rigidLid;

%%
[~,Ginv_k,h_k] = im_k.ModesAtWavenumber(k);
Ginv_k = Ginv_k(2:end-1,1:end-1);
Q_k = max(abs(Ginv_k),[],1);
QGinv_k = Ginv_k./Q_k;
QG_k = inv(QGinv_k);

%%
[cond(QGinv), cond(QGinv_k)]
figure, plot(QGinv(:,12),z_q(2:end-1)), hold on, plot(QGinv_k(:,12),z_k(2:end-1))

%%
K=8;
spline = BSpline(K,BSpline.knotPointsForDataPoints(z_q,K=K));
% these matrices map coefficients to interpolated points
Z_q = BSpline.Spline( z_q, spline.t_knot, spline.K );
Z_k = BSpline.Spline( z_k, spline.t_knot, spline.K );
% Pr maps from the z_k to the z_points
Pr = Z_k*inv(Z_q);

%% plot
iMode = 50;
[~,Ginv_zq_k] = im.ModesAtWavenumber(k);
mode = Ginv_zq_k(:,iMode)./Q_k(iMode);
figure, plot(QGinv_k(:,iMode),z_k(2:end-1)), hold on, plot(Pr*mode,z_k)

% take all the k-modes on the zq grid, and project:
A = Pr*Ginv_zq_k(:,1:end-1);
A = A(2:end-1,:);

error = QG_k*(Ginv_k./Q_k - A./Q_k);

%%
figure
plot(1:size(error,2),diag(error)), yscale('log'), hold on
plot(1:size(error,2),vecnorm(error))

%%
condNumber = zeros(size(Ginv_zq_k,2),1);
for iMode=1:size(Ginv_zq_k,2)
    condNumber(iMode) = cond(Ginv_zq_k(:,1:iMode));
end
figure, plot(1:length(condNumber),condNumber), yscale('log')

%%
r = 60;
[~,Ginv_zq_k] = im.ModesAtWavenumber(k);
Ginv_zq_k = Ginv_zq_k(2:end-1,1:end-1);
Q_zq_k = max(abs(Ginv_zq_k),[],1);
QGinv_zq_k = Ginv_zq_k./Q_zq_k;

QGinv_zq_k_r = QGinv_zq_k(:,1:r);
M = (N2(z_q(2:end-1)) - im.f0*im.f0)/g;
WG = dz(2:end-1) .* M .* QGinv_zq_k_r;
A = (QGinv_zq_k_r.' * WG) \ (WG.');


A = cat(1,A,zeros(size(Ginv_k,1)-r,size(A,2)));

% compare the projection

% [~,Ginv_zq_k] = im.ModesAtWavenumber(k);
error = QG_k*(Ginv_k./Q_k) - A*QGinv_zq_k;
figure
plot(1:size(error,2),diag(abs(error))), yscale('log'), hold on
plot(1:size(error,2),vecnorm(error))

% The important conclusion here is that we get a damn good reproduction of
% the spectrum with this approach and that the condition number is, in
% fact, telling us something very useful about how well it's doing (or will do!).