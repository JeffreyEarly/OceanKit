%% Set stratification and determine quadrature points
Lz = 4000;
N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
N2 = @(z) N0*N0*exp(2*z/L_gm);
% N2 = @(z) N0*N0*ones(size(z));
Nz = 512;
latitude = 31;
g = 9.81;

z = WVStratification.quadraturePointsForStratifiedFlow(Lz,Nz,N2=N2,latitude=latitude);

%% Now init the internal modes at those points

im = InternalModesWKBSpectral(N2=N2,zIn=[-Lz 0],zOut=z,latitude=latitude,nEVP=max(256,floor(2.1*Nz)),nModes = Nz-1, g=g);
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

%% This estimates the quadrature weights just using the grid
z_mid = z(1:end-1) + diff(z)/2;
z2 = cat(1,z(1),z_mid,z(end));
dz = diff(z2);

% Now apply the approximated inverse
M = (N2(z(2:end-1)) - im.f0*im.f0)/g;
Icheck = Ginv.' * ( dz(2:end-1) .* M .* Ginv);
T = eye(Nz-2);
disp(norm(Icheck - T, 'fro')/norm(T,'fro'))

%% This determines the weights using the relation that \int F dz = 0
b = zeros(Nz,1);
b(1) = Lz;
w = (PFinv.')\b;

% Now apply the approximated inverse
M = (N2(z(2:end-1)) - im.f0*im.f0)/g;
Icheck = Ginv.' * ( w(2:end-1) .* M .* Ginv);
T = eye(Nz-2);
disp(norm(Icheck - T, 'fro')/norm(T,'fro'))

%% This determines the weight using the orthogonality relation
M = (N2(z(2:end-1)) - im.f0*im.f0)/g;
S = (M .* Ginv) * Ginv.';
w2 = 1 ./ diag(S);

Icheck = Ginv.' * ( w2 .* M .* Ginv);
T = eye(Nz-2);
disp(norm(Icheck - T, 'fro')/norm(T,'fro'))

%% This determines the weight using the orthogonality relation, from least-squares
G = Ginv;
M = (N2(z(2:end-1)) - im.f0*im.f0)/g;
N = Nz-2;

% Build target H (either eye(N) if orthonormalized, else your computed H)
T = eye(N);  % or eye(N)

% Build system A*W = b from upper triangle (m<=n)
idx = 0;
A = zeros(N*(N+1)/2, N);
b = zeros(N*(N+1)/2, 1);

for m = 1:N
  for n = m:N
    idx = idx + 1;
    A(idx,:) = ( M.* G(:,m) .* G(:,n)).';   % row vector length N
    b(idx)   = T(m,n);
  end
end

w3 = A \ b;                % or lsqminnorm(A,b) if near rank-deficient

Icheck = G.' * (M.* w3 .* G);
disp(norm(Icheck - T, 'fro')/norm(T,'fro'))

% CONCLUSION: They all produce the same relative error, but the first and
% second methods are better because they works for the first and last point
% of the F modes.

%%

M = (N2(z(2:end-1)) - im.f0*im.f0)/g;
G = ( w(2:end-1) .* M .* Ginv).';
I = G*Ginv;

QG2 = ( (w(2:end-1) .* M .* Ginv) .* Q).';
disp(norm(QG2*QGinv - T, 'fro')/norm(T,'fro'))
disp(norm(QG*QGinv - T, 'fro')/norm(T,'fro'))

a1 = vecnorm(QG*QGinv-T,Inf,1);
a2 = vecnorm(QG2*QGinv-T,Inf,1);

% [cond(PFinv), cond(QGinv), cond(PF), cond(QG)]

%%

tic
WG = w(2:end-1) .* M .* QGinv;
QG2 = (QGinv.' * WG) \ (WG.');
toc

tic
QG = inv(QGinv);
toc


% disp(norm(QG2*QGinv - T, 'fro')/norm(T,'fro'))