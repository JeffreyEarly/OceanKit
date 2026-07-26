function zIsopycnal = placeParticlesOnIsopycnal(wvt,x,y,zNoMotion)
% places Lagrangian particles along a specified isopycnal
%
% Given particle position (x,y), `zNoMotion` is used to determine the target
% isopycnal using the no-motion density,
%
% ```matlab
% targetRho = wvt.rhoFunction(zNoMotion);
% ```
%
% and a minimization algorithm is used to find zIsopycnal such that
%
% ```matlab
% zIsopycnal = rho(targetRho);
% ```
%
% where rho is the current total density field of the fluid.
%
% Note that the density is not necessarily monotonic, so the answer is not
% necessarily unique.
%
% - Topic: Utilities
% - Declaration: zIsopycnal = placeParticlesOnIsopycnal(x,y,z)
% - Parameter x: array of particle x positions
% - Parameter y: array of particle y positions
% - Parameter z: array of z location of target isopycnal in the no-motion profile
% - Returns zIsopycnal: particle depth 
arguments
    wvt 
    x 
    y 
    zNoMotion 
end

% the target isopycnal will be taken as the isopycnal at that depth in the
% no-motion state.
targetRho = wvt.rho0 - wvt.rhoFunction(zNoMotion);

% Now create an interpolat for the current total density field, and pad it
% so that particles at the boundaries get the correct interpolation.
U = wvt.rho0 - wvt.rho_total;
S = 4;
ix = mod((-S:(wvt.Nx+S-1)), wvt.Nx) + 1;
iy = mod((-S:(wvt.Ny+S-1)), wvt.Ny) + 1;
U_padded = U(ix, iy, :);

dx = wvt.x(2)-wvt.x(1);
dy = wvt.y(2)-wvt.y(1);
x_padded = (wvt.x(1) - S*dx : dx : wvt.x(end) + S*dx).';
y_padded = (wvt.y(1) - S*dy : dy : wvt.y(end) + S*dy).';
[X,Y,Z] = ndgrid(x_padded,y_padded,wvt.z);

rho_interpolant = griddedInterpolant(X,Y,Z,U_padded);

% make sure the particles are wrapped to the correct domain
x = mod(x,wvt.Lx);
y = mod(y,wvt.Ly);

% now use a simple bisection algorithm to zero-in on the right location.
% Our density is not necessarily monotonic, so the answer may not be
% unique.
tol = 1e-12;
a = -wvt.Lz*ones(size(x));
b = zeros(size(x));
c = (a + b)/2;

while ( norm(b - a, inf) >= tol )
    vals = rho_interpolant(x,y,c);
    % Bisection:
    I1 = ((vals-targetRho) <= -tol);
    I2 = ((vals-targetRho) >= tol);
    I3 = ~I1 & ~I2;
    a = I1.*c + I2.*a + I3.*c;
    b = I1.*b + I2.*c + I3.*c;
    c = (a+b)/2;
end

zIsopycnal = c;
% rhoIsopycnal = vals;

end