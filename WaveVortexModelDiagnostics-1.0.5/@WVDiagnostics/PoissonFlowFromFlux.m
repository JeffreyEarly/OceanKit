function [X,Y,U,V] = PoissonFlowFromFlux(wvd, flux)
% We will treat the first dimension as `x' and the second
% dimension as `y'. This means that the flux in the usual form,
% which is j by kRadial, might need to be transposed to get
% what you want.
%
% [X,Y,U,V] = WVDiagnostics.PoissonFlowFromFlux(wvt.kRadial,jWavenumber,flux.');
% quiver(X,Y,10*U,10*V,'off',Color=0*[1 1 1])

% compute flux in mode space
x = wvd.kRadial/wvd.wvt.dk + 0*1/2; % this gives horizontal mode number corresponding to kRadial
y = wvd.j + 0*1/2;

% move the first row and column half an increment to finite wavenumber for display.
x(1) = x(2)/2;
y(1) = y(2)/2;

[X,Y,U,V] = WVDiagnostics.PoissonFlowFromFluxWithAxes(x,y,flux);

% scale for different aspect ratio in X/Y direction:
V = V.* ((max(wvd.kRadial)/max(wvd.jWavenumber))/(max(x)/max(y)));


end