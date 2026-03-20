function [X,Y,U,V] = PoissonFlowFromFlux(wvd, flux)
% We will treat the first dimension as `x' and the second
% dimension as `y'. This means that the flux in the usual form,
% which is j by kRadial, might need to be transposed to get
% what you want.
%
% [X,Y,U,V] = WVDiagnostics.PoissonFlowFromFlux(wvt.kRadial,jWavenumber,flux.');
% quiver(X,Y,10*U,10*V,'off',Color=0*[1 1 1])
% For the DCT2/DST2 we use a half-shift grid
x = wvd.kRadial + 0*(wvd.kRadial(2)-wvd.kRadial(1))/2;
y = wvd.jWavenumber + 0*(wvd.jWavenumber(2)-wvd.jWavenumber(1))/2;
[X,Y,U,V] = WVDiagnostics.PoissonFlowFromFluxWithAxes(x,y,flux);

U = U/wvd.kRadial(2);
V = V/wvd.jWavenumber(2);
end