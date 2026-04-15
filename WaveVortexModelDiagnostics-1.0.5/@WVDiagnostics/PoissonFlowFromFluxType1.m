function [X,Y,U,V] = PoissonFlowFromFluxType1(wvd, flux)
% We will treat the first dimension as `x' and the second
% dimension as `y'. This means that the flux in the usual form,
% which is j by kRadial, might need to be transposed to get
% what you want.
%
% [X,Y,U,V] = WVDiagnostics.PoissonFlowFromFlux(wvt.kRadial,jWavenumber,flux.');
% quiver(X,Y,10*U,10*V,'off',Color=0*[1 1 1])
% For the DCT2/DST2 we use a half-shift grid
x = wvd.kRadial;
y = wvd.jWavenumber;

N = length(x);
M = length(y);
dk = 1/(2*(N-1)*(x(2)-x(1)));
dl = 1/(2*(M-1)*(y(2)-y(1)));

k=dk*(0:(N-1)).';
l=dl*(0:(M-1)).';

DCTx = WVDiagnostics.DCT1(N);
DCTy = WVDiagnostics.DCT1(M);

flux_ky = DCTx*flux;
flux_kl = shiftdim(DCTy*shiftdim(flux_ky,1),1);

[K,L] = ndgrid(k,l);

D = -((2*pi*K).^2 + (2*pi*L).^2);
D(1,1) = Inf;

UFactor = 2*pi*K./D;
VFactor = 2*pi*L./D;

iDCTx = WVDiagnostics.iDCT1(N);
iDSTy = WVDiagnostics.iDST1(M);
iDSTy = cat(2,zeros(M,1),iDSTy,zeros(M,1));

V_xl = iDCTx*(VFactor.*flux_kl);
V = shiftdim(iDSTy*shiftdim(V_xl,1),1);

iDCTy = WVDiagnostics.iDCT1(M);
iDSTx = WVDiagnostics.iDST1(N);
iDSTx = cat(2,zeros(N,1),iDSTx,zeros(N,1));

U_xl = iDSTx*(UFactor.*flux_kl);
U = shiftdim(iDCTy*shiftdim(U_xl,1),1);

[X,Y] = ndgrid(x,y);
end