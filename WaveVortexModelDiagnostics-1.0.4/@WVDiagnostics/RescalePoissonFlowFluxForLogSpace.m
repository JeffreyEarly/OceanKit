function [logX,logY,Uprime,Vprime] = RescalePoissonFlowFluxForLogSpace(wvd,X,Y,U,V,options)
% Rescale Poisson Flow Flux For Log Space.
%
% RescalePoissonFlowFluxForLogSpace is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Diagnostics — Flux diagnostics — General — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
% - Declaration: [logX,logY,Uprime,Vprime] = RescalePoissonFlowFluxForLogSpace(wvd,X,Y,U,V,options)
% - Parameter wvd: WVDiagnostics object
% - Parameter X: input argument `X`
% - Parameter Y: input argument `Y`
% - Parameter U: input argument `U`
% - Parameter V: input argument `V`
% - Parameter shouldOnlyRescaleDirection: (optional) input argument `shouldOnlyRescaleDirection` (default: true)
% - Returns logX: output value `logX`
% - Returns logY: output value `logY`
% - Returns Uprime: output value `Uprime`
% - Returns Vprime: output value `Vprime`
arguments
    wvd
    X
    Y
    U
    V
    options.shouldOnlyRescaleDirection logical = true
end
logX = log10(X);
logY = log10(Y);
if options.shouldOnlyRescaleDirection == true
    r = 1;
else
    r = sqrt( (1./log10(X)).^2 + (1./log10(Y)).^2 );
end
theta = atan2(abs(log10(X)),abs(log10(Y)));
Uprime = U.*r.*cos(theta);
Vprime = V.*r.*sin(theta);
end