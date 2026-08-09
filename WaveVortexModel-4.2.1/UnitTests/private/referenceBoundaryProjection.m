function [directFp,directFm,fourierFp,fourierFm,gBottom,gBottomFourier] = referenceBoundaryProjection(wvt,topographicHeight,barotropicVelocity)
%REFERENCEBOUNDARYPROJECTION Evaluate the bottom-pressure projection directly.

gBottom = barotropicVelocity(1)*wvt.diffX(topographicHeight)+barotropicVelocity(2)*wvt.diffY(topographicHeight);
gBottomTransform = wvt.transformFromSpatialDomainWithFourier(repmat(gBottom,1,1,wvt.Nz));
gBottomFourier = gBottomTransform(1,:);
directFp = complex(zeros(size(wvt.Ap)));
directFm = complex(zeros(size(wvt.Am)));
fourierFp = complex(zeros(size(wvt.Ap)));
fourierFm = complex(zeros(size(wvt.Am)));
[~,iBottom] = min(wvt.z);
x = reshape(wvt.x,[],1);
y = reshape(wvt.y,1,[]);

for index = reshape(find(wvt.waveComponent.maskAp),1,[])
    [~,iHorizontal] = ind2sub(size(wvt.Ap),index);
    coefficient = complex(zeros(size(wvt.Ap)));
    coefficient(index) = wvt.NAp(index);
    pressure = wvt.g*wvt.transformToSpatialDomainWithF(Apm=coefficient);
    pressureFourier = wvt.transformFromSpatialDomainWithFourier(pressure);
    pressureBottom = pressureFourier(iBottom,iHorizontal)*wvt.phase(index);
    pressurePlane = pressureBottom*exp(1i*(wvt.K(index)*x+wvt.L(index)*y));
    directFp(index) = mean(conj(pressurePlane).*gBottom,"all")/wvt.Apm_TE_factor(index);
    fourierFp(index) = conj(pressureBottom)*gBottomFourier(iHorizontal)/wvt.Apm_TE_factor(index);
end

for index = reshape(find(wvt.waveComponent.maskAm),1,[])
    [~,iHorizontal] = ind2sub(size(wvt.Am),index);
    coefficient = complex(zeros(size(wvt.Am)));
    coefficient(index) = wvt.NAm(index);
    pressure = wvt.g*wvt.transformToSpatialDomainWithF(Apm=coefficient);
    pressureFourier = wvt.transformFromSpatialDomainWithFourier(pressure);
    pressureBottom = pressureFourier(iBottom,iHorizontal)*wvt.conjPhase(index);
    pressurePlane = pressureBottom*exp(1i*(wvt.K(index)*x+wvt.L(index)*y));
    directFm(index) = mean(conj(pressurePlane).*gBottom,"all")/wvt.Apm_TE_factor(index);
    fourierFm(index) = conj(pressureBottom)*gBottomFourier(iHorizontal)/wvt.Apm_TE_factor(index);
end
end
