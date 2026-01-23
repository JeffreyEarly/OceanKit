basedir = "/Users/Shared/CimRuns_June2025/output/";
basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";

%%
runNumber=1; runName = "hydrostatic: geostrophic";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","256") + ".nc");

%%
runNumber=9; runName = "hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","256") + ".nc");

%%
wvt = wvd.wvt;

U2 = wvd.spectrumWithFgTransform(wvt.u);
V2 = wvd.spectrumWithFgTransform(wvt.v);
N2 = wvd.spectrumWithGgTransform(wvt.eta);
% N2 = wvd.crossSpectrumWithGgTransform(wvt.eta,wvt.eta);

E = 0.5*(U2 + V2 + N2);

E2 = wvt.Apm_TE_factor.*( abs(wvt.Ap).^2 + abs(wvt.Am).^2 )  + wvt.A0_TE_factor.*( abs(wvt.A0).^2);

max(abs(E(:)-E2(:)))

%%
int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);
[u,v,ape,apv] = wvt.variableWithName('u','v','ape','apv');
ke = (u.^2 + v.^2)/2;
int_vol(ke + ape) - sum(E(:))
int_vol(ape) - sum(sum(0.5*N2(:,2:end)))
int_vol(ape) - sum(sum(0.5*N2(:,2:end)))

%%
rho_nm = wvt.chebfunForZArray(wvt.rho_nm)/wvt.rho0 - 1;
p_nm = - wvt.g * cumsum(rho_nm);
p_nm = p_nm - p_nm(0);
eta_true = wvt.eta_true;
Z = wvt.Z;
int_vol(wvt.g*eta_true.*rho_nm(Z - eta_true) + p_nm(Z) - p_nm(Z - eta_true))

%%
ape3 = zeros(wvt.spatialMatrixSize);
N = size(ape3,1)*size(ape3,2);
for iZ=1:length(wvt.z)
    zp = wvt.z(iZ);
    f = chebfun( @(eta) eta*N2t(zp-eta),[zp zp + wvt.Lz],'splitting','on');
    index_offset = N*(iZ-1);
    for index=1:N
        eta = wvt.eta_true(index + index_offset);
        ape3(index + index_offset) = sum(f,0,eta);
    end
end
% ape3 = shiftdim(ape3,1);
int_vol(ape3)

%%

N2t = -wvt.g*diff(rho_nm);
logN2t = log(N2t) - log(N2t(end));

ape2 = 0.5*N2t(wvt.Z).*logN2t(wvt.Z-wvt.eta_true).*wvt.eta_true.*wvt.eta_true;
int_vol(ape2)

dLogN2t = diff(logN2t);

zp = wvt.z(10);
f = chebfun( @(eta) eta*eta*dLogN2t(zp-eta),[zp zp + wvt.Lz],'splitting','on');

ape3 = zeros(wvt.spatialMatrixSize);
for index=1:numel(Z)
    z = Z(index);
    eta = eta_true(index);
    ape3(i) = sum(dLogN2t,0,);
end

%%

first_term = wvd.crossSpectrumWithGgTransform(wvt.g*eta_true./wvt.N2Function(Z),rho_nm(Z - eta_true));
second_term = wvd.crossSpectrumWithFgTransform(p_nm(Z) - p_nm(Z - eta_true),ones(size(Z)));
spectrum = first_term + second_term;
sum(spectrum(:))

int_vol(wvt.g*eta_true.*rho_nm(Z - eta_true)) - sum(first_term(:))
int_vol(p_nm(Z) - p_nm(Z - eta_true)) - sum(second_term(:))

%%
S_jk = wvt.transformToRadialWavenumber(spectrum);
figure
plot(wvt.kRadial,sum(S_jk,1))

%%
alt_first_term = 0.5*wvd.crossSpectrumWithGgTransform(eta_true.*wvt.N2Function(Z-eta_true)./wvt.N2Function(Z),eta_true);
sum(alt_first_term(:))
%%

prefactorJ = wvt.h_0; prefactorJ(1) = wvt.Lz;
prefactorK = 2*ones(1,wvt.Nkl); prefactorK(1) = 1;
prefactor = prefactorJ * prefactorK;

f_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(p_nm(Z) - p_nm(Z - eta_true)));
S_f = prefactor.*abs(f_bar).^2;
sum(S_f(:))


%%
runNumber=18; runName = "non-hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

%%
wvt = wvd.wvt;

U2 = wvd.spectrumWithFgTransform(wvt.u);
V2 = wvd.spectrumWithFgTransform(wvt.v);
W2 = wvd.spectrumWithGgTransform(wvt.w);
rho_nm = @(z) wvt.rhoFunction(z) - wvt.rho0;
rho_nm_z = @(z) -(wvt.rho0/wvt.g)*wvt.N2Function(z);
eta_nl = (wvt.g/wvt.rho0)*rho_nm(wvt.Z - wvt.eta_true) ./ wvt.N2Function(wvt.Z);
N2 = wvd.crossSpectrumWithGgTransform(wvt.eta_true,eta_nl);

E = 0.5*(U2 + V2 + W2 + N2);

E2 = wvt.Apm_TE_factor.*( abs(wvt.Ap).^2 + abs(wvt.Am).^2 )  + wvt.A0_TE_factor.*( abs(wvt.A0).^2);

sum(E(:))-sum(E2(:))