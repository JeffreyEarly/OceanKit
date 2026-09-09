%% Vertical spectra and component roles
% This example demonstrates the vertical spectrum API. Geostrophic F and G
% spectra are canonical because both components come from Sturm-Liouville
% projection operators. For nonzero-kappa IGW modes, G is canonical but F
% is numerical-only.

%% Build geostrophic modes
Lz = 4000;
N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);

latitude = 31;
g = 9.81;
nModes = 24;
nEVP = max(256,ceil(2.1*(nModes + 1)));
zDomain = [-Lz 0];
z = linspace(zDomain(1),zDomain(2),256).';

im = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=z,latitude=latitude,nEVP=nEVP,nModes=nModes,g=g);
im.normalization = Normalization.geostrophic;
im.upperBoundary = UpperBoundary.rigidLid;

geostrophicBasis = InternalModesBasis.fromSolverAtFrequency(im,0,nModes=nModes,g=g);
geostrophicTransform = geostrophicBasis.nativeTransform(component="both",projectionMethod="weightedPseudoinverse");

%% F and G spectra reproduce their Parseval weights
fCoefficients = exp(-(0:length(geostrophicTransform.spectralWeightsF)-1).'/6);
gCoefficients = exp(-(1:length(geostrophicTransform.spectralWeightsG)).'/6);

SF = geostrophicTransform.spectrum(fCoefficients,component="F");
SG = geostrophicTransform.spectrum(gCoefficients,component="G");

expectedSF = geostrophicTransform.spectralWeightsF .* abs(fCoefficients).^2;
expectedSG = geostrophicTransform.spectralWeightsG .* abs(gCoefficients).^2;

fprintf('Geostrophic F spectrum relative error: %.3e\n', norm(SF - expectedSF)/norm(expectedSF));
fprintf('Geostrophic G spectrum relative error: %.3e\n', norm(SG - expectedSG)/norm(expectedSG));
fprintf('Geostrophic F status: %s, component role: %s\n', geostrophicTransform.transformStatusF, geostrophicTransform.componentRoleF);
fprintf('Geostrophic G status: %s, component role: %s\n', geostrophicTransform.transformStatusG, geostrophicTransform.componentRoleG);

%% Wave G is canonical, wave F is numerical-only
kappa = 2*pi/150;
im.normalization = Normalization.kConstant;
waveBasis = InternalModesBasis.fromSolverAtWavenumber(im,kappa,nModes=nModes,g=g);
waveTransform = waveBasis.nativeTransform(component="both",projectionMethod="weightedPseudoinverse",allowNoncanonical=true);

waveGCoefficients = exp(-(1:length(waveTransform.spectralWeightsG)).'/5);
waveSG = waveTransform.spectrum(waveGCoefficients,component="G");
expectedWaveSG = waveTransform.spectralWeightsG .* abs(waveGCoefficients).^2;
fprintf('Wave G spectrum relative error: %.3e\n', norm(waveSG - expectedWaveSG)/norm(expectedWaveSG));
fprintf('Wave F status: %s, canonical flag: %d\n', waveTransform.transformStatusF, waveTransform.forwardProjectionAvailableF);

try
    waveFCoefficients = exp(-(1:length(waveTransform.spectralWeightsF)).'/5);
    waveTransform.spectrum(waveFCoefficients,component="F");
    error('VerticalSpectraAndComponentRoles:ExpectedFailure', 'Expected wave-F spectrum to fail.');
catch ME
    if ME.identifier ~= "InternalModesTransform:NoncanonicalSpectrum"
        rethrow(ME)
    end
    fprintf('Wave F canonical spectrum correctly unavailable: %s\n', ME.identifier);
end

%% Cross-spectra use the manuscript weights and horizontal multiplicity
otherGCoefficients = 0.5*gCoefficients;
SGCross = geostrophicTransform.crossSpectrum(gCoefficients,otherGCoefficients,component="G",horizontalMultiplicity=2);
expectedSGCross = 2*geostrophicTransform.spectralWeightsG .* real(gCoefficients .* conj(otherGCoefficients));
fprintf('Geostrophic G cross-spectrum relative error: %.3e\n', norm(SGCross - expectedSGCross)/norm(expectedSGCross));
