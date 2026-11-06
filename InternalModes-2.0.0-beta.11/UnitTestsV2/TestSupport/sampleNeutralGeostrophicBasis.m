function values = sampleNeutralGeostrophicBasis(basis,z)
% Differentiate independently sampled stored fields, without using their ODE.
% This helper is deliberately restricted to the spectral readiness study.
arguments
    basis
    z (:,1) double {mustBeReal,mustBeFinite}
end
if isstruct(basis)
    values = basis.evaluate(z);
    if ~isfield(values,"originalShape")
        values.originalShape = size(values.F);
    end
    return
end
solver = basis.solver;
if ~isa(solver,"IMSolverSpectral")
    error("NeutralGeostrophicReference:UnsupportedSolver","Readiness differentiation requires IMSolverSpectral.");
end
values.F = reshape(basis.F(z),numel(z),[]);
values.G = reshape(basis.G(z),numel(z),[]);
values.originalShape = size(basis.F(solver.zNative));
if isa(basis,"IMGeostrophicModesBasis")
    values.Fz = basis.uz(z);
    values.Gz = basis.Gz(z);
elseif isa(basis,"IMMeanDensityAnomalyModesBasis")
    % Integrating n coefficients stores n+1 pressure coefficients. Recover
    % that full polynomial; n native samples would alias T_n into T_(n-2).
    coefficients = sampledCoefficients(basis,"F",solver.nEVP+1);
    values.Fz = solver.evaluatePhysicalDerivative(coefficients,z,1);
    values.Gz = basis.Gz(z);
    factors = basis.normalizationFactors(basis.normalization);
    values.Gzz = solver.evaluatePhysicalDerivative(basis.nativeModes,z,2)./factors;
else
    % Legacy F-only zero-APV construction also integrates one additional
    % companion coefficient. The canonical paired solve has degree n-1.
    nSamples = solver.nEVP+double(isa(basis,"IMGeostrophicZeroAPVModesBasis"));
    FCoefficients = sampledCoefficients(basis,"F",nSamples);
    GCoefficients = sampledCoefficients(basis,"G",nSamples);
    values.Fz = solver.evaluatePhysicalDerivative(FCoefficients,z,1);
    values.Gz = solver.evaluatePhysicalDerivative(GCoefficients,z,1);
end
end

function coefficients = sampledCoefficients(basis,variable,nSamples)
solver = basis.solver;
if nSamples == solver.nEVP
    z = solver.zNative;
else
    bounds = [min(solver.xNative) max(solver.xNative)];
    x = mean(bounds)+(diff(bounds)/2)*cos((0:nSamples-1)'*pi/(nSamples-1));
    z = solver.zOfX(x);
    z(1) = solver.zDomain(2);
    z(end) = solver.zDomain(1);
end
if variable == "F"
    samples = basis.F(z);
else
    samples = basis.G(z);
end
samples = reshape(samples,nSamples,[]);
% The even-extension FFT is the DCT-I for the complete Lobatto polynomial.
extended = [samples;samples(end-1:-1:2,:)];
transformed = fft(extended,[],1)/(nSamples-1);
coefficients = transformed(1:nSamples,:);
if isreal(samples)
    coefficients = real(coefficients);
end
coefficients([1 end],:) = coefficients([1 end],:)/2;
end
