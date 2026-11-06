function [z, gridDesign] = pointsFromModeRoots(self, options)
% Return physical endpoints and roots of the next selected mode.
%
% Let `nModes=N` select the first $$N$$ columns of the ordered basis. The
% generating mode is the next selected column,
%
% $$
% u_{\mathrm{gen}}(z)=u_{\mathrm{selected},N+1}(z).
% $$
%
% This is a column index in the selected basis, not a physical mode label.
% In particular, `modeNumber(N+1)` need not equal $$N+1$$ when negative or
% zero modes precede the positive modes. The returned physical grid is
%
% $$
% \mathbf{z}_{\mathrm{root}}=\{z_b\}\cup
% \left\{z\in(z_b,z_s):u_{\mathrm{gen}}(z)=0\right\}\cup\{z_s\}.
% $$
%
% The interior roots are converted to physical $$z$$ coordinates, sorted,
% and deduplicated. Both physical endpoints are always included. If column
% $$N+1$$ is not already retained, the stored solver and EVP obtain one
% auxiliary mode and verify that the new solve reproduces the original
% selected prefix; the basis set itself is not changed.
%
% The result is a deterministic mode-root candidate grid, not a complete
% quadrature rule. Use `quadratureWeightsForPoints` to fit weights and assess
% how accurately the resulting rule reproduces the retained modal Gram
% matrix.
%
% ```matlab
% z = basisSet.pointsFromModeRoots(nModes=8);
% weights = basisSet.quadratureWeightsForPoints(z=z,nModes=8);
% transform = basisSet.discreteTransform(z=z,weights=weights,nModes=8);
% ```
%
% - Topic: Build discrete transforms
% - Declaration: z = pointsFromModeRoots(basisSet,options)
% - Parameter options.nModes: number of leading selected basis columns represented by the grid
% - Returns z: increasing physical grid containing both endpoints and the generating mode's interior roots
arguments
    self IMBasisSet
    options.nModes (1,1) double {mustBeInteger, mustBePositive} = size(self.nativeModes,2)
end

nModes = options.nModes;
nRetained = size(self.nativeModes,2);
if nModes > nRetained
    error("IMBasisSet:InvalidQuadratureModeCount", ...
        "The basis set contains %d retained modes, but nModes=%d was requested.", nRetained, nModes);
end

if nModes < nRetained
    generatingMode = self.nativeModes(:,nModes+1);
    generatingModeNumber = self.modeNumber(nModes+1);
else
    try
        auxiliaryBasis = self.solver.solveEVP(self.evp,nModes=nModes+1);
    catch exception
        wrapped = MException("IMBasisSet:AuxiliaryModeUnavailable", ...
            "The solver could not obtain the auxiliary mode required for a %d-mode quadrature grid.", nModes);
        wrapped = addCause(wrapped,exception);
        throwAsCaller(wrapped)
    end
    if size(auxiliaryBasis.nativeModes,2) < nModes+1
        error("IMBasisSet:AuxiliaryModeUnavailable", ...
            "The solver returned %d selected modes, but %d are required to generate this grid. Increase the solver resolution.", ...
            size(auxiliaryBasis.nativeModes,2), nModes+1);
    end

    candidateModeNumber = auxiliaryBasis.modeNumber(1:nModes);
    referenceEigenvalues = self.eigenvalues(1:nModes);
    candidateEigenvalues = auxiliaryBasis.eigenvalues(1:nModes);
    eigenvalueScale = max(1,max(abs([referenceEigenvalues candidateEigenvalues]),[],"all"));
    if ~isequal(candidateModeNumber,self.modeNumber(1:nModes)) ...
            || any(abs(candidateEigenvalues-referenceEigenvalues) > 1e-8*eigenvalueScale)
        error("IMBasisSet:AuxiliaryModeMismatch", ...
            "The auxiliary solve did not reproduce the retained mode prefix used by this basis set.");
    end
    generatingMode = auxiliaryBasis.nativeModes(:,nModes+1);
    generatingModeNumber = auxiliaryBasis.modeNumber(nModes+1);
end

zRoots = self.solver.rootsOfNativeMode(generatingMode);
zTolerance = max(100*eps(max(1,max(abs(self.zDomain)))), 1e-10*max(1,diff(self.zDomain)));
zRoots = sort(zRoots(:));
zRoots = zRoots(zRoots > self.zDomain(1)+zTolerance & zRoots < self.zDomain(2)-zTolerance);
if ~isempty(zRoots)
    zRoots = zRoots([true; diff(zRoots) > zTolerance]);
end
z = [self.zDomain(1); zRoots; self.zDomain(2)];
if length(z) < nModes
    error("IMBasisSet:InsufficientQuadraturePoints", ...
        "The auxiliary mode produced %d full-depth points for %d retained modes.", length(z), nModes);
end
gridDesign = IMDiscreteTransformTools.modeRootGridDesign(self,z,nModes,generatingModeNumber);
end
