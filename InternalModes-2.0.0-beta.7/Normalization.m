classdef Normalization
    % Enumerate internal-mode normalization conventions.
    %
    % `Normalization` names standard scaling conventions installed by
    % internal-mode basis sets. Each enum value resolves to a basis-set
    % normalization rule; `IMInternalModesBasis.normalizationFactors`
    % evaluates the selected rule and returns a per-mode factor $$s_j$$.
    % Evaluated internal-mode variables divide both $$F_j$$ and $$G_j$$ by
    % that same factor:
    % $$V_j^{\mathrm{out}}(z)=V_j^{\mathrm{raw}}(z)/s_j.$$
    % Generic canonical basis sets use string rule names directly, while
    % internal-mode basis sets may use either strings or these enum values.
    %
    % The internal-mode values are:
    %
    % - `Normalization.unity`: unit norm magnitude under the active internal-mode inner product; signed families retain a diagonal signature of `-1` or `+1`
    % - `Normalization.kConstant`: unit `G` inner-product norm for fixed-$$K$$ wave modes
    % - `Normalization.uMax`: scale by $$\max_z |F_j(z)|$$
    % - `Normalization.wMax`: scale by $$\max_z |G_j(z)|$$
    % - `Normalization.surfacePressure`: scale by the raw surface value of $$F$$
    % - `Normalization.geostrophic`: hydrostatic geostrophic normalization
    % - `Normalization.depth`: volume-only unit depth-mean-square `F` normalization, $$s_j=\sqrt{D^{-1}\int F_j^2\,dz}$$, with the same positive factor applied to `F` and `G` and no endpoint terms
    % - `Normalization.generalizedPotentialEnstrophy`: unit depth-mean generalized-potential-enstrophy normalization, including active endpoint terms
    %
    % ```matlab
    % evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-4000 0]);
    % solver = IMSolverSpectral(nEVP=128);
    % basisSet = solver.solveEVP(evp,nModes=4);
    % basisSet.normalization = Normalization.geostrophic;
    % G = basisSet.G(z);
    % ```
    %
    % - Topic: Configure normalization
    % - Declaration: classdef Normalization
    enumeration
        unity, kConstant, uMax, wMax, surfacePressure, geostrophic, depth, generalizedPotentialEnstrophy
    end
end
