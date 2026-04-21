classdef Normalization
    % Enumerate the supported modal normalization conventions.
    %
    % `Normalization` collects the normalization choices used throughout
    % the `internal-modes` hierarchy. Following Section 2.4 of Early,
    % Lelong, and Smith (2020), the main wave-mode norms are the
    % `kConstant` norm based on the fixed-$$K$$ orthogonality relation and
    % the `omegaConstant` norm based on the fixed-$$\omega$$ relation.
    %
    % The valid values are:
    %
    % - `Normalization.kConstant` for the manuscript's $$K$$-constant norm
    % - `Normalization.omegaConstant` for the manuscript's
    %   $$\omega$$-constant norm
    % - `Normalization.uMax` to scale each mode so $$\max F_j = 1$$
    % - `Normalization.wMax` to scale each mode so $$\max G_j = 1$$
    % - `Normalization.surfacePressure` to scale by the surface value of
    %   $$F$$
    % - `Normalization.geostrophic` for the near-geostrophic interior mode
    %   normalization used by some helper workflows
    %
    % ```matlab
    % im = InternalModes(rho, zIn, zOut, latitude);
    % im.normalization = Normalization.kConstant;
    % ```
    %
    % - Topic: Configure normalization and boundaries
    % - Declaration: classdef Normalization
    enumeration
        kConstant, omegaConstant, uMax, wMax, surfacePressure, geostrophic
    end    
end
