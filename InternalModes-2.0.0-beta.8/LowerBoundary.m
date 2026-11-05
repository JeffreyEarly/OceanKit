classdef LowerBoundary
    % Enumerate the supported lower boundary conditions.
    %
    % `LowerBoundary` collects the lower-boundary labels used by the
    % `internal-modes` solvers. In the notation of Section 2.3 of Early,
    % Lelong, and Smith (2020), the main user-facing values are
    % `freeSlip`, corresponding to $$G(-D) = 0$$, and `noSlip`,
    % corresponding to $$F(-D) = 0$$.
    %
    % The valid values are:
    %
    % - `LowerBoundary.freeSlip` for $$G(-D) = 0$$
    % - `LowerBoundary.noSlip` for $$F(-D) = 0$$
    % - `LowerBoundary.mda`, `LowerBoundary.buoyancyAnomaly`, and
    %   `LowerBoundary.custom` for specialized internal workflows
    % - `LowerBoundary.none` when no explicit lower condition should be
    %   imposed
    %
    % ```matlab
    % im = InternalModes(rho, zIn, zOut, latitude);
    % im.lowerBoundary = LowerBoundary.freeSlip;
    % ```
    %
    % - Topic: Configure normalization and boundaries
    % - Declaration: classdef LowerBoundary
    enumeration
        freeSlip, noSlip, mda, buoyancyAnomaly, custom, none
    end    
end
