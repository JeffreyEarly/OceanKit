classdef (Hidden) LegacyWKBSolver < IMSolverSpectral
    % Retain the legacy N stretch for this coordinate comparison only.
    %
    % The physical EVP keeps its rotation rate and latitude. Only the
    % protected coordinate rotation rate is zeroed, giving dx/dz = N.
    % This helper adds no legacy coordinate option to the production solver.
    %
    % ```matlab
    % solver = IMBandwidthComparison.LegacyWKBSolver(nEVP=64,coordinateKind="wkb");
    % ```
    %
    % - Topic: Compare coordinate choices
    % - Declaration: classdef LegacyWKBSolver < IMSolverSpectral
    % - Developer: true

    methods
        function self = LegacyWKBSolver(options)
            % Forward the spectral solver constructor options.
            %
            % - Topic: Compare coordinate choices
            % - Declaration: solver = LegacyWKBSolver(options)
            % - Parameter options.nEVP: number of EVP coefficients
            % - Parameter options.coordinateKind: native coordinate kind
            % - Returns self: study-only spectral solver
            % - Developer: true
            arguments (Input)
                options.nEVP (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.nEVP,4)} = 64
                options.coordinateKind {mustBeTextScalar,mustBeMember(options.coordinateKind,["z","wkb","density"])} = "z"
            end
            solverOptions = namedargs2cell(options);
            self@IMSolverSpectral(solverOptions{:});
        end
    end

    methods (Access = protected)
        function self = setupCoordinate(self)
            self.coordinateRotationRate = 0;
            self = setupCoordinate@IMSolverSpectral(self);
        end
    end
end
