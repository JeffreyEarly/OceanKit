classdef (Hidden) SmoothedBandwidthSolver < IMSolverSpectral
    % Compare a Helmholtz-smoothed bandwidth coordinate in this study only.
    %
    % The physical-depth filter solves $$Q-L_s^2 Q_{zz}=B$$ with homogeneous
    % Neumann conditions at both endpoints, then uses $$dx/dz=\sqrt{Q}$$.
    % Here $$B=\sqrt{(N^2+\tilde f^2-f_0^2)^2+4f_0^2\tilde f^2}$$.
    % Smoothing changes only the numerical coordinate; the physical EVP
    % retains its original stratification and rotation. A positive source
    % somewhere gives a positive coordinate through a neutral layer.
    % This monitor is an experimental numerical choice, not a WKB theorem.
    %
    % ```matlab
    % solver = IMBandwidthComparison.SmoothedBandwidthSolver(nEVP=128,smoothingLength=150);
    % ```
    %
    % - Topic: Compare coordinate choices
    % - Declaration: classdef SmoothedBandwidthSolver < IMSolverSpectral
    % - Developer: true

    properties (SetAccess = private)
        % Physical-depth smoothing length $$L_s$$ in metres.
        %
        % - Topic: Compare coordinate choices
        % - Developer: true
        smoothingLength

        % Filtered squared monitor $$Q$$ as a physical-depth Chebfun.
        %
        % Rebuilt whenever the solver is configured for a new physical
        % problem. This diagnostic permits independent checks of the
        % Helmholtz equation and its endpoint conditions.
        %
        % - Topic: Compare coordinate choices
        % - Developer: true
        smoothedBandwidthSquared = []
    end

    methods
        function self = SmoothedBandwidthSolver(options)
            % Create the study-only smoothed coordinate solver.
            %
            % - Topic: Compare coordinate choices
            % - Declaration: solver = SmoothedBandwidthSolver(options)
            % - Parameter options.nEVP: number of EVP coefficients
            % - Parameter options.smoothingLength: positive finite smoothing length in metres
            % - Returns self: study-only spectral solver
            % - Developer: true
            arguments (Input)
                options.nEVP (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.nEVP,4)} = 64
                options.smoothingLength (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 147.77523446738553
            end
            self@IMSolverSpectral(nEVP=options.nEVP,coordinateKind="wkb");
            self.smoothingLength = options.smoothingLength;
        end
    end

    methods (Access = protected)
        function self = setupCoordinate(self)
            physicalN2Function = self.N2Function;
            physicalRotationRate = self.coordinateRotationRate;
            physicalLatitude = self.coordinateLatitude;
            f0 = 2*physicalRotationRate*sind(physicalLatitude);
            horizontalCoriolis = 2*physicalRotationRate*cosd(physicalLatitude);

            % Splitting resolves the finite-smoothness neutral-layer joins
            % before the physical-depth differential filter is applied.
            bandwidthSquared = chebfun(@fullBandwidthSquared,self.zDomain,'splitting','on');
            if max(bandwidthSquared) <= 0
                error("IMSolverSpectral:InvalidCoordinate","The smoothed coordinate requires nonzero bandwidth somewhere in the domain.");
            end
            smoothingLengthSquared = self.smoothingLength^2;
            filter = chebop(self.zDomain);
            filter.op = @(z,Q) Q-smoothingLengthSquared*diff(Q,2);
            filter.lbc = @(Q) diff(Q);
            filter.rbc = @(Q) diff(Q);
            Q = filter\bandwidthSquared;
            if ~isreal(Q) || min(Q) <= 0
                error("IMSolverSpectral:InvalidCoordinate","The filtered squared coordinate derivative must remain finite, real, and positive.");
            end
            self.smoothedBandwidthSquared = Q;

            % Reuse the parent's single Chebfun map, inverse, and derivative
            % representation, then restore the physical problem metadata.
            self.N2Function = @(z) feval(Q,z);
            self.coordinateRotationRate = 0;
            self = setupCoordinate@IMSolverSpectral(self);
            self.N2Function = physicalN2Function;
            self.coordinateRotationRate = physicalRotationRate;
            self.coordinateLatitude = physicalLatitude;

            function B = fullBandwidthSquared(z)
                N2 = physicalN2Function(z);
                if ~isreal(N2) || any(~isfinite(N2(:))) || any(N2(:) < 0)
                    error("IMSolverSpectral:InvalidCoordinate","The smoothed coordinate requires finite real nonnegative physical N2.");
                end
                B = hypot(N2+horizontalCoriolis^2-f0^2,2*f0*horizontalCoriolis);
            end
        end
    end
end
