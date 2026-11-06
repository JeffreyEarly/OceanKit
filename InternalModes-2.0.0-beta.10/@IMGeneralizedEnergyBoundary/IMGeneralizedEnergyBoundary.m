classdef IMGeneralizedEnergyBoundary
    % Define the signed endpoint coefficients of a generalized-energy basis.
    %
    % Finite `g0` and `gd` values, including zero, activate the surface and
    % bottom coordinates respectively. Positive infinity removes that
    % endpoint coordinate. These basis choices are independent of the
    % physical surface convention stored by `IMBackground`.
    %
    % ```matlab
    % boundary = IMGeneralizedEnergyBoundary(g0=-0.02,gd=Inf);
    % ```
    %
    % - Topic: Create energy boundary coefficients
    % - Topic: Inspect energy boundary coefficients
    % - Declaration: classdef IMGeneralizedEnergyBoundary

    properties (SetAccess = private)
        % Signed surface coefficient in meters per second squared.
        %
        % - Topic: Inspect energy boundary coefficients
        g0

        % Signed bottom coefficient in meters per second squared.
        %
        % - Topic: Inspect energy boundary coefficients
        gd
    end

    properties (Dependent, SetAccess = private)
        % Finite-coefficient endpoints, ordered surface before bottom.
        %
        % - Topic: Inspect energy boundary coefficients
        activeEndpoints
    end

    methods
        function self = IMGeneralizedEnergyBoundary(options)
            % Create explicitly declared generalized-energy coefficients.
            %
            % - Topic: Create energy boundary coefficients
            % - Declaration: boundary = IMGeneralizedEnergyBoundary(options)
            % - Parameter options.g0: finite signed surface coefficient or positive infinity
            % - Parameter options.gd: finite signed bottom coefficient or positive infinity
            % - Returns self: read-only endpoint coefficients
            arguments
                options.g0 (1,1) double {mustBeReal}
                options.gd (1,1) double {mustBeReal}
            end
            if isnan(options.g0) || options.g0 == -Inf || isnan(options.gd) || options.gd == -Inf
                error("IMGeneralizedEnergyBoundary:InvalidCoefficient","g0 and gd must be finite signed values or positive infinity.");
            end
            self.g0 = options.g0;
            self.gd = options.gd;
        end

        function endpoints = get.activeEndpoints(self)
            endpoints = ["surface", "bottom"];
            endpoints = endpoints(isfinite([self.g0, self.gd]));
        end
    end
end
