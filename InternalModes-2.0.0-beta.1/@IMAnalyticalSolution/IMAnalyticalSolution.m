classdef (Abstract) IMAnalyticalSolution
    % Describe common state for a stratification family with exact solutions.
    %
    % `IMAnalyticalSolution` owns only physical state shared by analytical
    % stratification families. Concrete classes advertise exact solution
    % operations by implementing those construction methods directly.
    %
    % - Topic: Create analytical solutions
    % - Topic: Inspect analytical solutions
    % - Declaration: classdef (Abstract) IMAnalyticalSolution

    properties (SetAccess = protected)
        % Physical vertical domain.
        %
        % - Topic: Inspect analytical solutions
        zDomain

        % Coriolis parameter in radians per second.
        %
        % - Topic: Inspect analytical solutions
        f0

        % Gravitational acceleration in meters per second squared.
        %
        % - Topic: Inspect analytical solutions
        g
    end

    methods
        function self = IMAnalyticalSolution(options)
            % Initialize common analytical-stratification state.
            %
            % - Topic: Create analytical solutions
            % - Declaration: solution = IMAnalyticalSolution(options)
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Returns solution: analytical solution family
            arguments
                options.zDomain (1,2) double {mustBeReal, mustBeFinite} = [-1 0]
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
            end

            self.zDomain = sort(options.zDomain);
            self.f0 = options.f0;
            self.g = options.g;
        end

        function summarize(self)
            % Print common analytical-stratification state.
            %
            % Concrete families extend this summary with their parameters and
            % exact construction methods.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: summarize(solution)
            arguments
                self IMAnalyticalSolution
            end

            fprintf("%s\n",class(self));
            fprintf("  zDomain: [%g, %g]\n",self.zDomain(1),self.zDomain(2));
            fprintf("  f0: %g\n",self.f0);
            fprintf("  g: %g\n",self.g);
        end
    end

    methods (Abstract)
        % Evaluate buoyancy frequency squared $$N^2(z)$$.
        %
        % - Topic: Inspect analytical solutions
        % - Declaration: values = N2(solution,z)
        % - Parameter z: physical coordinate
        % - Returns values: buoyancy frequency squared
        values = N2(self,z)
    end
end
