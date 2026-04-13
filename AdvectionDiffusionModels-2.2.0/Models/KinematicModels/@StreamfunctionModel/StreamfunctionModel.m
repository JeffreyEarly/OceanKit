classdef StreamfunctionModel < KinematicModel
    % Base class for incompressible two-dimensional streamfunction models.
    %
    % Streamfunction models define a scalar field $$\psi(t,x,y)$$ such that
    %
    % $$ u(t,x,y) = -\frac{\partial \psi}{\partial y}, \qquad v(t,x,y) = \frac{\partial \psi}{\partial x}. $$
    %
    % Subclasses implement `psi(self, t, x, y)` and inherit the common
    % streamfunction plotting helper.
    %
    % - Topic: Evaluate streamfunctions
    % - Topic: Plot streamfunctions
    % - Declaration: classdef StreamfunctionModel < KinematicModel

    methods
        function self = StreamfunctionModel()
            % Create a streamfunction-model base instance for subclass setup.
            %
            % `StreamfunctionModel` is abstract, so this constructor is
            % used indirectly by subclasses that implement
            % `psi(self,t,x,y)`.
            %
            % - Topic: Evaluate streamfunctions
            % - Declaration: self = StreamfunctionModel()
            % - Returns self: `StreamfunctionModel` instance
            % - Developer: true
        end
    end

    methods (Abstract)
        % Evaluate the streamfunction.
        %
        % Output shape matches the input `x` and `y` arrays.
        %
        % - Topic: Evaluate streamfunctions
        % - Declaration: psiValue = psi(self,t,x,y)
        % - Parameter t: scalar evaluation time in seconds
        % - Parameter x: x-coordinate array in meters
        % - Parameter y: y-coordinate array in meters
        % - Returns psiValue: streamfunction values with the same shape as `x` and `y`
        psiValue = psi(self, t, x, y);
    end
end
