function [varargout] = variableAtPositionWithName(self,x,y,z,variableNames,options)
% Access dynamical variables at arbitrary positions in the domain.
%
% Computes (or retrieves from cache) any known state variables and computes
% their values at the requested positions (x,y,z). For two-dimensional
% variables, interpolation uses only (x,y), and z may be empty.
%
% The interpolation method may be `linear` or `spline`. Horizontal
% coordinates are wrapped periodically before interpolation.
%
% - Topic: State variables
% - Declaration: [varargout] = variableAtPositionWithName(x,y,z,variableNames,options)
% - Parameter x: array of x-positions
% - Parameter y: array of y-positions
% - Parameter z: array of z-positions, or empty for two-dimensional variables
% - Parameter variableNames: strings of variable names.
% - Parameter options.interpolationMethod: periodic `linear` or cubic `spline` interpolation; default `linear`
% - Returns varargout: interpolated arrays in the requested order and query shape
arguments
    self WVTransform {mustBeNonempty}
    x (1,:) double
    y (1,:) double
    z (1,:) double
end
arguments (Repeating)
    variableNames char
end
arguments
    options.interpolationMethod char {mustBeMember(options.interpolationMethod,["linear","spline"])} = "linear"
end

varargout = cell(size(variableNames));
[varargout{:}] = self.variableWithName(variableNames{:});
[varargout{:}] = self.interpolatedFieldAtPosition(x,y,z,options.interpolationMethod,varargout{:});
end
