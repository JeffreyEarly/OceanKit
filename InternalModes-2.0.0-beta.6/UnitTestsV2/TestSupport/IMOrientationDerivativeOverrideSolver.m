classdef IMOrientationDerivativeOverrideSolver < IMSolverSpectral
    methods
        function self = IMOrientationDerivativeOverrideSolver(options)
            arguments
                options.nEVP (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.nEVP,4)} = 33
            end
            self@IMSolverSpectral(nEVP=options.nEVP);
        end

        function values = differentiateGridValues(self,values,derivativeOrder)
            values = differentiateGridValues@IMSolverSpectral(self,values,derivativeOrder);
            if derivativeOrder == 1
                values = -values;
            end
        end
    end
end
