classdef IMInstrumentedGeostrophicZeroAPVSolver < IMSolverSpectral
    properties (SetAccess = private)
        solveCounter
    end

    methods
        function self = IMInstrumentedGeostrophicZeroAPVSolver(options)
            arguments
                options.nEVP (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.nEVP,4)} = 64
                options.solveCounter (1,1) IMGeostrophicZeroAPVTestCounter
            end
            self@IMSolverSpectral(nEVP=options.nEVP);
            self.solveCounter = options.solveCounter;
        end
    end

    methods (Access = protected)
        function values = solveBoundaryValueSystems(self,matrix,rightHandSides)
            self.solveCounter.record(size(rightHandSides,2));
            matrixFactorization = decomposition(matrix);
            values = matrixFactorization\rightHandSides;
        end
    end
end
