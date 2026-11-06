classdef IMScalarTransformCountingSolver < IMSolverSpectral
    properties (SetAccess = private)
        testCounter (1,1) IMScalarTransformTestCounter
    end

    methods
        function self = IMScalarTransformCountingSolver(counter)
            arguments
                counter (1,1) IMScalarTransformTestCounter
            end
            self@IMSolverSpectral(nEVP=64,coordinateKind="z");
            self.testCounter = counter;
        end
    end

    methods (Access = protected)
        function [basis,costs,selection] = solveAssembledEVP(self,varargin)
            self.testCounter.solveCount = self.testCounter.solveCount+1;
            [basis,costs,selection] = solveAssembledEVP@IMSolverSpectral(self,varargin{:});
        end
    end
end
