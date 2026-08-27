classdef TestRandomFlow < matlab.unittest.TestCase
    properties
        wvt
        solutionGroup
    end

    properties (ClassSetupParameter)
        Lxyz = struct('Lxyz',[1000, 500, 500]);
        Nxyz = struct('Nx16Ny8Nz9',[16 8 9]);
        transform = {'constant'};
    end

    methods (TestClassSetup)
        function classSetup(testCase,Lxyz,Nxyz,transform)
            switch transform
                case 'constant'
                    testCase.wvt = WVTransformConstantStratification(Lxyz, Nxyz, latitude=33, isHydrostatic=0);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
            end
        end
    end

    methods (TestParameterDefinition,Static)
        function [flowComponent] = initializeProperty(Lxyz,Nxyz,transform)
            % If you want to dynamically adjust the test parameters, you
            % have to do it here.
            switch transform
                case 'constant'
                    tmpwvt = WVTransformConstantStratification(Lxyz, Nxyz, latitude=33, isHydrostatic=0);
                case 'hydrostatic'
                    tmpwvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
                case 'boussinesq'
                    tmpwvt = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
            end
            
            flowComponent = cellstr(tmpwvt.flowComponentNames);
        end
    end

    properties (TestParameter)
        flowComponent
    end

    methods (Test, TestTags = "full")
        function testSolution(self,flowComponent)
            seedRandomNumberGenerator(self,17431);
            self.wvt.removeAll();
            self.wvt.initWithRandomFlow(flowComponent,uvMax=0.1);
            self.verifyEqual(self.wvt.totalEnergy,self.wvt.totalEnergySpatiallyIntegrated, "AbsTol",1e-7,"RelTol",1e-7);
        end

    end

end
