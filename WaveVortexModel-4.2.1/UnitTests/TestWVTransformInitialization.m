classdef TestWVTransformInitialization < matlab.unittest.TestCase
    properties (TestParameter)
        transform = {'constant','hydrostatic','boussinesq','stratified-qg','barotropic-qg'}
        latitude = struct('southOutside',-86,'southUpperBoundary',-85,'southLowerBoundary',-5,'southInsideLowerBoundary',-4,'equator',0,'northInsideLowerBoundary',4,'northLowerBoundary',5,'northUpperBoundary',85,'northOutside',86)
    end

    methods (Test, TestTags = "full")
        function testInitWithLatitude(testCase,transform,latitude)
            constructor = @()TestWVTransformInitialization.transformAtLatitude(transform,latitude);
            if abs(latitude) < 5
                testCase.verifyError(constructor,'MATLAB:validators:mustBeGreaterThanOrEqual')
            elseif abs(latitude) > 85
                testCase.verifyError(constructor,'MATLAB:validators:mustBeLessThanOrEqual')
            else
                testCase.verifyWarningFree(constructor)
            end
        end
    end

    methods (Static, Access=private)
        function wvt = transformAtLatitude(transform,latitude)
            Lxyz = [15e3 15e3 5000];
            Nxyz = [8 8 5];
            N2 = @(z)(5.2e-3)^2*ones(size(z));
            switch transform
                case 'constant'
                    wvt = WVTransformConstantStratification(Lxyz,Nxyz,latitude=latitude,shouldAntialias=false);
                case 'hydrostatic'
                    wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,latitude=latitude,shouldAntialias=false);
                case 'boussinesq'
                    wvt = WVTransformBoussinesq(Lxyz,Nxyz,N2=N2,latitude=latitude,shouldAntialias=false);
                case 'stratified-qg'
                    wvt = WVTransformStratifiedQG(Lxyz,Nxyz,N2=N2,latitude=latitude,shouldAntialias=false);
                case 'barotropic-qg'
                    wvt = WVTransformBarotropicQG(Lxyz(1:2),Nxyz(1:2),latitude=latitude,shouldAntialias=false);
                otherwise
                    error("TestWVTransformInitialization:UnknownTransform","Unknown transform parameter %s.",transform)
            end
        end
    end
end
