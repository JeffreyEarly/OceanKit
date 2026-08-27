classdef TestVerticalCalculus < matlab.unittest.TestCase
    properties
        wvt
        isConstant
    end

    properties (ClassSetupParameter)
        transformType = struct( ...
            constantHydrostatic="constantHydrostatic", ...
            constantNonhydrostatic="constantNonhydrostatic", ...
            hydrostatic="hydrostatic", ...
            boussinesq="boussinesq", ...
            stratifiedQG="stratifiedQG")
    end

    methods (TestClassSetup)
        function createTransform(testCase,transformType)
            Lxyz = [6e3 4e3 1e3];
            Nxyz = [6 6 9];
            N2 = @(z) 2e-5*exp(z/4000);
            switch transformType
                case "constantHydrostatic"
                    testCase.wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=true,shouldAntialias=false);
                    testCase.isConstant = true;
                case "constantNonhydrostatic"
                    testCase.wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=false);
                    testCase.isConstant = true;
                case "hydrostatic"
                    testCase.wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false);
                    testCase.isConstant = false;
                case "boussinesq"
                    testCase.wvt = WVTransformBoussinesq(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false);
                    testCase.isConstant = false;
                case "stratifiedQG"
                    testCase.wvt = WVTransformStratifiedQG(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false);
                    testCase.isConstant = false;
            end
        end
    end

    methods (Test, TestTags = "full")
        function constantModesDifferentiateAnalytically(testCase)
            if ~testCase.isConstant
                return
            end

            wvt = testCase.wvt;
            [X,Y,Z] = wvt.xyzGrid;
            horizontal = cos(2*pi*X/wvt.Lx).*cos(2*pi*Y/wvt.Ly);
            mode = 2;
            m = mode*pi/wvt.Lz;
            theta = m*(Z+wvt.Lz);
            F = horizontal.*cos(theta);
            G = horizontal.*sin(theta);
            expectedF = { ...
                -m*horizontal.*sin(theta), ...
                -(m^2)*horizontal.*cos(theta), ...
                (m^3)*horizontal.*sin(theta), ...
                (m^4)*horizontal.*cos(theta)};
            expectedG = { ...
                m*horizontal.*cos(theta), ...
                -(m^2)*horizontal.*sin(theta), ...
                -(m^3)*horizontal.*cos(theta), ...
                (m^4)*horizontal.*sin(theta)};

            for n = 1:4
                testCase.verifyRelative(wvt.diffZF(F,n=n),expectedF{n},2e-11)
                testCase.verifyRelative(wvt.diffZG(G,n=n),expectedG{n},2e-11)
            end
        end

        function variableModesSatisfyDefiningIdentities(testCase)
            if testCase.isConstant
                return
            end

            wvt = testCase.wvt;
            [X,Y,~] = wvt.xyzGrid;
            horizontal = cos(2*pi*X/wvt.Lx).*cos(2*pi*Y/wvt.Ly);
            modeIndex = 2;
            F = horizontal.*reshape(wvt.FinvMatrix(:,modeIndex),1,1,wvt.Nz);
            G = horizontal.*reshape(wvt.GinvMatrix(:,modeIndex),1,1,wvt.Nz);
            expectedDzF = -reshape(wvt.N2,1,1,wvt.Nz).*G/wvt.g;
            expectedDzG = F/wvt.h_0(modeIndex);

            testCase.verifyRelative(wvt.diffZF(F),expectedDzF,2e-11)
            testCase.verifyRelative(wvt.diffZG(G),expectedDzG,2e-11)
        end

        function directAndRepeatedDerivativesAgree(testCase)
            wvt = testCase.wvt;
            [F,G] = testCase.modeFields;

            repeatedF = wvt.diffZF(F);
            repeatedG = wvt.diffZG(G);
            testCase.verifyRelative(wvt.diffZF(F,n=1),repeatedF,2e-11)
            testCase.verifyRelative(wvt.diffZG(G,n=1),repeatedG,2e-11)
            for n = 2:4
                if mod(n,2) == 0
                    repeatedF = wvt.diffZG(repeatedF);
                    repeatedG = wvt.diffZF(repeatedG);
                else
                    repeatedF = wvt.diffZF(repeatedF);
                    repeatedG = wvt.diffZG(repeatedG);
                end
                testCase.verifyRelative(wvt.diffZF(F,n=n),repeatedF,2e-11)
                testCase.verifyRelative(wvt.diffZG(G,n=n),repeatedG,2e-11)
            end
        end

        function firstAntiderivativesSatisfyContracts(testCase)
            wvt = testCase.wvt;
            [F,G] = testCase.modeFields;

            intF = wvt.intZF(F);
            intG = wvt.intZG(G);
            testCase.verifyRelative(wvt.diffZG(intF),F,2e-11)
            testCase.verifyRelative(wvt.diffZF(intG),G,2e-11)
            testCase.verifyLessThanOrEqual(max(abs(intF(:,:,[1 end])),[],"all"),2e-11*max(norm(intF(:)),1))
            testCase.verifyLessThanOrEqual(max(abs(intG(:,:,1)),[],"all"),2e-11*max(norm(intG(:)),1))

            barotropicF = wvt.FinvMatrix(:,1);
            testCase.verifyLessThanOrEqual(norm(wvt.intZF(barotropicF)),2e-11*max(norm(barotropicF),1))
        end

        function integrationPreservesSupportedLayouts(testCase)
            wvt = testCase.wvt;
            [F,G] = testCase.modeFields;
            FMatrix = reshape(permute(F,[3 1 2]),wvt.Nz,[]);
            GMatrix = reshape(permute(G,[3 1 2]),wvt.Nz,[]);

            intFGrid = wvt.intZF(F);
            intGGrid = wvt.intZG(G);
            intFMatrix = wvt.intZF(FMatrix);
            intGMatrix = wvt.intZG(GMatrix);
            testCase.verifySize(intFGrid,[wvt.Nx wvt.Ny wvt.Nz])
            testCase.verifySize(intGGrid,[wvt.Nx wvt.Ny wvt.Nz])
            testCase.verifySize(intFMatrix,size(FMatrix))
            testCase.verifySize(intGMatrix,size(GMatrix))
            testCase.verifyEqual(intFMatrix,reshape(permute(intFGrid,[3 1 2]),wvt.Nz,[]),RelTol=2e-13,AbsTol=2e-13)
            testCase.verifyEqual(intGMatrix,reshape(permute(intGGrid,[3 1 2]),wvt.Nz,[]),RelTol=2e-13,AbsTol=2e-13)

            FColumn = wvt.FinvMatrix(:,2);
            GColumn = wvt.GinvMatrix(:,2);
            testCase.verifySize(wvt.intZF(FColumn),[wvt.Nz 1])
            testCase.verifySize(wvt.intZG(GColumn),[wvt.Nz 1])
        end

        function invalidOrdersAndLayoutsUseBuiltinValidation(testCase)
            wvt = testCase.wvt;
            grid = zeros(wvt.Nx,wvt.Ny,wvt.Nz);
            testCase.verifyError(@()wvt.diffZF(grid,n=0),"MATLAB:validators:mustBeMember")
            testCase.verifyError(@()wvt.diffZG(grid,n=5),"MATLAB:validators:mustBeMember")
            testCase.verifyError(@()wvt.intZF(grid,n=2),"MATLAB:validators:mustBeMember")
            testCase.verifyError(@()wvt.intZG(grid,n=2),"MATLAB:validators:mustBeMember")
            testCase.verifyError(@()wvt.diffZF(zeros(wvt.Nx,wvt.Ny)),"MATLAB:validators:mustBeMember")
            testCase.verifyError(@()wvt.diffZG(zeros(wvt.Nx+1,wvt.Ny,wvt.Nz)),"MATLAB:validators:mustBeMember")
            testCase.verifyError(@()wvt.intZF(zeros(wvt.Nz-1,2)),"MATLAB:validators:mustBeMember")
            testCase.verifyError(@()wvt.intZG(zeros(wvt.Nz,2,1,2)),"MATLAB:validators:mustBeMatrix")
        end
    end

    methods (Access = private)
        function [F,G] = modeFields(testCase)
            wvt = testCase.wvt;
            [X,Y,~] = wvt.xyzGrid;
            horizontal = cos(2*pi*X/wvt.Lx).*cos(2*pi*Y/wvt.Ly);
            FProfile = wvt.FinvMatrix(:,2) + 0.3*wvt.FinvMatrix(:,3);
            GProfile = wvt.GinvMatrix(:,2) - 0.2*wvt.GinvMatrix(:,3);
            F = horizontal.*reshape(FProfile,1,1,wvt.Nz);
            G = horizontal.*reshape(GProfile,1,1,wvt.Nz);
        end
    end

    methods (Access = private)
        function verifyRelative(testCase,actual,expected,tolerance)
            errorValue = norm(actual(:)-expected(:))/max(norm(expected(:)),eps);
            testCase.verifyLessThanOrEqual(errorValue,tolerance)
        end
    end
end
