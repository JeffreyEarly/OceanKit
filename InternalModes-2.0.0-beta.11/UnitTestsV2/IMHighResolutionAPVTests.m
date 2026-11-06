classdef IMHighResolutionAPVTests < matlab.unittest.TestCase
    % Refine the solve without changing the retained physical mode band.
    methods (Test)
        function highResolutionRetainsAnalyticalEigendepths(testCase)
            N0 = 5.2e-3;
            b = 1300;
            D = 4000;
            g = 9.81;
            N2 = @(z) N0^2*exp(2*z/b);
            integratedN2 = N0^2*b*(1-exp(-2*D/b))/2;
            evp = IMInternalModes.geostrophicAPVModes(rotationRate=0, ...
                N2=N2,zDomain=[-D 0],g=g,g0=-integratedN2, ...
                gd=integratedN2,surfaceBoundary="freeSurface");
            exact = IMExponentialStratificationSolution(rotationRate=0, ...
                N0=N0,b=b,zDomain=[-D 0],g=g).internalModes(evp,nModes=84);
            z = linspace(-D,0,1601).';
            exactF = exact.F(z);
            for nEVP = [783 1551]
                solved = IMSolverSpectral(nEVP=nEVP).solveEVP(evp,nModes=84);
                testCase.verifyEqual(solved.modeNumber,exact.modeNumber)
                testCase.verifyLessThan(max(abs(solved.h-exact.h)./abs(exact.h)),1e-6)
                F = solved.F(z);
                orientation = sign(sum(F.*exactF,1));
                % Compare each scientific mode independently; no mixing or
                % fitting is allowed to conceal a changed resolved basis.
                testCase.verifyLessThan(max(vecnorm(F.*orientation-exactF)./vecnorm(exactF)),1e-5)
            end
        end

        function linuxRefinementRetainsFiniteSignedAPVModes(testCase)
            previousThreads = maxNumCompThreads(1);
            testCase.addTeardown(@maxNumCompThreads,previousThreads);
            N2 = @(z) 1e-4*exp(2*z/700);
            zDomain = [-1000 0];
            % Preserve the exact endpoint coefficient from the Linux reproducer.
            integratedN2 = hex2num('3fa0e40bf5376f01');
            background = IMBackground(N2=N2,zDomain=zDomain,g=9.81,latitude=30, ...
                rotationRate=7.2921e-5,rho0=1025,surfaceCondition="freeSurface");
            evp = IMInternalModes.geostrophicAPVModes(background,g0=-integratedN2,gd=integratedN2);
            exact = IMExponentialStratificationSolution(N0=1e-2,b=700,zDomain=zDomain, ...
                g=background.g,latitude=background.latitude,rotationRate=background.rotationRate).internalModes(evp,nModes=3);
            nEVP = 455;
            solver = IMSolverSpectral(nEVP=nEVP).configuredForEVP(evp);
            solved = solver.solveEVP(evp,nModes=101);
            testCase.verifyEqual(solved.modeNumber(1:3),[-1 1 2])
            testCase.verifyNumElements(find(solved.h<0),1)
            testCase.verifyEqual(solved.h(1:3),exact.h,RelTol=1e-9)

            grid = IMSolverSpectral(nEVP=97,coordinateKind="wkb").configuredForEVP(evp);
            [z,weights] = grid.nativeDifferentiationRule(zDomain);
            weights = weights*(background.D/sum(weights));
            [transform,assessment] = solved.discreteTransform(z=z,weights=weights,variables=["F","G"],gramTolerance=1e-2);
            testCase.assertGreaterThanOrEqual(assessment.retainedModeCount,3)
            testCase.verifyEqual(transform.modeNumber(1:3),[-1 1 2])
            testCase.verifyLessThan(max(assessment.prefixDiagnostics.gramError(1:3)),1e-9)

            [A,B,map] = evp.prepareSolvePencil(solver,[],[]);
            testCase.assertTrue(isfield(map,"reduceEquilibratedPencil"))
            physicalScale = map.nativeColumnScale;
            degreeScale = repmat(max(1,0:nEVP-1).^(-2),1,2);
            rowScale = max(max(abs(A.*physicalScale),[],2),max(abs(B.*physicalScale),[],2));
            rowScale(rowScale==0) = 1;
            [finiteA,finiteB] = map.reduceEquilibratedPencil( ...
                (A./rowScale).*(physicalScale.*degreeScale),(B./rowScale).*(physicalScale.*degreeScale));
            testCase.verifySize(finiteA,[nEVP-1 nEVP-1])
            testCase.verifySize(finiteB,[nEVP-1 nEVP-1])
            testCase.verifyEqual(rank(finiteB),nEVP-1)
        end
    end
end
