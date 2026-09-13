classdef TestNoMotionProfileSolver < matlab.unittest.TestCase
    methods (Test, TestTags="full")
        function discreteMomentsUseTheActualVolumeWeights(testCase)
            % End planes contain equal volumes of densities 1000 and 1004.
            % The middle plane contains density 1002. Their total masses
            % give moments .35 + .3*(1/2)^k independently of array order.
            rho = reshape([1000 1004 1002 1002 1004 1000],2,1,3);
            actual = WVNoMotionProfileOperation.moments_from_rho_tot(rho,1000,1004,[.2;.3;.5],1);
            testCase.verifyEqual(actual,[.5;.425;.3875],AbsTol=2e-16);
        end

        function changedStableRestProfileNeedsNoOptimization(testCase)
            z = linspace(-1,0,17).';
            reference = 1025-z;
            changed = 1025-2*z-.3*z.^2;
            weights = ones(size(z))/16;
            weights([1 end]) = weights([1 end])/2;
            density = repmat(reshape(changed,1,1,[]),4,3,1);
            [actual,exitflag,output] = WVNoMotionProfileOperation.find_rho_nm(weights,1,density,reference);
            testCase.verifyEqual(actual,changed);
            testCase.verifyGreaterThan(exitflag,0);
            testCase.verifyEqual(output.iterations,0);
            profile = WVNoMotionProfile(z,actual);
            testCase.verifyEqual(profile.inverse(changed),z,AbsTol=2e-15);
            testCase.verifyEqual(profile.availablePotentialEnergy(z,z,9.81,1025),zeros(size(z)));
        end

        function profileFitRejectsInconsistentDimensions(testCase)
            testCase.verifyError(@()WVNoMotionProfileOperation.find_rho_nm([.25;.5;.25],1,ones(2,2,4),[3;2;1]),'WVNoMotionProfileOperation:InvalidDimensions');
        end

        function dampedSolverQualifiesBothCapturedMomentFixtures(testCase)
            root = fileparts(mfilename('fullpath'));
            for day = [3000,3250]
                fixture = jsondecode(fileread(fullfile(root,'fixtures','density-run18-day'+string(day)+'.json')));
                rho0 = fixture.rho_nm0(end);
                rhoD = fixture.rho_nm0(1);
                u0 = WVNoMotionProfileOperation.rho_to_u(fixture.rho_nm0,rho0,rhoD);
                [u,exitflag,output] = WVNoMotionProfileOperation.solveMoments(u0,flip(fixture.z_int/fixture.Lz),fixture.moments);
                profile = WVNoMotionProfileOperation.u_to_rho(u,rho0,rhoD);
                testCase.verifyEqual(exitflag,1);
                testCase.verifyLessThan(output.maximumResidual,1e-9);
                testCase.verifyLessThanOrEqual(output.gradientNorm,1e-12);
                testCase.verifyLessThan(max(diff(profile)),0);
                testCase.verifyEqual(profile([1 end]),[rhoD;rho0]);
                testCase.verifyGreaterThan(output.acceptedSteps,0);
                testCase.verifyLessThan(output.finalCost,output.initialCost/1e8);
            end
        end

        function exhaustedSolverBudgetsReturnFailureStatus(testCase)
            weights = [.25;.5;.25];
            u0 = WVNoMotionProfileOperation.rho_to_u([3;2;1],1,3);
            target = [.55;.4;.33];
            [u,exitflag,output] = WVNoMotionProfileOperation.solveMoments(u0,weights,target,maximumIterations=0);
            testCase.verifyEqual(u,u0);
            testCase.verifyEqual(exitflag,0);
            testCase.verifyEqual(output.reason,"iteration-limit");
            testCase.verifyEqual(output.iterations,0);
            testCase.verifyEqual(output.evaluations,1);
            [u,exitflag,output] = WVNoMotionProfileOperation.solveMoments(u0,weights,target,maximumEvaluations=1);
            testCase.verifyEqual(u,u0);
            testCase.verifyEqual(exitflag,0);
            testCase.verifyEqual(output.reason,"evaluation-limit");
            testCase.verifyEqual(output.evaluations,1);
            testCase.verifyError(@()WVNoMotionProfileOperation.solveMoments(u0,[.5;.5],target),'WVNoMotionProfileOperation:InvalidDimensions');
        end

        function rejectedFitIsReportedAndNeverCached(testCase)
            wvt = WVTransformConstantStratification([4000 4000 1000],[8 8 5],shouldAntialias=false);
            operation = WVNoMotionProfileOperation(solver="dampedLeastSquares");
            wvt.addOperation(operation,shouldOverwriteExisting=true,shouldSuppressWarning=true);
            % A single outlying parcel has less volume than the endpoint
            % quadrature weight. The fixed-node profile cannot represent its
            % distribution: a small step must not certify the fit.
            density = repmat(reshape(wvt.rho_nm0,1,1,[]),wvt.Nx,wvt.Ny,1);
            density(1,1,3) = max(density,[],"all")+10;
            wvt.addToVariableCache('rho_total',density);
            testCase.verifyError(@()wvt.performOperationWithName('rho_nm'),'WVNoMotionProfileOperation:UnqualifiedFit');
            testCase.verifyFalse(isKey(wvt.variableCache,'rho_nm'));
            testCase.verifyGreaterThan(operation.lastSolverOutput.maximumResidual,1e-8);
            testCase.verifyEqual(operation.lastSolverOutput.solver,"dampedLeastSquares");
        end

        function changedExtremaAreRecoveredAfterEqualVolumeRearrangement(testCase)
            z = linspace(-1,0,17).';
            reference = 1025-z;
            changed = 1026-2*z;
            weights = ones(size(z))/16;
            weights([1 end]) = weights([1 end])/2;
            density = repmat(reshape(changed,1,1,[]),4,3,1);
            density(1,1,[3 12]) = density(1,1,[12 3]);
            [actual,exitflag,output] = WVNoMotionProfileOperation.find_rho_nm(weights,1,density,reference);
            testCase.verifyGreaterThan(exitflag,0);
            testCase.verifyLessThan(output.maximumResidual,1e-12);
            testCase.verifyEqual(actual,changed,AbsTol=2e-12);
            testCase.verifyEqual(actual([1 end]),[max(density,[],"all");min(density,[],"all")]);
        end

        function constantDensityHasNoUniqueMaterialHeight(testCase)
            testCase.verifyError(@()WVNoMotionProfileOperation.find_rho_nm([.25;.5;.25],1,ones(2,2,3),[3;2;1]),'WVNoMotionProfileOperation:NonInvertibleDistribution');
        end

        function qualifiedParcelRearrangementReportsItsSolver(testCase)
            % Interior levels of the constant-stratification grid carry
            % equal quadrature volumes. Swapping those parcels in just one
            % column preserves the complete discrete distribution.
            wvt = WVTransformConstantStratification([4000 4000 1000],[8 8 5],shouldAntialias=false);
            operation = WVNoMotionProfileOperation(solver="dampedLeastSquares");
            wvt.addOperation(operation,shouldOverwriteExisting=true,shouldSuppressWarning=true);
            density = repmat(reshape(wvt.rho_nm0,1,1,[]),wvt.Nx,wvt.Ny,1);
            testCase.assertEqual(wvt.z_int(2),wvt.z_int(4));
            density(1,1,[2 4]) = density(1,1,[4 2]);
            wvt.addToVariableCache('rho_total',density);
            actual = wvt.rho_nm;
            testCase.verifyEqual(actual,wvt.rho_nm0,AbsTol=2e-12);
            testCase.verifyEqual(operation.lastSolverOutput.exitflag,1);
            testCase.verifyLessThan(operation.lastSolverOutput.maximumResidual,1e-12);
            testCase.verifyTrue(isKey(wvt.variableCache,'rho_nm'));
            profile = WVNoMotionProfile(wvt.z,actual);
            expectedHeight = wvt.Z;
            expectedHeight(1,1,[2 4]) = expectedHeight(1,1,[4 2]);
            actualHeight = profile.inverse(density);
            testCase.verifyEqual(actualHeight,expectedHeight,AbsTol=1e-8);
            expectedAPE = wvt.N2(1)*(wvt.Z-expectedHeight).^2/2;
            testCase.verifyEqual(profile.availablePotentialEnergy(wvt.Z,actualHeight,wvt.g,wvt.rho0),expectedAPE,AbsTol=1e-10);
        end
    end
end
