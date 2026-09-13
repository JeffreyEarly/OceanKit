classdef TestDensityDiagnosticReference < matlab.unittest.TestCase
    methods (Test, TestTags="full")
        function actualDefaultPipelineConvergesForContinuousRearrangement(testCase)
            % Initialize physical velocity and displacement through the
            % transform, then use its actual default fitter and diagnostics.
            % The independent polar map supplies continuum truth; finite
            % spectral projection error is measured separately.
            root = fileparts(fileparts(mfilename('fullpath')));
            originalPath = path;
            pathCleanup = onCleanup(@()path(originalPath));
            addpath(fullfile(root,'tools','density-diagnostics'));
            report = qualifyContinuousDensityDiagnostics();
            errors = zeros(3,3);
            spectralErrors = zeros(3,1);
            for index = 1:3
                result = report.cases{index};
                testCase.assertEqual(result.status,"returned");
                testCase.verifyTrue(result.actualReferenceDefault);
                testCase.verifyEqual(result.solver.solver,"dampedLeastSquares");
                testCase.verifyGreaterThan(result.solver.exitflag,0);
                testCase.verifyLessThanOrEqual(result.solver.maximumResidual,1e-8);
                testCase.verifyLessThanOrEqual(result.densityClosureMaximum,8*eps(result.densityMagnitude));
                testCase.verifyGreaterThanOrEqual(result.minimumAPE,0);
                testCase.verifyTrue(result.apvAllFinite);
                errors(index,:) = [result.etaRMSError,result.apeRMSError,result.apvRMSError];
                spectralErrors(index) = result.spectralDensityMaximumError;
            end
            testCase.verifyGreaterThan(errors(1,:),[1e-5,1e-8,1e-4]);
            testCase.verifyLessThan(errors(3,:),errors(1,:)./[1000,1000,20]);
            testCase.verifyLessThan(errors(3,1:2),[1e-7,1e-9]);
            testCase.verifyLessThan(report.cases{3}.apvRelativeRMSError,1e-3);
            testCase.verifyLessThan(spectralErrors(3),spectralErrors(1)/1000);
            fprintf('Default density pipeline eta/APE/APV RMS errors by grid: %s\n',mat2str(errors,9));
        end

        function apvConvergesToAnalyticMaterialHeightVorticity(testCase)
            % Isolate APV's displacement-input boundary from the separately
            % qualified default fitter.
            % The reference is (curl(u)+f*e_z).grad(materialHeight)-f.
            errors = zeros(2,3);
            referenceRMS = zeros(2,3);
            for axis = 1:2
                for level = 1:3
                    n = 16*2^(level-1);
                    wvt = WVTransformConstantStratification([2 2 2],[n n n+1],N0=.2,shouldAntialias=false);
                    X = wvt.X;
                    Y = wvt.Y;
                    Z = wvt.Z;
                    u = .03*sin(pi*Y).*cos(pi*Z/2);
                    v = .02*sin(pi*X).*cos(pi*Z/2);
                    wvt.initWithUVEta(u,v,zeros(size(Z)));
                    zetaX = .02*(pi/2)*sin(pi*X).*sin(pi*Z/2);
                    zetaY = -.03*(pi/2)*sin(pi*Y).*sin(pi*Z/2);
                    zetaZ = pi*(.02*cos(pi*X)-.03*cos(pi*Y)).*cos(pi*Z/2);
                    if axis == 1
                        coordinate = X;
                        horizontalVorticity = zetaX;
                    else
                        coordinate = Y;
                        horizontalVorticity = zetaY;
                    end
                    radius = .7;
                    angle = 3;
                    dx = coordinate-1;
                    dz = Z+1;
                    q = max(0,1-(dx.^2+dz.^2)/radius^2);
                    theta = angle*q.^3;
                    thetaH = -6*angle*dx.*q.^2/radius^2;
                    thetaZ = -6*angle*dz.*q.^2/radius^2;
                    horizontalArm = cos(theta).*dx+sin(theta).*dz;
                    gradientH = -sin(theta)-horizontalArm.*thetaH;
                    gradientZ = cos(theta)-horizontalArm.*thetaZ;
                    [~,materialHeight] = DensityDiagnosticReference.inversePolarTwist(coordinate,Z,[1 -1],radius,angle);
                    profile = WVNoMotionProfile(wvt.z,1025-4*wvt.z);
                    recoveredHeight = profile.inverse(1025-4*materialHeight);
                    testCase.verifyEqual(recoveredHeight,materialHeight,AbsTol=6e-14);
                    wvt.addToVariableCache('eta_true',Z-recoveredHeight);
                    actual = wvt.apv;
                    expected = horizontalVorticity.*gradientH+(zetaZ+wvt.f).*gradientZ-wvt.f;
                    errors(axis,level) = sqrt(mean((actual-expected).^2,'all'));
                    referenceRMS(axis,level) = sqrt(mean(expected.^2,'all'));
                    testCase.verifyLessThan(max(abs(wvt.zeta_x-zetaX),[],'all'),1e-12);
                    testCase.verifyLessThan(max(abs(wvt.zeta_y-zetaY),[],'all'),1e-12);
                    testCase.verifyLessThan(max(abs(wvt.zeta_z-zetaZ),[],'all'),1e-12);
                    testCase.verifyLessThan(min(gradientZ,[],'all'),0);
                    testCase.verifyGreaterThan(max(abs(expected-zetaZ),[],'all'),.01);
                end
            end
            testCase.verifyLessThan(errors(:,3),errors(:,1)/20);
            % The compact-support map is C2, so its differentiated field
            % converges algebraically. Require sub-0.1% relative RMS APV
            % error at the finest grid as well as the refinement factor.
            testCase.verifyLessThan(errors(:,3)./referenceRMS(:,3),1e-3);
            fprintf('APV x-z RMS errors: %.9g %.9g %.9g; relative finest %.9g\n',errors(1,:),errors(1,3)/referenceRMS(1,3));
            fprintf('APV y-z RMS errors: %.9g %.9g %.9g; relative finest %.9g\n',errors(2,:),errors(2,3)/referenceRMS(2,3));
        end

        function unequalParcelVolumesDefineTheEmpiricalDistribution(testCase)
            z = [-1;-.5;.8];
            weights = DensityDiagnosticReference.verticalCellWeights(z,[-1 1]);
            testCase.verifyEqual(weights,[.25;.9;.85],AbsTol=2e-15);
            values = [3;1;3];
            distribution = DensityDiagnosticReference.empiricalDistribution(values,weights);
            testCase.verifyEqual(distribution.values,[1;3]);
            testCase.verifyEqual(distribution.masses,[.9;1.1]/2,AbsTol=2e-15);
            testCase.verifyEqual(distribution.cumulativeMass,[.45;1],AbsTol=2e-15);
            testCase.verifyEqual(distribution.totalVolume,2);
            testCase.verifyEqual(sum(distribution.values.*distribution.masses),sum(values.*weights)/sum(weights),AbsTol=2e-15);
            unweighted = DensityDiagnosticReference.empiricalDistribution(values,ones(size(weights)));
            testCase.verifyGreaterThan(abs(distribution.cumulativeMass(1)-unweighted.cumulativeMass(1)),.1);
            permuted = DensityDiagnosticReference.empiricalDistribution(values([3 1 2]),weights([3 1 2]));
            testCase.verifyEqual(permuted,distribution);
        end

        function zeroVolumesAreExcludedAndInvalidWeightsRejected(testCase)
            distribution = DensityDiagnosticReference.empiricalDistribution([1 2 3],[1 0 3]);
            testCase.verifyEqual(distribution.values,[1;3]);
            testCase.verifyEqual(distribution.masses,[.25;.75]);
            testCase.verifyError(@()DensityDiagnosticReference.empiricalDistribution([1 2],[0 0]),'DensityDiagnosticReference:InvalidWeights');
            testCase.verifyError(@()DensityDiagnosticReference.empiricalDistribution([1 2],1),'DensityDiagnosticReference:InvalidWeights');
            testCase.verifyError(@()DensityDiagnosticReference.verticalCellWeights([0;-1],[-1 1]),'DensityDiagnosticReference:InvalidVerticalGrid');
        end

        function interiorTwistHasAnalyticInverseAndFixedExterior(testCase)
            [x,z] = ndgrid(linspace(-1,1,31));
            center = [.13 -.17];
            radius = .7;
            [mappedX,mappedZ] = DensityDiagnosticReference.polarTwist(x,z,center,radius,1.3);
            [materialX,materialHeight] = DensityDiagnosticReference.inversePolarTwist(mappedX,mappedZ,center,radius,1.3);
            testCase.verifyEqual(materialX,x,AbsTol=6e-16);
            testCase.verifyEqual(materialHeight,z,AbsTol=6e-16);
            outside = (x-center(1)).^2+(z-center(2)).^2 >= radius^2;
            testCase.verifyEqual(mappedX(outside),x(outside));
            testCase.verifyEqual(mappedZ(outside),z(outside));
            testCase.verifyGreaterThan(max(abs(mappedZ-z),[],"all"),.1);
            testCase.verifyEqual((mappedX-center(1)).^2+(mappedZ-center(2)).^2, ...
                (x-center(1)).^2+(z-center(2)).^2,AbsTol=1e-15);
        end

        function interiorTwistPreservesTheContinuousJacobian(testCase)
            x = [.1;.25;-.3;.8];
            z = [-.2;.15;-.25;.8];
            h = 1e-5;
            [xp,zp] = DensityDiagnosticReference.polarTwist(x+h,z,[.13 -.17],.7,1.3);
            [xm,zm] = DensityDiagnosticReference.polarTwist(x-h,z,[.13 -.17],.7,1.3);
            dXdx = (xp-xm)/(2*h);
            dZdx = (zp-zm)/(2*h);
            [xp,zp] = DensityDiagnosticReference.polarTwist(x,z+h,[.13 -.17],.7,1.3);
            [xm,zm] = DensityDiagnosticReference.polarTwist(x,z-h,[.13 -.17],.7,1.3);
            determinant = dXdx.*((zp-zm)/(2*h))-dZdx.*((xp-xm)/(2*h));
            testCase.verifyEqual(determinant,ones(size(x)),AbsTol=5e-9);
        end

        function linearStableProfileHasAnalyticDisplacementAndEnergy(testCase)
            [x,z] = ndgrid(linspace(-1,1,25));
            [~,materialHeight] = DensityDiagnosticReference.inversePolarTwist(x,z,[.13 -.17],.7,1.3);
            N2 = .04;
            [eta,ape] = DensityDiagnosticReference.linearReference(z,materialHeight,N2);
            % In a linear profile b=N2*z, the displaced buoyancy anomaly is
            % -N2*eta. Integrating that restoring force gives N2*eta^2/2.
            buoyancyAnomaly = N2*materialHeight-N2*z;
            testCase.verifyEqual(eta,-buoyancyAnomaly/N2,AbsTol=3e-16);
            referenceEnergy = N2*(z.^2/2-materialHeight.*z+materialHeight.^2/2);
            testCase.verifyEqual(ape,referenceEnergy,AbsTol=1e-17);
            testCase.verifyGreaterThan(max(ape,[],"all"),1e-4);
            testCase.verifyGreaterThanOrEqual(min(ape,[],"all"),0);
        end

        function continuousMomentsConvergeButArbitraryDisplacementDoesNot(testCase)
            sizes = [16 32 64 128];
            errors = zeros(numel(sizes),2);
            negativeErrors = zeros(size(sizes));
            for index = 1:numel(sizes)
                n = sizes(index);
                grid = -1+((0:n-1)+.5)*2/n;
                [x,z] = ndgrid(grid);
                [~,materialHeight] = DensityDiagnosticReference.inversePolarTwist(x,z,[.13 -.17],.7,1.3);
                distribution = DensityDiagnosticReference.empiricalDistribution(materialHeight,ones(size(z))*(2/n)^2);
                moments = [sum(distribution.masses.*distribution.values.^2) sum(distribution.masses.*distribution.values.^4)];
                errors(index,:) = abs(moments-[1/3 1/5]);
                % This vertical-only deformation fixes both walls but has
                % a nonunit Jacobian: it is not a parcel rearrangement.
                arbitraryHeight = z-.3*sin(pi*x).*(1-z.^2);
                negativeErrors(index) = abs(mean(arbitraryHeight.^2,"all")-1/3);
                if index == 1
                    % A continuous area-preserving map is NOT an exact
                    % permutation of this finite set of sampled parcels.
                    testCase.verifyGreaterThan(abs(moments(2)-mean(z.^4,"all")),1e-8);
                end
            end
            testCase.verifyLessThan(errors(end,:),errors(1,:)/32);
            testCase.verifyLessThan(errors(end,:),[3e-5 5e-5]);
            testCase.verifyGreaterThan(min(negativeErrors),.02);
            % Exact bias: .3^2 * mean(sin(pi*x)^2) * mean((1-z^2)^2).
            testCase.verifyEqual(negativeErrors(end),.3^2*.5*(8/15),AbsTol=3e-5);
        end

        function changedNonlinearRestProfilesHaveIndependentAnalyticInverses(testCase)
            original = DensityDiagnosticReference.nonlinearRestProfile(N2=.04,curvature=.3);
            changed = DensityDiagnosticReference.nonlinearRestProfile(N2=.07,curvature=-.4);
            z = linspace(-1,1,41).';
            for profile = [original changed]
                density = profile.density(z);
                testCase.verifyLessThan(max(diff(density)),0);
                testCase.verifyGreaterThan(min(profile.localN2(z)),0);
                testCase.verifyEqual(profile.heightForDensity(density),z,AbsTol=2e-13);
                testCase.verifyEqual(profile.buoyancy(z),profile.N2*(z+profile.curvature*z.^2/2),AbsTol=1e-16);
            end
            testCase.verifyGreaterThan(max(abs(changed.density(z)-original.density(z))),1);
            testCase.verifyGreaterThan(abs(original.heightForDensity(changed.density(.5))-.5),.2);
            testCase.verifyError(@()DensityDiagnosticReference.nonlinearRestProfile(curvature=1),'DensityDiagnosticReference:UnstableProfile');
        end

        function productionLinearProfileResolvesBothDisplacementSignsAndTinyEnergy(testCase)
            N2 = .04;
            knots = [-1;-.4;0;.4;1];
            profile = WVNoMotionProfile(knots,1-N2*knots);
            materialHeight = [-.7 -.1 .3;.5 -.25 .8];
            eta = [.2 1e-7 -.2;-.3 -1e-7 .01];
            z = materialHeight+eta;
            testCase.verifyEqual(profile.inverse(1-N2*materialHeight),materialHeight,AbsTol=1e-14);
            actual = profile.availablePotentialEnergy(z,materialHeight,1,1);
            expected = N2*(z-materialHeight).^2/2;
            testCase.verifyEqual(actual,expected,RelTol=1e-13);
            testCase.verifyEqual(profile.availablePotentialEnergy(z,z,1,1),zeros(size(z)));
            testCase.verifyEqual(size(profile.inverse(profile.density(z))),size(z));
            % Use actual representable z-s, not the requested increment:
            % the smallest displacement has no resolvable density anomaly.
            materialHeight = [-.25 -.25 -.25 -.25 -.4 -.4 .4 .4 0 0 -eps(1) eps(1)];
            z = [materialHeight(1:4)+[1e-12 -1e-12 eps(.25) -eps(.25)], ...
                -.4-eps(.4),-.4+eps(.4),.4-eps(.4),.4+eps(.4),eps(1),-eps(1),eps(1),-eps(1)];
            actual = profile.availablePotentialEnergy(z,materialHeight,1,1);
            testCase.verifyGreaterThan(actual,zeros(size(actual)));
            testCase.verifyEqual(actual,N2*(z-materialHeight).^2/2,RelTol=1e-13);
        end

        function productionTwoNodeProfileUsesLinearPolynomial(testCase)
            N2 = .04;
            profile = WVNoMotionProfile([-1;1],1-N2*[-1;1]);
            materialHeight = [-.75 -.2 .2 .75 -.25 -.25 -.25 -.25];
            z = materialHeight+[.3 1e-7 -1e-7 -.3 1e-12 -1e-12 eps(.25) -eps(.25)];
            testCase.verifyEqual(profile.density(z),1-N2*z,AbsTol=3e-16);
            testCase.verifyEqual(profile.inverse(1-N2*materialHeight),materialHeight,AbsTol=1e-14);
            actual = profile.availablePotentialEnergy(z,materialHeight,1,1);
            testCase.verifyGreaterThan(actual,zeros(size(actual)));
            testCase.verifyEqual(actual,N2*(z-materialHeight).^2/2,RelTol=1e-13);
        end

        function productionSharpProfileInvertsAZeroEndpointSlope(testCase)
            profile = WVNoMotionProfile([-1;0;1],[4;3;0]);
            % PCHIP's first endpoint slope is zero for these secants. Its
            % first segment has the explicit density 4-1.5*t^2+.5*t^3.
            t = [0 1e-4 .1 .75 1];
            density = 4-1.5*t.^2+.5*t.^3;
            testCase.verifyEqual(profile.density(-1+t),density,AbsTol=1e-15);
            testCase.verifyEqual(profile.inverse(density),-1+t,AbsTol=2e-11);
            testCase.verifyEqual(profile.inverse([4 0]),[-1 1]);
        end

        function productionCompositionPreservesThreeDimensionalShapes(testCase)
            N2 = .04;
            profile = WVNoMotionProfile([-1;-.4;0;.4;1],1-N2*[-1;-.4;0;.4;1]);
            materialHeight = reshape(linspace(-.7,.7,8),[2 2 2]);
            z = materialHeight+reshape([.2 -.1 .15 -.2 .25 -.3 .1 -.15],[2 2 2]);
            density = profile.density(materialHeight);
            recovered = profile.inverse(density);
            actual = profile.availablePotentialEnergy(z,recovered,1,1);
            testCase.verifyEqual(size(density),[2 2 2]);
            testCase.verifyEqual(size(recovered),[2 2 2]);
            testCase.verifyEqual(size(actual),[2 2 2]);
            testCase.verifyEqual(recovered,materialHeight,AbsTol=1e-14);
            testCase.verifyEqual(actual,N2*(z-materialHeight).^2/2,RelTol=1e-12);
        end

        function productionNonlinearReconstructionConvergesToAnalyticTruth(testCase)
            profiles = [DensityDiagnosticReference.nonlinearRestProfile(N2=.04,curvature=.3), ...
                DensityDiagnosticReference.nonlinearRestProfile(N2=.07,curvature=-.4)];
            [x,z] = ndgrid(linspace(-.8,.8,19));
            [~,materialHeight] = DensityDiagnosticReference.inversePolarTwist(x,z,[.13 -.17],.7,1.3);
            for reference = profiles
                errors = zeros(3,2);
                counts = [17 65 257];
                for index = 1:numel(counts)
                    knots = linspace(-1,1,counts(index)).';
                    profile = WVNoMotionProfile(knots,reference.density(knots));
                    recovered = profile.inverse(reference.density(materialHeight));
                    errors(index,1) = max(abs(recovered-materialHeight),[],"all");
                    actual = profile.availablePotentialEnergy(z,materialHeight,reference.gravity,reference.rho0);
                    expected = reference.availablePotentialEnergy(z,materialHeight);
                    errors(index,2) = max(abs(actual-expected),[],"all");
                    testCase.verifyGreaterThanOrEqual(min(actual,[],"all"),0);
                end
                testCase.verifyLessThan(errors(end,:),errors(1,:)/100);
                testCase.verifyLessThan(errors(end,:),[1e-8 1e-9]);
            end
        end

        function productionWeakGradientAccuracyAccountsForDensityResolution(testCase)
            rho0 = 1025;
            g = 9.81;
            N2 = 1e-10;
            knots = [-1000;-731;-207;0];
            density = @(z) rho0*(1-N2*z/g);
            profile = WVNoMotionProfile(knots,density(knots));
            materialHeight = [-850 -500 -100];
            z = materialHeight+[25 -10 2];
            densityHeightResolution = eps(rho0)*g/(rho0*N2);
            testCase.verifyEqual(profile.inverse(density(materialHeight)),materialHeight,AbsTol=4*densityHeightResolution);
            actual = profile.availablePotentialEnergy(z,materialHeight,g,rho0);
            testCase.verifyEqual(actual,N2*(z-materialHeight).^2/2,RelTol=2e-6);
        end

        function productionLargeArrayCrossesChunkBoundariesAndManyIntervals(testCase)
            N2 = .04;
            knots = linspace(-1,1,33).';
            profile = WVNoMotionProfile(knots,1-N2*knots);
            shape = [17 17 241];
            index = (1:prod(shape)).';
            materialHeight = -1+2*mod(index,257)/256;
            z = -materialHeight;
            z(mod(index,3)==0) = materialHeight(mod(index,3)==0);
            z(mod(index,3)==1) = knots(1+mod(index(mod(index,3)==1),numel(knots)));
            % Distinct cases around the 65536-element chunk boundary and
            % final tail include zero, upward, downward, and domain-wide moves.
            selected = [1 65535 65536 65537 numel(index)];
            materialHeight(selected) = [-1 .5 -.5 .25 1];
            z(selected) = [1 .5 .75 -.75 -1];
            materialHeight = reshape(materialHeight,shape);
            z = reshape(z,shape);
            actual = profile.availablePotentialEnergy(z,materialHeight,1,1);
            expected = N2*(z-materialHeight).^2/2;
            testCase.verifyEqual(size(actual),shape);
            testCase.verifyEqual(actual,expected,RelTol=2e-13);
            testCase.verifyEqual(actual(z==materialHeight),zeros(nnz(z==materialHeight),1));
            testCase.verifyGreaterThan(actual(z~=materialHeight),zeros(nnz(z~=materialHeight),1));
        end

        function productionIrregularIntervalsMatchIndependentQuadrature(testCase)
            knots = [-2;-1.7;-.8;-.79;-.1;.25;1];
            density = [9;8.99;8.7;6;5.9;3;1];
            profile = WVNoMotionProfile(knots,density);
            % Differentiate MATLAB's density PCHIP directly and integrate
            % numerically; do not use the production primitive or its
            % normalized coefficients as the reference integration method.
            pp = pchip(knots,density);
            derivative = mkpp(pp.breaks,pp.coefs(:,1:end-1).*(pp.order-1:-1:1));
            materialHeight = reshape([-2 -1.7 -.8 -.795 -.79 -.1 .25 1 -.25 .8 -.1 -.1],3,4);
            z = reshape([1 -.1 .25 -.794 -2 .25 -.8 -2 .9 -.7 -.1 -.8],3,4);
            g = 9.81;
            rho0 = 1025;
            expected = zeros(size(z));
            for index = 1:numel(z)
                lower = min(z(index),materialHeight(index));
                upper = max(z(index),materialHeight(index));
                if lower == upper
                    continue
                end
                waypoints = knots(knots>lower & knots<upper);
                integrand = @(r) (r-z(index)).*ppval(derivative,r);
                expected(index) = (g/rho0)*sign(z(index)-materialHeight(index)) * ...
                    integral(integrand,lower,upper,Waypoints=waypoints,AbsTol=1e-12,RelTol=1e-12);
            end
            actual = profile.availablePotentialEnergy(z,materialHeight,g,rho0);
            testCase.verifyEqual(actual,expected,AbsTol=5e-13,RelTol=5e-12);
            testCase.verifyGreaterThanOrEqual(actual,zeros(size(actual)));
        end

        function productionLargeInverseUsesAnalyticNonlinearDensityAcrossChunks(testCase)
            profile = WVNoMotionProfile([-1;0;1],[4;3;0]);
            % These two cubic formulas independently describe this PCHIP.
            % The first has zero derivative at -1, so near-endpoint queries
            % take a different inversion path from exact knots and interiors.
            densityForHeight = @(h) (h<0).*(4-1.5*(h+1).^2+.5*(h+1).^3) + ...
                (h>=0).*(3-1.5*h-2*h.^2+.5*h.^3);
            shape = [129 129 65];
            index = (1:prod(shape)).';
            expectedHeight = -1+2*mod(73*index,1025)/1024;
            selected = [1 1048575 1048576 1048577 numel(index)];
            expectedHeight(selected) = [-1 -1+1e-4 0 1 -.2];
            expectedHeight = reshape(expectedHeight,shape);
            density = densityForHeight(expectedHeight);
            actual = profile.inverse(density);
            testCase.verifyEqual(size(actual),shape);
            testCase.verifyEqual(actual,expectedHeight,AbsTol=2e-11);
            testCase.verifyEqual(densityForHeight(actual),density,AbsTol=1e-13);
            testCase.verifyEqual(actual(selected([1 3 4])),[-1 0 1]);
            % Grouping invariance supplements the analytic-height oracle;
            % matching a scalar call alone would not establish correctness.
            for element = selected
                testCase.verifyEqual(actual(element),profile.inverse(density(element)),AbsTol=2e-14);
            end
        end

        function productionInversePreservesRowAndEmptyShapes(testCase)
            knots = [-2;-1.7;-.8;-.79;-.1;.25;1];
            densities = [9;8.99;8.7;6;5.9;3;1];
            profile = WVNoMotionProfile(knots,densities);
            expectedHeight = [-2 -1.99 -1.7 -.795 -.1 .07 .25 .9 1];
            reference = pchip(knots,densities);
            density = ppval(reference,expectedHeight);
            actual = profile.inverse(density);
            testCase.verifyEqual(size(actual),size(expectedHeight));
            testCase.verifyEqual(actual,expectedHeight,AbsTol=2e-10);
            testCase.verifyEqual(ppval(reference,actual),density,AbsTol=1e-13);
            for shape = {[0 3],[2 0 4],[1 0]}
                emptyDensity = zeros(shape{1});
                testCase.verifyEqual(profile.inverse(emptyDensity),emptyDensity);
            end
        end

        function productionInverseEndpointRoundoffRemainsInsideProfile(testCase)
            % Both exact endpoint queries previously returned one ulp above
            % the original upper knot, making a valid APE composition fail.
            profiles = {WVNoMotionProfile([-.1;.3],[2;1]), ...
                WVNoMotionProfile([-1000;-1000/9;0],[1;.010699594101397314;0])};
            for index = 1:numel(profiles)
                profile = profiles{index};
                densitySpacing = eps(max(abs(profile.rho)));
                density = [profile.rho(1),profile.rho(end),profile.rho(end)-4*densitySpacing];
                materialHeight = profile.inverse(density);
                testCase.verifyEqual(materialHeight,[profile.z(1),profile.z(end),profile.z(end)]);
                testCase.verifyGreaterThanOrEqual(materialHeight,profile.z(1));
                testCase.verifyLessThanOrEqual(materialHeight,profile.z(end));
                height = mean(profile.z([1 end]))*ones(size(materialHeight));
                ape = profile.availablePotentialEnergy(height,materialHeight,1,1);
                testCase.verifyTrue(all(isfinite(ape)));
                testCase.verifyGreaterThanOrEqual(ape,0);
                if index == 1
                    expected = (profile.rho(1)-profile.rho(end))/(2*diff(profile.z))*(height-materialHeight).^2;
                    testCase.verifyEqual(ape,expected,RelTol=2e-14);
                end
                testCase.verifyError(@()profile.inverse(profile.rho(end)-9*densitySpacing), ...
                    'WVNoMotionProfile:DensityOutsideProfile');
                outside = profile.z(end)+eps(max(abs(profile.z)));
                testCase.verifyError(@()profile.availablePotentialEnergy(height(1),outside,1,1), ...
                    'WVNoMotionProfile:HeightOutsideProfile');
            end
        end

        function productionRejectsPlateausAndOutOfRangeQueries(testCase)
            testCase.verifyError(@()WVNoMotionProfile([-1;0;1],[3;2;2]),'WVNoMotionProfile:NonInvertibleDensity');
            testCase.verifyError(@()WVNoMotionProfile([-1;0;1],[1;2;3]),'WVNoMotionProfile:NonInvertibleDensity');
            profile = WVNoMotionProfile([-1;0;1],[3;2;1]);
            testCase.verifyError(@()profile.inverse(.9),'WVNoMotionProfile:DensityOutsideProfile');
            testCase.verifyError(@()profile.inverse(3.1),'WVNoMotionProfile:DensityOutsideProfile');
            testCase.verifyError(@()profile.density(1.01),'WVNoMotionProfile:HeightOutsideProfile');
            testCase.verifyError(@()profile.availablePotentialEnergy(0,-1.01,1,1),'WVNoMotionProfile:HeightOutsideProfile');
            testCase.verifyEqual(profile.inverse([3 1]),[-1 1]);
        end
    end
end
