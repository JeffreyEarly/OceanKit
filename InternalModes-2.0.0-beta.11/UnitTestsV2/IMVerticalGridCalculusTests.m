classdef IMVerticalGridCalculusTests < matlab.unittest.TestCase
    % Check sampled calculus against independently known physical functions.

    methods (TestClassSetup)
        function addProviderPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(repoRoot));
        end
    end

    methods (Test)
        function physicalDepthCalculusPreservesComplexColumns(testCase)
            grid = IMVerticalGrid.physicalDepth([-8 -2],17);
            calculus = grid.calculus;
            zQuery = [-8;-7.25;-4;-2];
            s = @(z) (z+8)/6;
            profiles = @(t) [1+2*t-t.^3,3-1i*t.^2];
            derivatives = @(t) [(2-3*t.^2)/6,-1i*t/3];
            values = profiles(s(grid.z));

            testCase.verifyEqual(calculus.interpolationMatrix(zQuery)*values,profiles(s(zQuery)),AbsTol=1e-12);
            testCase.verifyEqual(calculus.firstDerivativeMatrix()*values,derivatives(s(grid.z)),AbsTol=1e-11);
            testCase.verifyEqual(calculus.endpointTraces(values),profiles([0;1]));
            testCase.verifyEqual(calculus.endpointTraces(values,derivativeOrder=1),derivatives([0;1]),AbsTol=1e-11);

            t = s(zQuery);
            expectedIntegral = 6*[7/4-t-t.^2+t.^4/4,3*(1-t)-1i*(1-t.^3)/3];
            testCase.verifyEqual(calculus.integrateFromSurface(values,zQuery),expectedIntegral,AbsTol=1e-11);
        end

        function exponentialWKBCalculusUsesThePhysicalJacobian(testCase)
            zDomain = [-1350 -150];
            b = 300;
            exponential = @(z) exp((z-zDomain(2))/b);
            background = IMBackground(N2=@(z) 1e-4*exponential(z).^2,zDomain=zDomain,surfaceCondition="freeSurface");
            grid = IMVerticalGrid.wkb(background,33);
            calculus = grid.calculus;
            zQuery = [zDomain(1);-1160;-725;-281;zDomain(2)];
            profiles = @(e) [e.^2+1i*e,1-2*e];
            derivatives = @(e) [2*e.^2+1i*e,-2*e]/b;
            values = profiles(exponential(grid.z));

            testCase.verifyEqual(calculus.firstDerivativeMatrix()*values,derivatives(exponential(grid.z)),AbsTol=5e-12);
            testCase.verifyEqual(calculus.interpolationMatrix(zQuery)*values,profiles(exponential(zQuery)),AbsTol=5e-12);
            testCase.verifyEqual(calculus.endpointTraces(values,derivativeOrder=1),derivatives(exponential(zDomain(:))),AbsTol=5e-12);
            e = exponential(zQuery);
            expectedIntegral = [b*(1-e.^2)/2+1i*b*(1-e),zDomain(2)-zQuery-2*b*(1-e)];
            testCase.verifyEqual(calculus.integrateFromSurface(values,zQuery),expectedIntegral,AbsTol=2e-7);

            [zQuad,wQuad] = calculus.quadratureRule(9);
            eBottom = exponential(zDomain(1));
            xQuad = 2*(exponential(zQuad)-eBottom)/(1-eBottom)-1;
            polynomial = legendre(9,xQuad.');
            testCase.verifyLessThan(max(abs(polynomial(1,:))),2e-11);
            testCase.verifyGreaterThan(min(wQuad),0);
            testCase.verifyEqual(wQuad.'*exponential(zQuad),b*(1-eBottom),AbsTol=2e-9);
        end

        function strongSmoothedBandwidthMatchesIndependentFilterAndWeakPairing(testCase)
            zDomain = [-1000 0];
            D = diff(zDomain);
            Ls = 90;
            k = pi/D;
            meanN2 = 5e-5;
            amplitude = 0.995*meanN2;
            N2 = @(z) meanN2+amplitude*cos(k*(z-zDomain(1)));
            background = IMBackground(N2=N2,zDomain=zDomain,rotationRate=0,surfaceCondition="freeSurface");
            grid = IMVerticalGrid.smoothedBandwidth(background,65,smoothingLength=Ls);
            calculus = grid.calculus;

            % Zero rotation gives B=N2. This cosine has an exact Neumann-filter response.
            filteredAmplitude = amplitude/(1+(Ls*k)^2);
            metric = @(z) sqrt(meanN2+filteredAmplitude*cos(k*(z-zDomain(1))));
            coordinateLength = integral(metric,zDomain(1),zDomain(2),RelTol=1e-13,AbsTol=1e-13);
            coordinate = @(z) arrayfun(@(depth) integral(metric,zDomain(1),depth,RelTol=1e-13,AbsTol=1e-13),z)/coordinateLength;
            zQuery = [-998;-849;-513;-191;-4];
            s = coordinate(grid.z);
            sQuery = coordinate(zQuery);
            values = [1+s,s.^2];
            expectedDerivatives = (metric(grid.z)/coordinateLength).*[ones(size(s)),2*s];
            Dz = calculus.firstDerivativeMatrix();

            testCase.verifyEqual(Dz*values,expectedDerivatives,AbsTol=2e-10);
            testCase.verifyEqual(calculus.interpolationMatrix(zQuery)*values,[1+sQuery,sQuery.^2],AbsTol=2e-11);
            testCase.verifyEqual(calculus.integrateFromSurface(ones(size(s)),zQuery),-zQuery,AbsTol=2e-6);

            % Verify the physical integration-by-parts identity on resolved profiles.
            u = values(:,1);
            v = values(:,2);
            weakBoundary = sum(grid.weights.*((Dz*u).*v+u.*(Dz*v)));
            testCase.verifyEqual(weakBoundary,u(end)*v(end)-u(1)*v(1),AbsTol=2e-8);

            [zQuad,wQuad] = calculus.quadratureRule(11);
            polynomial = legendre(11,2*coordinate(zQuad).'-1);
            testCase.verifyLessThan(max(abs(polynomial(1,:))),2e-10);
            testCase.verifyEqual(wQuad.'*metric(zQuad),coordinateLength,AbsTol=2e-10);

            nFine = 41;
            [zFine,wFine,DzFine] = calculus.lobattoRule(nFine);
            sFine = coordinate(zFine);
            expectedNodes = (1-cos(pi*(0:nFine-1).'/(nFine-1)))/2;
            testCase.verifyEqual(zFine([1 end]).',zDomain);
            testCase.verifyEqual(sFine,expectedNodes,AbsTol=2e-11);
            testCase.verifyEqual(sum(wFine),D,AbsTol=5e-10);
            testCase.verifyEqual(DzFine*sFine,metric(zFine)/coordinateLength,AbsTol=2e-10);

            restored = IMLobattoCalculus.fromSnapshot(calculus.snapshot());
            testCase.verifyEqual(restored.firstDerivativeMatrix(),Dz);
            testCase.verifyEqual(restored.interpolationMatrix(zQuery),calculus.interpolationMatrix(zQuery));
            [restoredZ,restoredWeights] = restored.quadratureRule(11);
            testCase.verifyEqual(restoredZ,zQuad);
            testCase.verifyEqual(restoredWeights,wQuad);
        end

        function calculusRequiresAnExplicitRepresentation(testCase)
            grid = IMVerticalGrid.physicalDepth([-1 0],9);
            geometryOnly = IMVerticalGrid(z=grid.z,weights=grid.weights,kind=grid.kind,design=grid.design);
            testCase.verifyEmpty(geometryOnly.calculus);

            calculus = grid.calculus;
            snapshot = calculus.snapshot();
            snapshot.version = [];
            testCase.verifyError(@() IMLobattoCalculus.fromSnapshot(snapshot),"IMLobattoCalculus:InvalidSnapshot");
            testCase.verifyError(@() calculus.interpolationMatrix([-2;-0.5]),"IMLobattoCalculus:OutsideDomain");
        end
    end
end
