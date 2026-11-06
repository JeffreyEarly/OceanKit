classdef IMHydrostaticRootTransformTests < matlab.unittest.TestCase
    % Independent checks of IMInternalModesDiscreteTransform.hydrostaticRoot (issue #55).

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot);
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function constantStratificationIsScaledDCTAndDST(testCase)
            N0 = 5e-3; D = 1000;
            background = IMBackground.constantStratification(N0=N0,zDomain=[-D 0],surfaceCondition="rigidLid");
            for nPoints = [17 33]
                [transform,grid,assessment] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,nPoints);
                r = nPoints-2;
                dz = D/(nPoints-1);
                % Uniform root points and trapezoid weights.
                testCase.verifyEqual(grid.z,linspace(-D,0,nPoints).',AbsTol=1e-9*D)
                testCase.verifyEqual(grid.weights,dz*[0.5;ones(nPoints-2,1);0.5],RelTol=1e-9)
                testCase.verifyEqual(transform.weights,grid.weights)
                % Full active band, aligned identities and masks.
                testCase.verifyEqual(assessment.retainedInternalModeCount,r)
                testCase.verifyEqual(assessment.limitingCriterion,"maximumRootGridBand")
                testCase.verifyEqual(transform.modeNumber,0:r)
                testCase.verifyEqual(transform.activeModeMask(variable="F"),true(1,r+1))
                testCase.verifyEqual(transform.activeModeMask(variable="G"),[false true(1,r)])
                testCase.verifyLessThan(assessment.final.energyErrorF,1e-10)
                testCase.verifyLessThan(assessment.final.energyErrorG,1e-10)
                % Analytical shapes, equivalent depths and F = h G'.
                g = background.g; z = grid.z; j = 1:r;
                h = N0^2*D^2./(g*(j*pi).^2);
                amplitude = sqrt(2*g./(N0^2*h));
                Gexact = amplitude.*sin((z+D)*(j*pi/D));
                Fexact = h.*amplitude.*(j*pi/D).*cos((z+D)*(j*pi/D));
                F = transform.inverseMatrix(variable="F");
                G = transform.inverseMatrix(variable="G");
                testCase.verifyEqual(transform.h(2:end),h,RelTol=1e-10)
                testCase.verifyEqual(F(:,1),ones(nPoints,1),AbsTol=1e-12)
                testCase.verifyEqual(G(:,1),zeros(nPoints,1))
                % One shared orientation per mode: F and G scale together.
                orientation = sign(sum(G(:,2:end).*Gexact,1));
                testCase.verifyEqual(abs(orientation),ones(1,r))
                testCase.verifyEqual(F(:,2:end),Fexact.*orientation,AbsTol=1e-9*max(abs(Fexact(:))))
                testCase.verifyEqual(G(:,2:end),Gexact.*orientation,AbsTol=1e-9*max(amplitude))
                testCase.verifyEqual(G([1 end],:),zeros(2,r+1))
            end
        end

        function exponentialStratificationMatchesAnalyticalModes(testCase)
            background = IMBackground.exponentialStratification(N0=5.2e-3,b=1300,zDomain=[-4000 0],surfaceCondition="rigidLid");
            [transform,grid,assessment] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,33);
            r = assessment.retainedInternalModeCount;
            solution = background.analyticalSolution();
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) solution.N2(z),zDomain=background.zDomain, ...
                rotationRate=background.rotationRate,latitude=background.latitude,g=background.g,surfaceBoundary=IMBoundaryCondition.dirichlet());
            exact = solution.internalModes(evp,nModes=r,normalization="rmsF");
            testCase.verifyEqual(transform.h(2:end),exact.h(1:r),RelTol=1e-8)
            F = transform.inverseMatrix(variable="F"); G = transform.inverseMatrix(variable="G");
            Fexact = exact.F(grid.z,normalization="rmsF"); Gexact = exact.G(grid.z,normalization="rmsF");
            orientation = sign(sum(G(:,2:end).*Gexact,1));
            testCase.verifyEqual(F(:,2:end),Fexact.*orientation,AbsTol=1e-6*max(abs(Fexact(:))))
            testCase.verifyEqual(G(:,2:end),Gexact.*orientation,AbsTol=1e-6*max(abs(Gexact(:))))
            testCase.verifyLessThanOrEqual(assessment.final.energyErrorF,assessment.gramTolerance)
            testCase.verifyLessThanOrEqual(assessment.final.energyErrorG,assessment.gramTolerance)
        end

        function rootPointsReproduceV4Grids(testCase)
            fixture = readtable(fullfile(fileparts(mfilename("fullpath")),"Fixtures","IMHydrostaticRootTransformV4Points.csv"),TextType="string");
            profiles = struct(gmExponential=@(z) (5.2e-3)^2*exp(2*z/1300),pycnocline=@(z) 1e-6+4e-5*exp(-((z+300)/80).^2));
            cases = unique(fixture(:,["profile","depth","nPoints"]),"rows");
            for k = 1:height(cases)
                rows = sortrows(fixture(fixture.profile==cases.profile(k) & fixture.nPoints==cases.nPoints(k),:),"node");
                background = IMBackground(N2=profiles.(cases.profile(k)),zDomain=[-cases.depth(k) 0],latitude=33,surfaceCondition="rigidLid");
                [~,grid] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,cases.nPoints(k));
                % v4 roots carry their own unqualified resolution error, so
                % agreement is required to the declared root tolerance.
                testCase.verifyEqual(grid.z,rows.z,cases.profile(k)+" "+cases.nPoints(k),AbsTol=1e-6*cases.depth(k))
            end
        end

        function stronglyVaryingProfileIsIndependentOfResolution(testCase)
            background = IMBackground(N2=@(z) 1e-6+4e-5*exp(-((z+300)/80).^2),zDomain=[-2000 0],latitude=33,surfaceCondition="rigidLid");
            [T1,grid1,A1] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,33);
            % The default path raises resolution until an independent refinement agrees.
            testCase.verifyGreaterThan(height(A1.convergence.attempts),1)
            testCase.verifyEqual(A1.convergence.attempts.nEVP(end),A1.convergence.candidateNEVP)
            [T2,grid2,A2] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,33, ...
                solver=IMSolverSpectral(nEVP=2*A1.convergence.candidateNEVP,coordinateKind=A1.convergence.coordinateKind));
            testCase.verifyEqual(height(A2.convergence.attempts),1)
            testCase.verifyEqual(A2.retainedInternalModeCount,A1.retainedInternalModeCount)
            testCase.verifyEqual(grid2.z,grid1.z,AbsTol=1e-8*2000)
            testCase.verifyEqual(T2.h,T1.h,RelTol=1e-8)
            testCase.verifyEqual(grid2.weights,grid1.weights,AbsTol=1e-6*max(grid1.weights))
            testCase.verifyLessThanOrEqual(A1.convergence.maximumShapeError,A1.convergence.tolerance)
            testCase.verifyLessThanOrEqual(A1.convergence.rootDisplacement,A1.convergence.tolerance)
            % A supplied solver is used as given and must pass the check itself.
            testCase.verifyError(@() IMInternalModesDiscreteTransform.hydrostaticRoot(background,33,solver=IMSolverSpectral(nEVP=128,coordinateKind="wkb")), ...
                "IMInternalModesDiscreteTransform:UnconvergedHydrostaticRootModes")
        end

        function automaticSelectionRecordsSearchAndLimit(testCase)
            background = IMBackground(N2=@(z) 1e-4*exp(2*z/400),zDomain=[-1000 0],surfaceCondition="rigidLid");
            [transform,grid,assessment] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,33);
            search = assessment.search;
            r = assessment.retainedInternalModeCount;
            testCase.verifyEqual(assessment.selection,"automatic")
            testCase.verifyEqual(search.internalModeCount.',31:-1:r)
            testCase.verifyEqual(search.accepted.',[false(1,height(search)-1) true])
            testCase.verifyLessThan(r,31)
            testCase.verifyEqual(assessment.limitingCriterion,search.rejectionReason(end-1))
            testCase.verifyEqual(assessment.final,search(end,:))
            testCase.verifyEqual(numel(transform.modeNumber),r+1)
            testCase.verifyEqual(grid.kind,"modeRoots")
            testCase.verifyEqual(grid.design.retainedInternalModeCount,r)
            % The returned weights are the final band's fit, not a full-band fit.
            Phi = transform.inverseMatrix(variable="F");
            gram = Phi.'*(grid.weights.*Phi);
            H = 1000*ones(size(transform.h));
            epsF = norm((gram-diag(H))./sqrt(H.'*H),2);
            testCase.verifyEqual(epsF,assessment.final.energyErrorF,RelTol=1e-8)
        end

        function strictRequestsAreExact(testCase)
            background = IMBackground(N2=@(z) 1e-4*exp(2*z/400),zDomain=[-1000 0],surfaceCondition="rigidLid");
            [transform,~,assessment] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,33,nInternalModes=20);
            testCase.verifyEqual(assessment.selection,"strict")
            testCase.verifyEqual(assessment.retainedInternalModeCount,20)
            testCase.verifyEqual(height(assessment.search),1)
            testCase.verifyEqual(transform.modeNumber,0:20)
            testCase.verifyError(@() IMInternalModesDiscreteTransform.hydrostaticRoot(background,33,nInternalModes=31),"IMInternalModesDiscreteTransform:HydrostaticRootBandRejected")
            testCase.verifyError(@() IMInternalModesDiscreteTransform.hydrostaticRoot(background,33,nInternalModes=32),"IMInternalModesDiscreteTransform:InvalidHydrostaticRootBand")
        end

        function galerkinRecoveryDiffersFromPrescribedDual(testCase)
            background = IMBackground(N2=@(z) 1e-4*exp(2*z/400),zDomain=[-1000 0],surfaceCondition="rigidLid");
            [transform,grid] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,33,forwardOperator="galerkin");
            testCase.verifyEqual(transform.forwardOperator(variable="F"),"galerkin")
            Phi = transform.inverseMatrix(variable="F");
            coefficients = linspace(1,-1,size(Phi,2)).';
            values = Phi*coefficients;
            galerkin = transform.transformForward(values,variable="F");
            H = 1000*ones(numel(transform.h),1);
            prescribed = transform.modeProjectionFunctional(values,variable="F")./H;
            testCase.verifyEqual(galerkin,coefficients,AbsTol=1e-10*norm(coefficients))
            testCase.verifyGreaterThan(norm(prescribed-coefficients),1e3*norm(galerkin-coefficients))
            % A nonmodal input: Galerkin and prescribed-dual coefficients differ.
            x = exp(grid.z/200);
            testCase.verifyGreaterThan(norm(transform.transformForward(x,variable="F")-transform.modeProjectionFunctional(x,variable="F")./H),0)
        end

        function truncatedCollocationRejectsModesWithinTheRootSpan(testCase)
            background = IMBackground(N2=@(z) 1e-6+4e-5*exp(-((z+300)/80).^2),zDomain=[-2000 0],latitude=33,surfaceCondition="rigidLid");
            nPoints = 33;
            [collocation,grid,assessment] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,nPoints);
            galerkin = IMInternalModesDiscreteTransform.hydrostaticRoot(background,nPoints,forwardOperator="galerkin");
            r = assessment.retainedInternalModeCount;
            testCase.verifyEqual(assessment.forwardOperator,"truncatedCollocation")
            testCase.verifyEqual(collocation.forwardOperator(variable="F"),"truncatedCollocation")
            testCase.verifyEqual(collocation.forwardOperator(variable="G"),"truncatedCollocation")
            % Same band, weights and energy diagnostics; only the operator differs.
            testCase.verifyEqual(galerkin.weights,collocation.weights)
            testCase.verifyEqual(collocation.relativeGramOperatorError(variable="F"),galerkin.relativeGramOperatorError(variable="F"),RelTol=1e-12)
            % Independent modes over the whole span, at the converged resolution.
            solver = IMSolverSpectral(nEVP=assessment.convergence.candidateNEVP,coordinateKind=assessment.convergence.coordinateKind);
            evp = IMInternalModes.hydrostaticGModes(N2=background.N2,zDomain=background.zDomain,latitude=33,surfaceBoundary=IMBoundaryCondition.dirichlet());
            basis = solver.solveEVP(evp,nModes=nPoints-1);
            F = basis.F(grid.z,normalization="rmsF");
            G = basis.G(grid.z,normalization="rmsF"); G([1 end],:) = 0;
            withinSpanF = F(:,r+1:nPoints-1);
            withinSpanG = G(:,r+1:nPoints-2);
            HF = 2000*ones(r+1,1);
            leakF = @(T) norm(sqrt(HF).*(T.transformForward(withinSpanF,variable="F"))./sqrt(2000),2);
            leakG = @(T) norm(sqrt([0 2000./collocation.h(2:end)]).'.*T.transformForward(withinSpanG,variable="G")./sqrt(2000./basis.h(r+1:nPoints-2)),2);
            testCase.verifyLessThan(leakF(collocation),1e-9)
            testCase.verifyLessThan(leakG(collocation),1e-9)
            testCase.verifyGreaterThan(leakF(galerkin),1e-3)
            % Retained modes round-trip exactly under both operators.
            retained = [ones(nPoints,1) F(:,1:r)];
            testCase.verifyEqual(collocation.transformForward(retained,variable="F"),eye(r+1),AbsTol=1e-9)
            testCase.verifyEqual(galerkin.transformForward(retained,variable="F"),eye(r+1),AbsTol=1e-9)
            testCase.verifyLessThan(collocation.roundTripError(variable="G"),1e-9)
            % Recorded diagnostics use the selected operator.
            testCase.verifyLessThan(assessment.final.omittedGainF,1e-9)
            testCase.verifyGreaterThan(assessment.collocation.conditionG,1)
            % Prefixes keep the leading collocation rows.
            prefix = collocation.prefixTransform(3);
            testCase.verifyEqual(prefix.forwardOperator(variable="F"),"truncatedCollocation")
            fullForward = collocation.forwardMatrix(variable="F");
            testCase.verifyEqual(prefix.forwardMatrix(variable="F"),fullForward(1:3,:))
        end

        function returnedContinuousBasisAndFixedGridBridgeAgree(testCase)
            background = IMBackground.constantStratification(N0=5e-3,zDomain=[-1000 0],surfaceCondition="rigidLid");
            [T,grid,~,basis] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,17,nInternalModes=8);
            testCase.verifyClass(basis,"IMHydrostaticModesBasis")
            testCase.verifyEqual(basis.family,"hydrostaticFModes")
            testCase.verifyEqual(basis.normalization,"rmsF")
            testCase.verifyEqual(basis.columnLabels,["barotropic:0" "positive:"+string(1:8)])
            testCase.verifyEqual(T.basis,basis)
            z = linspace(-1000,0,73).';
            F = basis.F(z); Gz = basis.Gz(z);
            testCase.verifyEqual(F(:,2:end),basis.h(2:end).*Gz(:,2:end),AbsTol=1e-9)
            testCase.verifyEqual(Gz(:,1),zeros(73,1))
            for variable = ["F","G"]
                testCase.verifyEqual(T.inverseMatrix(variable=variable),basis.(variable)(grid.z))
                testCase.verifyEqual(T.targetGramMatrix(variable=variable),basis.gramMatrix(variable=variable))
                testCase.verifyEqual(T.targetMajorantGramMatrix(variable=variable),basis.majorantGramMatrix(variable=variable))
                bridge = basis.discreteTransform(grid);
                collection = IMBasisCollection({basis});
                projected = collection.projectionOnGrid(grid.z,grid.weights,variable=variable,prepareRecovery=true);
                testCase.verifyEqual(projected.recovery.forwardMatrix,bridge.forwardMatrix(variable=variable),AbsTol=1e-10)
                testCase.verifyEqual(bridge.inverseMatrix(variable=variable),T.inverseMatrix(variable=variable))
                testCase.verifyLessThan(bridge.roundTripError(variable=variable),1e-10)
                context = T.quadratureFitContext(variable=variable);
                testCase.verifyEqual(context.interiorWeight,ones(17,1)*(variable=="F")+background.N2(grid.z)/background.g*(variable=="G"))
                changed = T.transformWithWeights(1.1*T.weights);
                testCase.verifyEqual(changed.basis,basis)
                testCase.verifyEqual(changed.forwardMatrix(variable=variable),T.forwardMatrix(variable=variable))
                testCase.verifyEqual(changed.gramMatrix(variable=variable),1.1*T.gramMatrix(variable=variable),AbsTol=1e-8)
            end
            selected = basis.selectColumns([5 1 3]);
            testCase.verifyEqual(selected.F(z),F(:,[5 1 3]))
            sampled = selected.discreteTransform(grid);
            testCase.verifyEqual(sampled.modeNumber,[4 0 2])
            testCase.verifyEqual(sampled.activeModeMask(variable="G"),[true false true])
            prefix = T.prefixTransform(3);
            testCase.verifyEqual(prefix.basis.columnLabels,basis.columnLabels(1:3))
            testCase.verifyEqual(prefix.problemMetadata.columnLabels,prefix.basis.columnLabels)
        end

        function maxAbsNormalizationHasConsistentContinuousMetrics(testCase)
            background = IMBackground.constantStratification(N0=5e-3,zDomain=[-1000 0],surfaceCondition="rigidLid");
            grid = IMVerticalGrid.physicalDepth(background.zDomain,65);
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=64)}
                if isa(solver{1},"IMSolverAnalytical")
                    problem = IMInternalModes.hydrostaticFModes(background);
                    basis = solver{1}.solveEVP(problem,nModes=6).withNormalization("maxAbsF");
                else
                    problem = IMInternalModes.hydrostaticGModes(background);
                    basis = IMHydrostaticModesBasis(solver{1}.solveEVP(problem,nModes=5),normalization="maxAbsF");
                end
                T = basis.discreteTransform(grid);
                collection = IMBasisCollection({basis});
                for variable = ["F","G"]
                    projected = collection.projectionOnGrid(grid.z,grid.weights,variable=variable,prepareRecovery=true);
                    testCase.verifyEqual(projected.recovery.forwardMatrix,T.forwardMatrix(variable=variable),AbsTol=1e-10)
                    testCase.verifyEqual(projected.targetGramMatrix,T.targetGramMatrix(variable=variable),AbsTol=1e-8)
                    testCase.verifyLessThan(T.roundTripError(variable=variable),1e-10)
                end
                testCase.verifyEqual(diag(T.targetGramMatrix(variable="F")).',[1000 500*ones(1,5)],RelTol=1e-10)
            end
        end

        function customScalingPreservesPhysicalMetricsAndOffSpanCoefficients(testCase)
            background = IMBackground.constantStratification(N0=5e-3,zDomain=[-1000 0],surfaceCondition="rigidLid");
            grid = IMVerticalGrid.physicalDepth(background.zDomain,65);
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=64)}
                if isa(solver{1},"IMSolverAnalytical")
                    basis = solver{1}.solveEVP(IMInternalModes.hydrostaticFModes(background),nModes=6);
                else
                    basis = IMHydrostaticModesBasis(solver{1}.solveEVP(IMInternalModes.hydrostaticGModes(background),nModes=5));
                end
                factors = basis.normalizationFactors("rmsF");
                scaled = basis.addNormalization("halfAmplitude",@(~,j) 2*factors(j)).withNormalization("halfAmplitude");
                original = basis.discreteTransform(grid);
                changed = scaled.discreteTransform(grid);
                input = exp(grid.z/170)+1i*sin(grid.z/213);
                for variable = ["F","G"]
                    testCase.verifyEqual(changed.inverseMatrix(variable=variable),0.5*original.inverseMatrix(variable=variable),AbsTol=1e-12);
                    testCase.verifyEqual(changed.targetGramMatrix(variable=variable),0.25*original.targetGramMatrix(variable=variable),RelTol=1e-11,AbsTol=1e-10);
                    testCase.verifyEqual(changed.targetGramMatrix(variable=variable),scaled.gramMatrix(variable=variable),RelTol=1e-11,AbsTol=1e-8);
                    physical = scaled.projectionOnGrid(grid.z,grid.weights,variable=variable);
                    reference = basis.projectionOnGrid(grid.z,grid.weights,variable=variable);
                    testCase.verifyEqual(physical.projectPhysical(input),2*reference.projectPhysical(input),AbsTol=1e-10);
                    testCase.verifyEqual(changed.transformForward(input,variable=variable),2*original.transformForward(input,variable=variable),AbsTol=1e-10);
                end
            end
        end

        function finiteDifferenceBridgePreservesPhysicalGramAndOffSpanDual(testCase)
            background = IMBackground.constantStratification(N0=5e-3,zDomain=[-1000 0],surfaceCondition="rigidLid");
            solver = IMSolverFiniteDifference(z=linspace(-1000,0,17).');
            basis = solver.solveEVP(IMInternalModes.hydrostaticFModes(background),nModes=5);
            grid = IMVerticalGrid.physicalDepth(background.zDomain,33);
            transform = basis.discreteTransform(grid);
            input = exp(grid.z/170)+1i*sin(grid.z/213);
            for variable = ["F","G"]
                projected = basis.projectionOnGrid(grid.z,grid.weights,variable=variable);
                physical = transform.physicalProjection(variable=variable);
                testCase.verifyEqual(transform.targetGramMatrix(variable=variable),projected.targetGramMatrix,RelTol=1e-12,AbsTol=1e-10);
                testCase.verifyEqual(transform.targetMajorantGramMatrix(variable=variable),projected.majorantGramMatrix,RelTol=1e-12,AbsTol=1e-10);
                testCase.verifyEqual(physical.projectPhysical(input),projected.projectPhysical(input),RelTol=1e-12,AbsTol=1e-10);
            end
            testCase.verifyEqual(transform.activeModeMask(variable="G"),[false true(1,4)]);
            physicalG = transform.physicalProjection(variable="G");
            testCase.verifyEqual(physicalG.physicalDualMatrix(1,:),zeros(1,grid.nPoints));
        end

        function weakStratificationRetainsTheAnalyticalBand(testCase)
            for N0 = [1e-6 1e-7]
                background = IMBackground.constantStratification(N0=N0,zDomain=[-1000 0],surfaceCondition="rigidLid");
                for operator = ["truncatedCollocation","galerkin"]
                    [T,~,A] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,33,nInternalModes=31,forwardOperator=operator);
                    testCase.verifyEqual(A.retainedInternalModeCount,31)
                    testCase.verifyEqual(A.final.rankF,32)
                    testCase.verifyEqual(A.final.rankG,31)
                    testCase.verifyLessThan(T.roundTripError(variable="F"),1e-10)
                    testCase.verifyLessThan(T.roundTripError(variable="G"),1e-10)
                end
            end
        end

        function rootCountFailuresRefineWithCoordinateFallback(testCase)
            background = IMBackground.exponentialStratification(N0=5e-3,b=50,zDomain=[-1000 0],surfaceCondition="rigidLid");
            [T,grid,A] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,17);
            testCase.verifyEqual(A.retainedInternalModeCount,15)
            testCase.verifyEqual(A.convergence.coordinateKind,"z")
            testCase.verifyTrue(any(A.convergence.attempts.status=="underresolved"))
            testCase.verifyEqual(grid.nPoints,17)
            testCase.verifyLessThan(T.roundTripError(variable="F"),1e-9)
            testCase.verifyError(@() IMInternalModesDiscreteTransform.hydrostaticRoot(background,17,solver=IMSolverSpectral(nEVP=128,coordinateKind="wkb")), ...
                "IMInternalModesDiscreteTransform:UnconvergedHydrostaticRootModes")
        end

        function barotropicOnlyIsAnExplicitValidBand(testCase)
            background = IMBackground.constantStratification(N0=5e-3,zDomain=[-1000 0],surfaceCondition="rigidLid");
            for operator = ["truncatedCollocation","galerkin"]
                [T,grid,A,basis] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,17,nInternalModes=0,forwardOperator=operator);
                testCase.verifyEqual(A.retainedInternalModeCount,0)
                testCase.verifyEqual(basis.columnLabels,"barotropic:0")
                testCase.verifyEqual(T.activeModeMask(variable="G"),false)
                testCase.verifyEqual(T.forwardMatrix(variable="G"),zeros(1,17))
                testCase.verifyEqual(T.transformForward(ones(17,1),variable="F"),1,AbsTol=1e-12)
                testCase.verifyEqual(sum(grid.weights),1000,AbsTol=1e-9)
                testCase.verifyEqual(A.final.energyErrorG,0)
            end
        end

        function invalidRequestsAreRejected(testCase)
            free = IMBackground.constantStratification(N0=5e-3,zDomain=[-1000 0],surfaceCondition="freeSurface");
            testCase.verifyError(@() IMInternalModesDiscreteTransform.hydrostaticRoot(free,17),"IMInternalModesDiscreteTransform:HydrostaticRootRequiresRigidLid")
            neutral = IMBackground(N2=@(z) 1e-5*(z < -200),zDomain=[-1000 0],surfaceCondition="rigidLid");
            testCase.verifyError(@() IMInternalModesDiscreteTransform.hydrostaticRoot(neutral,17),"IMInternalModesDiscreteTransform:HydrostaticRootRequiresPositiveStratification")
            rigid = IMBackground.constantStratification(N0=5e-3,zDomain=[-1000 0],surfaceCondition="rigidLid");
            testCase.verifyError(@() IMInternalModesDiscreteTransform.hydrostaticRoot(rigid,17,solver=IMSolverFiniteDifference(z=linspace(-1000,0,129).')),"IMInternalModesDiscreteTransform:HydrostaticRootRequiresSpectralSolver")
        end

        function diagonalMetricForwardUsesEquivalentWeightedQR(testCase)
            z = linspace(-1,0,40).';
            Phi = [ones(40,1) cos(pi*(z+1)*(1:9))];
            w = 0.5+rand(RandStream("twister",Seed=3),40,1); w = w/sum(w);
            projection = IMProjection(Phi,diag(w),diag([1 0.5*ones(1,9)]));
            normal = (Phi.'*(w.*Phi))\(Phi.'.*w.');
            testCase.verifyEqual(projection.recovery.forwardMatrix,normal,AbsTol=1e-10*norm(normal))
            testCase.verifyEqual(projection.recovery.forwardMatrix*Phi,eye(10),AbsTol=1e-12)
        end
    end
end
