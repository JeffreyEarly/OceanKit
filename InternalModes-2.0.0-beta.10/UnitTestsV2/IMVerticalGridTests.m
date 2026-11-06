classdef IMVerticalGridTests < matlab.unittest.TestCase
    % Regression and independent checks for IMVerticalGrid sampling rules.

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
        function physicalDepthIsClenshawCurtis(testCase)
            grid = IMVerticalGrid.physicalDepth([-1000 0],17);
            testCase.verifyEqual(grid.kind,"physicalDepth")
            testCase.verifyEqual(grid.nPoints,17)
            testCase.verifyEqual(grid.z([1 end]).',[-1000 0])
            % CC with 17 points integrates polynomials of degree 16 exactly.
            s = (grid.z+500)/500;
            for k = 0:16
                exact = 500*(1-(-1)^(k+1))/(k+1);
                testCase.verifyEqual(grid.integrate(s.^k),exact,AbsTol=1e-10*500)
            end
        end

        function constantStratificationReducesToPhysicalDepth(testCase)
            background = IMBackground.constantStratification(N0=5e-3,zDomain=[-2000 0],surfaceCondition="freeSurface");
            smoothed = IMVerticalGrid.smoothedBandwidth(background,33);
            depth = IMVerticalGrid.physicalDepth(background,33);
            testCase.verifyEqual(smoothed.z,depth.z,AbsTol=1e-9)
            testCase.verifyEqual(smoothed.weights,depth.weights,RelTol=1e-9)
        end

        function defaultSmoothingLengthScalesWithGridSpacing(testCase)
            background = IMVerticalGridTests.neutralBackground();
            for nPoints = [33 65]
                grid = IMVerticalGrid.smoothedBandwidth(background,nPoints,breakpoints=[-400 -200]);
                testCase.verifyEqual(grid.design.smoothingLength,2.5*1000/(nPoints-1),RelTol=1e-14)
                testCase.verifyEqual(grid.design.smoothingLengthInGridSpacings,2.5,RelTol=1e-14)
            end
        end

        function reproducesCoordinateGridStudyNodes(testCase)
            fixture = readtable(fullfile(fileparts(mfilename("fullpath")),"Fixtures","IMVerticalGridNeutralSmoothedNodes.csv"));
            background = IMVerticalGridTests.neutralBackground();
            cases = unique(fixture(:,["Nz","smoothingLength"]),"rows");
            for j = 1:height(cases)
                rows = fixture(fixture.Nz==cases.Nz(j) & fixture.smoothingLength==cases.smoothingLength(j),:);
                rows = sortrows(rows,"node");
                grid = IMVerticalGrid.smoothedBandwidth(background,cases.Nz(j), ...
                    smoothingLength=cases.smoothingLength(j),breakpoints=[-400 -200]);
                label = sprintf("Nz=%d, Ls=%g",cases.Nz(j),cases.smoothingLength(j));
                testCase.verifyEqual(grid.z,rows.z,"Nodes differ for "+label,AbsTol=1e-8)
                testCase.verifyEqual(grid.weights,rows.weight,"Weights differ for "+label,RelTol=1e-8)
            end
        end

        function smoothedGridIntegratesSmoothFunctionsAccurately(testCase)
            N2 = @(z) 1e-6+5e-5*exp(z/300);
            background = IMBackground(N2=N2,zDomain=[-4000 0],latitude=33,surfaceCondition="freeSurface");
            grid = IMVerticalGrid.smoothedBandwidth(background,65);
            testCase.verifyEqual(sum(grid.weights),4000,RelTol=1e-12)
            testCase.verifyGreaterThan(min(grid.weights),0)
            exact = integral(N2,-4000,0,AbsTol=1e-14,RelTol=1e-13);
            testCase.verifyEqual(grid.integrate(N2(grid.z)),exact,RelTol=1e-6)
            % Points concentrate where stratification is strong.
            testCase.verifyLessThan(grid.z(end)-grid.z(end-1),grid.z(2)-grid.z(1))
        end

        function smoothingLengthLimitsBracketTheRule(testCase)
            background = IMVerticalGridTests.neutralBackground();
            wide = IMVerticalGrid.smoothedBandwidth(background,33,smoothingLength=1e5,breakpoints=[-400 -200]);
            depth = IMVerticalGrid.physicalDepth(background,33);
            testCase.verifyEqual(wide.z,depth.z,AbsTol=1)
        end

        function modeRootGridHasRequestedCountAndPositiveWeights(testCase)
            solver = IMSolverSpectral(nEVP=128);
            N2 = @(z) 1e-4*exp(2*z/400);
            evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-1000 0],latitude=33,surfaceBoundary=IMBoundaryCondition.neumann());
            basisSet = solver.solveEVP(evp,nModes=40);
            grid = IMVerticalGrid.modeRoots(basisSet,25);
            testCase.verifyEqual(grid.kind,"modeRoots")
            testCase.verifyEqual(grid.nPoints,25)
            testCase.verifyEqual(grid.z([1 end]).',[-1000 0],AbsTol=1e-12)
            testCase.verifyGreaterThanOrEqual(min(grid.weights),0)
            testCase.verifyEqual(sum(grid.weights),1000,RelTol=1e-12)
        end

        function explicitGridValidatesItsRule(testCase)
            grid = IMVerticalGrid(z=[-1;-0.5;0],weights=[0.25;0.5;0.25]);
            testCase.verifyEqual(grid.zDomain,[-1 0])
            testCase.verifyEqual(grid.kind,"custom")
            testCase.verifyError(@() IMVerticalGrid(z=[0;-1],weights=[0.5;0.5]),"IMVerticalGrid:InvalidPoints")
            testCase.verifyError(@() IMVerticalGrid(z=[-1;0],weights=[1.5;-0.5]),"IMVerticalGrid:InvalidWeights")
            testCase.verifyError(@() IMVerticalGrid(z=[-1;0],weights=[0.4;0.4]),"IMVerticalGrid:InvalidWeights")
            testCase.verifyError(@() IMVerticalGrid(z=[-.75;-.25],weights=[.5;.5],zDomain=[-1 0]),"IMVerticalGrid:MissingEndpoint")
            testCase.verifyError(@() IMVerticalGrid(z=[-1;-.25],weights=[.5;.25],zDomain=[-1 0]),"IMVerticalGrid:MissingEndpoint")
            background = IMVerticalGridTests.neutralBackground();
            testCase.verifyError(@() IMVerticalGrid.smoothedBandwidth(background,33,smoothingLength=-1),"IMVerticalGrid:InvalidSmoothingLength")
        end

        function nearEndpointsAreStoredExactly(testCase)
            % Downstream endpoint traces match zDomain to machine precision.
            grid = IMVerticalGrid(z=[-1000+1e-8;-500;-1e-8],weights=[250;500;250],zDomain=[-1000 0]);
            testCase.verifyTrue(grid.z(1) == -1000)
            testCase.verifyTrue(grid.z(end) == 0)
            testCase.verifyEqual(grid.z(2),-500)
            endpointTolerance = 100*eps(max(1,max(abs(grid.zDomain))));
            testCase.verifyNotEmpty(find(abs(grid.z-grid.zDomain(1)) <= endpointTolerance,1))
            testCase.verifyNotEmpty(find(abs(grid.z-grid.zDomain(2)) <= endpointTolerance,1))
            % Every constructor returns exact endpoints.
            background = IMVerticalGridTests.neutralBackground();
            grids = {IMVerticalGrid.physicalDepth(background,33), ...
                IMVerticalGrid.smoothedBandwidth(background,33,breakpoints=[-400 -200])};
            for j = 1:numel(grids)
                testCase.verifyTrue(isequal(grids{j}.z([1 end]).',background.zDomain))
            end
        end

        function stratificationJumpWithBreakpointIsAccepted(testCase)
            % A step in N2 makes Q_zz jump; residuals must be one-sided there.
            background = IMBackground(N2=@(z) 2.5e-5*(z < -200),zDomain=[-1000 0],surfaceCondition="freeSurface");
            grid = IMVerticalGrid.smoothedBandwidth(background,65,breakpoints=-200);
            testCase.verifyEqual(grid.nPoints,65)
            testCase.verifyLessThan(grid.design.filterRelativeResidual,1e-9)
            testCase.verifyGreaterThan(min(grid.weights),0)
            testCase.verifyEqual(sum(grid.weights),1000,RelTol=1e-12)
            testCase.verifyLessThan(grid.design.mapResidual,1e-10)
        end
    end

    methods (Static, Access = private)
        function background = neutralBackground()
            % Neutral 200 m surface layer over a smooth quintic transition.
            N2 = @(z) 2.5e-5*IMVerticalGridTests.ramp(min(max((-z-200)/200,0),1));
            background = IMBackground(N2=N2,zDomain=[-1000 0],latitude=45,rotationRate=7.2921e-5,surfaceCondition="freeSurface");
        end

        function value = ramp(t)
            value = 6*t.^5-15*t.^4+10*t.^3;
        end
    end
end
