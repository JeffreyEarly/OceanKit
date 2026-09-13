classdef TestTotalFlowComponent < matlab.unittest.TestCase
    properties (TestParameter)
        transformType = struct( ...
            constantHydrostatic="constantHydrostatic", ...
            constantNonhydrostatic="constantNonhydrostatic", ...
            hydrostatic="hydrostatic", ...
            boussinesq="boussinesq", ...
            stratifiedQG="stratifiedQG", ...
            barotropicQG="barotropicQG")
    end

    methods (Test, TestTags = "full")
        function totalCountMatchesSortedComponentRanges(testCase,transformType)
            wvt = testCase.newTransform(transformType);
            [components,nModesByComponent] = testCase.sortedComponents(wvt);

            testCase.verifyEqual(wvt.totalFlowComponent.nModes,sum(nModesByComponent))
            if wvt.hasWaveComponent
                testCase.verifyEqual(string({components.shortName}),["geostrophic" "inertial" "mda" "wave"])
            else
                testCase.verifyEqual(string({components.shortName}),"geostrophic")
            end
        end

        function componentBoundariesDelegateInCallerOrder(testCase,transformType)
            wvt = testCase.newTransform(transformType);
            indices = testCase.componentBoundaryIndices(wvt);

            seedRandomNumberGenerator(testCase,2901);
            seededState = rng;
            actual = wvt.totalFlowComponent.solutionForModeAtIndex(indices,amplitude='random');
            rng(seededState)
            expected = testCase.directSolutions(wvt,indices,'random');

            testCase.verifySize(actual,size(indices))
            testCase.verifyEquivalentSolutions(actual,expected,wvt)
        end

        function currentCoefficientsAndOutputShapesArePreserved(testCase,transformType)
            wvt = testCase.newTransform(transformType);
            seedRandomNumberGenerator(testCase,2902);
            wvt.initWithRandomFlow(uvMax=0.01);

            indices = flipud(testCase.componentBoundaryIndices(wvt));
            actual = wvt.totalFlowComponent.solutionForModeAtIndex(indices,amplitude='wvt');
            expected = testCase.directSolutions(wvt,indices,'wvt');

            testCase.verifySize(actual,size(indices))
            testCase.verifyEquivalentSolutions(actual,expected,wvt)
            for iSolution = 1:numel(actual)
                testCase.verifySolutionMatchesCurrentCoefficient(actual(iSolution),wvt)
            end

            scalarSolution = wvt.totalFlowComponent.solutionForModeAtIndex(indices(1),amplitude='wvt');
            testCase.verifySize(scalarSolution,[1 1])
            testCase.verifySolutionMatchesCurrentCoefficient(scalarSolution,wvt)
        end

        function invalidIndicesUseBuiltInValidation(testCase,transformType)
            wvt = testCase.newTransform(transformType);
            lookup = @(index) wvt.totalFlowComponent.solutionForModeAtIndex(index);

            testCase.verifyError(@() lookup(0),'MATLAB:validators:mustBePositive')
            testCase.verifyError(@() lookup(-1),'MATLAB:validators:mustBePositive')
            testCase.verifyError(@() lookup(1.5),'MATLAB:validators:mustBeInteger')
            testCase.verifyError(@() lookup(wvt.totalFlowComponent.nModes+1),'MATLAB:validators:mustBeLessThanOrEqual')
        end

        function barotropicModeIndexWrappersRoundTrip(testCase)
            wvt = testCase.newTransform("barotropicQG");
            spectralIndices = find(wvt.geostrophicComponent.maskA0);

            for iIndex = 1:numel(spectralIndices)
                [kMode,lMode,jMode] = wvt.modeNumberFromIndex(spectralIndices(iIndex));
                testCase.verifyEqual(jMode,wvt.j)
                testCase.verifyEqual(wvt.indexFromModeNumber(kMode,lMode,jMode),spectralIndices(iIndex))
            end
            [kMode,lMode] = wvt.modeNumberFromIndex(spectralIndices(1));
            testCase.verifyError(@() wvt.indexFromModeNumber(kMode,lMode,1-wvt.j),'MATLAB:validators:mustBeMember')

            seedRandomNumberGenerator(testCase,2903);
            wvt.initWithRandomFlow(uvMax=0.01);
            localIndices = [1; wvt.geostrophicComponent.nModes];
            solutions = wvt.geostrophicComponent.solutionForModeAtIndex(localIndices,amplitude='wvt');
            testCase.verifySize(solutions,size(localIndices))
            for iSolution = 1:numel(solutions)
                testCase.verifySolutionMatchesCurrentCoefficient(solutions(iSolution),wvt)
            end
        end
    end

    methods (Access = private)
        function verifyEquivalentSolutions(testCase,actual,expected,wvt)
            testCase.verifyEqual(numel(actual),numel(expected))
            if isa(wvt,"WVTransformBarotropicQG")
                fieldNames = ["u" "v" "eta" "qgpv" "psi"];
            else
                fieldNames = ["u" "v" "w" "eta" "rho_e" "p" "qgpv" "psi"];
            end
            x = 0.17*wvt.Lx;
            y = 0.23*wvt.Ly;
            z = -0.37*wvt.Lz;
            t = 123;
            for iSolution = 1:numel(actual)
                testCase.verifyEqual(actual(iSolution).kMode,expected(iSolution).kMode)
                testCase.verifyEqual(actual(iSolution).lMode,expected(iSolution).lMode)
                testCase.verifyEqual(actual(iSolution).jMode,expected(iSolution).jMode)
                testCase.verifyEqual(actual(iSolution).amplitude,expected(iSolution).amplitude,RelTol=1e-12,AbsTol=1e-12)
                testCase.verifyEqual(actual(iSolution).phase,expected(iSolution).phase,RelTol=1e-12,AbsTol=1e-12)
                testCase.verifyEqual(actual(iSolution).coefficientMatrix,expected(iSolution).coefficientMatrix)
                testCase.verifyEqual(actual(iSolution).coefficientMatrixIndex,expected(iSolution).coefficientMatrixIndex)
                testCase.verifyEqual(actual(iSolution).coefficientMatrixAmplitude,expected(iSolution).coefficientMatrixAmplitude,RelTol=1e-12,AbsTol=1e-12)
                for fieldName = fieldNames
                    actualFunction = actual(iSolution).(fieldName);
                    expectedFunction = expected(iSolution).(fieldName);
                    testCase.verifyEqual(actualFunction(x,y,z,t),expectedFunction(x,y,z,t),RelTol=1e-12,AbsTol=1e-12)
                end
                if ~isa(wvt,"WVTransformBarotropicQG")
                    testCase.verifyEqual(actual(iSolution).ssh(x,y,t),expected(iSolution).ssh(x,y,t),RelTol=1e-12,AbsTol=1e-12)
                end
            end
        end

        function verifySolutionMatchesCurrentCoefficient(testCase,solution,wvt)
            switch solution.coefficientMatrix
                case WVCoefficientMatrix.Ap
                    coefficient = wvt.Ap(solution.coefficientMatrixIndex);
                case WVCoefficientMatrix.Am
                    coefficient = wvt.Am(solution.coefficientMatrixIndex);
                case WVCoefficientMatrix.A0
                    coefficient = wvt.A0(solution.coefficientMatrixIndex);
            end
            testCase.verifyEqual(solution.coefficientMatrixAmplitude,coefficient,RelTol=1e-12,AbsTol=1e-12)
        end
    end

    methods (Static, Access = private)
        function [components,nModesByComponent,lastIndexByComponent] = sortedComponents(wvt)
            components = wvt.primaryFlowComponents;
            [~,componentOrder] = sort(string({components.shortName}));
            components = components(componentOrder);
            nModesByComponent = arrayfun(@(component) component.nModes,components);
            lastIndexByComponent = cumsum(nModesByComponent);
        end

        function indices = componentBoundaryIndices(wvt)
            [~,~,lastIndexByComponent] = TestTotalFlowComponent.sortedComponents(wvt);
            firstIndexByComponent = [1 lastIndexByComponent(1:end-1)+1];
            indices = unique([firstIndexByComponent(:); lastIndexByComponent(:)]);
        end

        function solutions = directSolutions(wvt,indices,amplitude)
            [components,nModesByComponent,lastIndexByComponent] = TestTotalFlowComponent.sortedComponents(wvt);
            solutions = WVOrthogonalSolution.empty(length(indices),0);
            for iSolution = 1:length(indices)
                iComponent = find(indices(iSolution) <= lastIndexByComponent,1);
                firstIndexForComponent = lastIndexByComponent(iComponent) - nModesByComponent(iComponent) + 1;
                localIndex = indices(iSolution) - firstIndexForComponent + 1;
                solutions(iSolution) = components(iComponent).solutionForModeAtIndex(localIndex,amplitude=amplitude);
            end
        end

        function wvt = newTransform(transformType)
            Lxyz = [8e3 6e3 1e3];
            Nxyz = [8 6 9];
            N2 = @(z) 2e-5*exp(z/4000);
            switch transformType
                case "constantHydrostatic"
                    wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=true,shouldAntialias=false);
                case "constantNonhydrostatic"
                    wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=false);
                case "hydrostatic"
                    wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false);
                case "boussinesq"
                    wvt = WVTransformBoussinesq(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false);
                case "stratifiedQG"
                    wvt = WVTransformStratifiedQG(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false);
                case "barotropicQG"
                    wvt = WVTransformBarotropicQG(Lxyz(1:2),Nxyz(1:2),latitude=45,shouldAntialias=false);
            end
        end
    end
end
