classdef TestDegreesOfFreedomSummary < matlab.unittest.TestCase
    properties (TestParameter)
        transformType = struct( ...
            constantHydrostatic="constantHydrostatic", ...
            constantNonhydrostatic="constantNonhydrostatic", ...
            hydrostatic="hydrostatic", ...
            boussinesq="boussinesq", ...
            stratifiedQG="stratifiedQG", ...
            barotropicQG="barotropicQG")
        shouldAntialias = struct(disabled=false,enabled=true)
    end

    methods (Test, TestTags = "full")
        function summaryMatchesActivePrimaryMasks(testCase,transformType,shouldAntialias)
            wvt = testCase.newTransform(transformType,shouldAntialias);
            [expectedOutput,expectedData] = testCase.expectedOutput(wvt,transformType);

            lastwarn('')
            actualOutput = evalc('wvt.summarizeDegreesOfFreedom();');
            [warningMessage,warningID] = lastwarn;
            repeatedOutput = evalc('wvt.summarizeDegreesOfFreedom();');

            testCase.verifyEmpty(warningMessage)
            testCase.verifyEmpty(warningID)
            testCase.verifyEqual(actualOutput,repeatedOutput)
            testCase.verifyEqual(actualOutput,expectedOutput)

            components = testCase.sortedComponents(wvt);
            testCase.verifyEqual(string({components.name}),expectedData.names)
            testCase.verifyEqual(string({components.shortName}),expectedData.shortNames)
            testCase.verifyEqual(arrayfun(@(component) component.degreesOfFreedomPerMode,components),expectedData.degreesOfFreedomPerMode)
            testCase.verifyEqual(arrayfun(@(component) component.nModes,components),expectedData.nModes)

            obsoleteText = ["Not implemented" "Nx-1" "discarded" "aliased" "Spatial domain" "four major" "equivalent"];
            testCase.verifyFalse(any(contains(string(actualOutput),obsoleteText,IgnoreCase=true)))
        end

        function antialiasingChangesOnlyReportedSpectralCounts(testCase,transformType)
            fullTransform = testCase.newTransform(transformType,false);
            antialiasedTransform = testCase.newTransform(transformType,true);

            [expectedDimensionNames,expectedGridSize] = testCase.expectedGrid(transformType);
            testCase.verifyEqual(string(fullTransform.spatialDimensionNames()),expectedDimensionNames)
            testCase.verifyEqual(string(antialiasedTransform.spatialDimensionNames()),expectedDimensionNames)
            testCase.verifyEqual(fullTransform.spatialMatrixSize,expectedGridSize)
            testCase.verifyEqual(antialiasedTransform.spatialMatrixSize,expectedGridSize)

            fullTotal = testCase.maskDerivedTotal(fullTransform);
            antialiasedTotal = testCase.maskDerivedTotal(antialiasedTransform);
            testCase.verifyLessThan(antialiasedTotal,fullTotal)

            fullOutput = evalc('fullTransform.summarizeDegreesOfFreedom();');
            antialiasedOutput = evalc('antialiasedTransform.summarizeDegreesOfFreedom();');
            gridLine = "Spatial grid: " + strjoin(expectedDimensionNames + "=" + string(expectedGridSize),", ");
            testCase.verifyTrue(contains(string(fullOutput),gridLine))
            testCase.verifyTrue(contains(string(antialiasedOutput),gridLine))
            testCase.verifyNotEqual(fullOutput,antialiasedOutput)
        end
    end

    methods (Static, Access = private)
        function [output,data] = expectedOutput(wvt,transformType)
            components = TestDegreesOfFreedomSummary.sortedComponents(wvt);
            shortNames = string({components.shortName});
            names = strings(size(shortNames));
            degreesOfFreedomPerMode = zeros(size(shortNames));
            nModes = zeros(size(shortNames));
            for iComponent = 1:numel(components)
                [names(iComponent),degreesOfFreedomPerMode(iComponent)] = TestDegreesOfFreedomSummary.expectedComponentProperties(shortNames(iComponent));
                nModes(iComponent) = TestDegreesOfFreedomSummary.primaryModeCount(components(iComponent));
            end
            subtotals = nModes.*degreesOfFreedomPerMode;

            [dimensionNames,gridSize] = TestDegreesOfFreedomSummary.expectedGrid(transformType);
            gridDescription = strjoin(dimensionNames + "=" + string(gridSize),", ");
            componentWidth = max([strlength("Component") strlength(names)]);
            shortNameWidth = max([strlength("Short name") strlength(shortNames)]);

            lines = strings(numel(components)+8,1);
            lines(1) = "Degrees of freedom summary";
            lines(2) = "Spatial grid: " + gridDescription;
            lines(3) = "";
            lines(4) = "Primary flow components (spectral):";
            lines(5) = string(sprintf('%-*s  %-*s  %9s  %9s  %9s',componentWidth,'Component',shortNameWidth,'Short name','Modes','DOF/mode','Subtotal'));
            lines(6) = string(sprintf('%-*s  %-*s  %9s  %9s  %9s',componentWidth,repmat('-',1,strlength("Component")),shortNameWidth,repmat('-',1,strlength("Short name")),repmat('-',1,strlength("Modes")),repmat('-',1,strlength("DOF/mode")),repmat('-',1,strlength("Subtotal"))));
            for iComponent = 1:numel(components)
                lines(iComponent+6) = string(sprintf('%-*s  %-*s  %9d  %9d  %9d',componentWidth,names(iComponent),shortNameWidth,shortNames(iComponent),nModes(iComponent),degreesOfFreedomPerMode(iComponent),subtotals(iComponent)));
            end
            lines(end-1) = "";
            lines(end) = "Total active spectral degrees of freedom: " + sum(subtotals);
            output = char(join(lines,newline) + string(newline));

            data = struct(names=names,shortNames=shortNames,nModes=nModes,degreesOfFreedomPerMode=degreesOfFreedomPerMode);
        end

        function components = sortedComponents(wvt)
            components = wvt.primaryFlowComponents;
            [~,componentOrder] = sort(string({components.shortName}));
            components = components(componentOrder);
        end

        function count = primaryModeCount(component)
            count = nnz(component.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.Ap)) ...
                + nnz(component.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.Am)) ...
                + nnz(component.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.A0));
        end

        function total = maskDerivedTotal(wvt)
            components = TestDegreesOfFreedomSummary.sortedComponents(wvt);
            total = 0;
            for iComponent = 1:numel(components)
                [~,degreesOfFreedomPerMode] = TestDegreesOfFreedomSummary.expectedComponentProperties(string(components(iComponent).shortName));
                total = total + TestDegreesOfFreedomSummary.primaryModeCount(components(iComponent))*degreesOfFreedomPerMode;
            end
        end

        function [name,degreesOfFreedomPerMode] = expectedComponentProperties(shortName)
            switch shortName
                case "geostrophic"
                    name = "geostrophic";
                    degreesOfFreedomPerMode = 2;
                case "inertial"
                    name = "inertial oscillation";
                    degreesOfFreedomPerMode = 2;
                case "mda"
                    name = "mean density anomaly";
                    degreesOfFreedomPerMode = 1;
                case "wave"
                    name = "internal gravity wave";
                    degreesOfFreedomPerMode = 2;
                otherwise
                    error("TestDegreesOfFreedomSummary:UnknownComponent","Unexpected primary component '%s'.",shortName)
            end
        end

        function [dimensionNames,gridSize] = expectedGrid(transformType)
            if transformType == "barotropicQG"
                dimensionNames = ["x" "y"];
                gridSize = [8 6];
            else
                dimensionNames = ["x" "y" "z"];
                gridSize = [8 6 9];
            end
        end

        function wvt = newTransform(transformType,shouldAntialias)
            Lxyz = [8e3 6e3 1e3];
            Nxyz = [8 6 9];
            N2 = @(z) 2e-5*exp(z/4000);
            switch transformType
                case "constantHydrostatic"
                    wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=true,shouldAntialias=shouldAntialias);
                case "constantNonhydrostatic"
                    wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=shouldAntialias);
                case "hydrostatic"
                    wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=shouldAntialias);
                case "boussinesq"
                    wvt = WVTransformBoussinesq(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=shouldAntialias);
                case "stratifiedQG"
                    wvt = WVTransformStratifiedQG(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=shouldAntialias);
                case "barotropicQG"
                    wvt = WVTransformBarotropicQG(Lxyz(1:2),Nxyz(1:2),latitude=45,shouldAntialias=shouldAntialias);
            end
        end
    end
end
