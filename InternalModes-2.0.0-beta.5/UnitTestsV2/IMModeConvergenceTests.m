classdef IMModeConvergenceTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addPackage(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename("fullpath")))));
        end
    end
    methods (Test)
        function singleModeStillReturnsColumnOrientedMeasurements(testCase)
            [a,z,w] = preparedModes();
            a = reorder(a,1);
            result = assessModeConvergence(a,a,z,w);
            testCase.verifySize(result.measurements,[8 5]);
            testCase.verifyEqual(result.measurements.columnLabel,repmat("0",8,1));
            testCase.verifyEqual(result.measurements.value,zeros(8,1));
        end

        function identicalModesRetainIdentityAndProvenance(testCase)
            [a,z,w] = preparedModes();
            original = a;
            result = assessModeConvergence(a,a,z,w);
            testCase.verifyEqual(result.measurements.value,zeros(height(result.measurements),1));
            testCase.verifyEqual(result.matches.columnLabel,["0";"1"]);
            testCase.verifyEqual(result.matches.referenceColumn,[1;2]);
            testCase.verifyEqual(result.provenance.candidate,a.provenance);
            testCase.verifyFalse(result.coverage.absoluteAccuracyGuarantee);
            testCase.verifyEqual(result.coverage.referenceAccuracy,"unverified");
            testCase.verifyEqual(a,original);
        end

        function reorderedReferenceAndCommonSignsPreservePhysicalModes(testCase)
            [a,z,w] = preparedModes();
            b = reorder(a,[2 1]);
            for field = ["F","G"]
                b.values.(field) = -b.values.(field);
                b.derivatives.(field) = -b.derivatives.(field);
            end
            result = assessModeConvergence(a,b,z,w);
            testCase.verifyEqual(result.matches.referenceColumn,[2;1]);
            testCase.verifyEqual(result.matches.orientation,[-1;-1]);
            testCase.verifyEqual(result.measurements.value,zeros(height(result.measurements),1));
        end

        function inconsistentPolarizationCannotBeSignAlignedAway(testCase)
            [a,z,w] = preparedModes();
            b = a;
            b.values.G = -b.values.G;
            b.derivatives.G = -b.derivatives.G;
            result = assessModeConvergence(a,b,z,w);
            shapes = result.measurements(result.measurements.quantity == "shape",:);
            testCase.verifyEqual(max(shapes.value),2,AbsTol=1e-14);
        end

        function missingAndAmbiguousModeLabelsAreInconclusive(testCase)
            [a,z,w] = preparedModes();
            b = reorder(a,1);
            result = assessModeConvergence(a,b,z,w);
            testCase.verifyEqual(result.matches.status,["matched";"inconclusive"]);
            testCase.verifyTrue(all(result.measurements.status(result.measurements.columnLabel=="1") == "inconclusive"));
            b = a;
            b.identity.columnLabels = ["0","0"];
            result = assessModeConvergence(a,b,z,w);
            testCase.verifyTrue(all(result.matches.status == "inconclusive"));
            result = assessModeConvergence(b,a,z,w);
            testCase.verifyTrue(all(result.matches.status == "inconclusive"));
        end

        function physicalIdentityAndGridMustAgree(testCase)
            [a,z,w] = preparedModes();
            b = a;
            b.identity.kappa = a.identity.kappa+eps(a.identity.kappa);
            testCase.verifyError(@() assessModeConvergence(a,b,z,w),"assessModeConvergence:IdentityMismatch");
            testCase.verifyError(@() assessModeConvergence(a,a,flipud(z),w),"assessModeConvergence:InvalidQuadrature");
            b = a;
            b.values.F = b.values.F(1:end-1,:);
            testCase.verifyError(@() assessModeConvergence(a,b,z,w),"assessModeConvergence:InvalidSamples");
            b = rmfield(a,"provenance");
            testCase.verifyError(@() assessModeConvergence(a,b,z,w),"assessModeConvergence:InvalidPreparedModes");
        end

        function h1KeepsConstantAndTinyDerivativesMeaningful(testCase)
            [a,z,w] = preparedModes();
            b = a;
            a.derivatives.F(:,1) = 1e-13;
            result = assessModeConvergence(a,b,z,w);
            row = result.measurements.columnLabel == "0" & result.measurements.variable == "F";
            testCase.verifyEqual(result.measurements.value(row & result.measurements.quantity=="derivative"),Inf);
            testCase.verifyEqual(result.measurements.value(row & result.measurements.quantity=="h1"),1e-13,RelTol=1e-13);
            b.derivatives.F(:,1) = 1e-15;
            a.derivatives.F(:,1) = 2e-15;
            result = assessModeConvergence(a,b,z,w);
            testCase.verifyEqual(result.measurements.value(row & result.measurements.quantity=="derivative"),1,AbsTol=1e-14);
        end

        function physicallyNullFieldsRequireNullDiscrepancy(testCase)
            [a,z,w] = preparedModes();
            a.values.G(:,1) = 0;
            a.derivatives.G(:,1) = 0;
            b = a;
            result = assessModeConvergence(a,b,z,w);
            row = result.measurements.columnLabel == "0" & result.measurements.variable == "G";
            testCase.verifyEqual(result.measurements.value(row),[0;0;0]);
            a.values.G(:,1) = eps;
            result = assessModeConvergence(a,b,z,w);
            testCase.verifyEqual(result.measurements.value(row),[Inf;0;Inf]);
        end

        function exceptionalScalarsAreExplicit(testCase)
            [a,z,w] = preparedModes();
            a.equivalentDepths = [Inf,0];
            b = a;
            result = assessModeConvergence(a,b,z,w);
            rows = result.measurements.quantity == "equivalentDepth";
            testCase.verifyEqual(result.measurements.value(rows),[0;0]);
            a.equivalentDepths = [1,1];
            result = assessModeConvergence(a,b,z,w);
            testCase.verifyEqual(result.measurements.status(rows),["inconclusive";"measured"]);
            testCase.verifyEqual(result.measurements.value(find(rows,1,"last")),Inf);
            a.eigenvalues(1) = NaN;
            result = assessModeConvergence(a,b,z,w);
            testCase.verifyEqual(result.measurements.status(1),"inconclusive");
        end

        function omittedMeasurementsStayUnrequested(testCase)
            [a,z,w] = preparedModes();
            a = rmfield(a,["derivatives","eigenvalues","equivalentDepths"]);
            result = assessModeConvergence(a,a,z,w);
            testCase.verifyTrue(all(result.measurements.status(result.measurements.quantity~="shape") == "notRequested"));
            testCase.verifyTrue(all(result.measurements.status(result.measurements.quantity=="shape") == "measured"));
        end

        function upperModeErrorSupportsAnExplicitContiguousPrefix(testCase)
            [a,z,w] = preparedModes();
            b = a;
            a.values.F(:,2) = 1.1*a.values.F(:,2);
            a.derivatives.F(:,2) = 1.1*a.derivatives.F(:,2);
            result = assessModeConvergence(a,b,z,w);
            rows = result.measurements.quantity=="h1" & result.measurements.variable=="F";
            errors = result.measurements.value(rows);
            testCase.verifyEqual(errors,[0;.1],AbsTol=1e-14);
            testCase.verifyEqual(sum(cumprod(errors<.01)),1);
            testCase.verifyEqual(result.coverage.requestedColumnCount,2);
        end

        function spectralRefinementComparesTheActualResolvedWaveModes(testCase)
            D = 1000;
            kappa = 2*pi/1e4;
            for profile = ["constant","exponential"]
                if profile == "constant"
                    N2 = @(z) 1e-4*ones(size(z));
                else
                    N2 = @(z) 1e-4*exp(2*z/700);
                end
                evp = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=[-D 0],k=kappa,f0=1e-4,surfaceBoundary=IMBoundaryCondition(a=0,b=1,c=1,d=0));
                basis = IMSolverSpectral(nEVP=64,coordinateKind="wkb").solveEVP(evp,nModes=6);
                refined = IMSolverSpectral(nEVP=128,coordinateKind="wkb").solveEVP(evp,nModes=6);
                rule = IMSolverSpectral(nEVP=257).configuredForEVP(evp);
                [z,w] = rule.nativeQuadratureRule([-D 0]);
                a = waveSamples(basis,z,kappa,64);
                b = waveSamples(refined,z,kappa,128);
                result = assessModeConvergence(a,b,z,w);
                rows = ismember(result.measurements.quantity,["equivalentDepth","h1"]);
                testCase.verifyTrue(all(result.measurements.status(rows)=="measured"));
                testCase.verifyLessThan(max(result.measurements.value(rows)),1e-7);
                if profile == "constant"
                    exact = IMConstantStratificationSolution(N0=.01,zDomain=[-D 0],f0=1e-4).internalModes(evp,nModes=6);
                    b = rmfield(b,"derivatives");
                    b.values = struct(F=exact.F(z),G=exact.G(z));
                    b.equivalentDepths = exact.h;
                    b.eigenvalues = exact.eigenvalues;
                    b.provenance = struct(solution="constantStratification");
                    result = assessModeConvergence(a,b,z,w);
                    rows = ismember(result.measurements.quantity,["equivalentDepth","shape"]);
                    testCase.verifyLessThan(max(result.measurements.value(rows)),1e-8);
                end
            end
        end
    end
end

function [a,z,w] = preparedModes()
z = linspace(-1,0,33).';
w = [0.5;ones(31,1);0.5]/32;
identity = struct(family="testWaves",columnLabels=["0","1"],normalization="fixed",zDomain=[-1 0],kappa=1e-4);
values = struct(F=[ones(size(z)),cos(pi*z)],G=[ones(size(z)),sin(pi*z)]);
derivatives = struct(F=[zeros(size(z)),-pi*sin(pi*z)],G=[zeros(size(z)),pi*cos(pi*z)]);
a = struct(identity=identity,values=values,derivatives=derivatives,eigenvalues=[0 1],equivalentDepths=[Inf 1],provenance=struct(nEVP=32,solverClass="preparedTest"));
end

function b = reorder(a,order)
b = a;
b.identity.columnLabels = a.identity.columnLabels(order);
b.eigenvalues = a.eigenvalues(order);
b.equivalentDepths = a.equivalentDepths(order);
for field = ["F","G"]
    b.values.(field) = a.values.(field)(:,order);
    b.derivatives.(field) = a.derivatives.(field)(:,order);
end
end

function samples = waveSamples(basis,z,kappa,nEVP)
factors = basis.normalizationFactors(basis.normalization);
values = struct(F=basis.F(z),G=basis.G(z));
derivatives = struct(F=basis.h.*basis.solver.evaluatePhysicalDerivative(basis.nativeModes,z,2)./factors,G=basis.solver.evaluatePhysicalDerivative(basis.nativeModes,z,1)./factors);
identity = struct(family=string(basis.evp.name),columnLabels=string(basis.modeNumber),normalization=string(basis.normalization),zDomain=basis.zDomain,kappa=kappa);
samples = struct(identity=identity,values=values,derivatives=derivatives,eigenvalues=basis.eigenvalues,equivalentDepths=basis.h,provenance=struct(nEVP=nEVP,coordinateKind="wkb",solverClass="IMSolverSpectral"));
end
