classdef IMScalarPhysicalGramTests < matlab.unittest.TestCase
    % Compare scalar physical forms with independently resolved reference rules.

    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function highestChebyshevProfileHasThePhysicalUnitDual(testCase)
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1);
            basis = polynomialBasis(evp,23);
            [z,w] = independentGaussRule(basis.zDomain,96);
            values = basis.u(z);
            projection = basis.projectionOnGrid(z,w,variable="u");
            testCase.verifyEqual(values.'*(w.*values),1,AbsTol=2e-13);
            testCase.verifyEqual(basis.gramMatrix(),1,AbsTol=2e-13);
            testCase.verifyEqual(projection.targetGramMatrix,1,AbsTol=2e-13);
            testCase.verifyEqual(projection.majorantGramMatrix,1,AbsTol=2e-13);
            testCase.verifyEqual(projection.projectPhysical(values),1,AbsTol=3e-13);
            testCase.verifyEmpty(projection.recovery);
        end

        function nearCapacityModesRetainTheFullRepresentedGram(testCase)
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1);
            basis = IMSolverSpectral(nEVP=24).solveEVP(evp,nModes=22);
            [z,w] = independentGaussRule(basis.zDomain,96);
            values = basis.u(z);
            reference = values.'*(w.*values);
            projection = basis.projectionOnGrid(z,w,variable="u");
            testCase.verifyGreaterThan(norm(reference-diag(diag(reference)),2),.1);
            testCase.verifyEqual(projection.targetGramMatrix,reference,AbsTol=8e-13);
            testCase.verifyEqual(projection.majorantGramMatrix,reference,AbsTol=8e-13);
            testCase.verifyEqual(projection.projectPhysical(values),eye(22),AbsTol=2e-12);
            transform = basis.discreteTransform(z=z,weights=w,nModes=22,gramTolerance=1e-10);
            testCase.verifyEqual(transform.targetGramMatrix,reference,AbsTol=8e-13);
            testCase.verifyEqual(transform.physicalDualMatrix,projection.physicalDualMatrix,AbsTol=2e-12);
            input = exp(.7*z)+1i*cos(2.3*z);
            expected = reference\(values.'*(w.*input));
            testCase.verifyEqual(projection.projectPhysical(input),expected,AbsTol=2e-12);
            % A changed sampled rule must remain visible in the physical map.
            changed = basis.projectionOnGrid(z,1.1*w,variable="u");
            testCase.verifyEqual(changed.targetGramMatrix,projection.targetGramMatrix);
            testCase.verifyEqual(changed.projectPhysical(input),1.1*expected,AbsTol=2e-12);
        end

        function signedEndpointsAndPartialBoundsKeepTheirContributions(testCase)
            surface = IMBoundaryCondition(a=0,b=1,c=-.4,d=0);
            bottom = IMBoundaryCondition(a=0,b=1,c=-.2,d=0);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=@(z) exp(.4*z),surfaceBoundary=surface,bottomBoundary=bottom);
            basis = polynomialBasis(evp,[0 1 23]);
            [z,w] = endpointRule(basis.zDomain,96);
            values = basis.u(z);
            interior = values.'*((w.*exp(.4*z)).*values);
            surfaceTerm = values(end,:).'*values(end,:);
            bottomTerm = values(1,:).'*values(1,:);
            reference = interior-.4*surfaceTerm+.2*bottomTerm;
            referenceMajorant = interior+.4*surfaceTerm+.2*bottomTerm;
            projection = basis.projectionOnGrid(z,w,variable="u");
            testCase.verifyLessThan(min(eig(reference)),0);
            testCase.verifyEqual(basis.gramMatrix(),reference,AbsTol=8e-13);
            testCase.verifyEqual(projection.targetGramMatrix,reference,AbsTol=8e-13);
            testCase.verifyEqual(projection.majorantGramMatrix,referenceMajorant,AbsTol=8e-13);
            testCase.verifyEqual(projection.projectPhysical(values),eye(3),AbsTol=2e-12);
            testCase.verifyEqual(basis.innerProductNormFactors(),sqrt(abs(diag(reference))).'.*basis.normalizationFactors(),AbsTol=6e-13);
            for bounds = {[-.8 -.2],[-1 -.2],[-.8 0]}
                interval = bounds{1};
                [zi,wi] = independentGaussRule(interval,96);
                sampled = basis.u(zi);
                expected = sampled.'*((wi.*exp(.4*zi)).*sampled);
                if interval(1) == -1, expected = expected+.2*bottomTerm; end
                if interval(2) == 0, expected = expected-.4*surfaceTerm; end
                testCase.verifyEqual(basis.gramMatrix(zBounds=interval),expected,AbsTol=8e-13);
            end
        end

        function normalizationAndSelectionPreservePhysicalForms(testCase)
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=@(z) exp(.4*z));
            basis = polynomialBasis(evp,[0 2 23]);
            scaled = basis.addNormalization("customScale",@(~,j) j+1).withNormalization("customScale");
            S = diag(basis.normalizationFactors()./scaled.normalizationFactors());
            [z,w] = independentGaussRule(basis.zDomain,96);
            physical = basis.projectionOnGrid(z,w,variable="u");
            changed = scaled.projectionOnGrid(z,w,variable="u");
            testCase.verifyEqual(changed.targetGramMatrix,S.'*physical.targetGramMatrix*S,AbsTol=3e-13);
            testCase.verifyEqual(changed.majorantGramMatrix,S.'*physical.majorantGramMatrix*S,AbsTol=3e-13);
            testCase.verifyEqual(changed.physicalDualMatrix,S\physical.physicalDualMatrix,AbsTol=3e-13);
            columns = [3 1];
            selected = scaled.projectionOnGrid(z,w,variable="u",columns=columns);
            testCase.verifyEqual(selected.targetGramMatrix,changed.targetGramMatrix(columns,columns),AbsTol=3e-13);
            values = scaled.u(z);
            reference = values.'*((w.*exp(.4*z)).*values);
            testCase.verifyEqual(selected.projectPhysical(values),reference(columns,columns)\reference(columns,:),AbsTol=6e-13);
            empty = scaled.projectionOnGrid(z,w,variable="u",columns=[]);
            testCase.verifySize(empty.targetGramMatrix,[0 0]);
            testCase.verifySize(empty.physicalDualMatrix,[0 numel(z)]);
        end

        function quadraticAssessmentUsesTheFullRepresentedGram(testCase)
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1);
            basis = polynomialBasis(evp,[0 2 23]);
            [z,w] = independentGaussRule(basis.zDomain,96);
            transform = basis.discreteTransform(z=z,weights=w,nModes=3,gramTolerance=1e-10,quadraticAliasingTolerance=1e-10);
            testCase.verifyTrue(transform.success);
            testCase.verifyEqual(transform.targetGramMatrix,basis.gramMatrix(),AbsTol=3e-13);
        end
    end
end

function basis = polynomialBasis(evp,degrees)
solver = IMSolverSpectral(nEVP=max(degrees)+1).configuredForEVP(evp);
native = zeros(solver.nEVP,numel(degrees));
for j = 1:numel(degrees), native(degrees(j)+1,j) = 1; end
basis = IMBasisSet(solver=solver,evp=evp,nativeModes=native,eigenvalues=1:numel(degrees),modeNumber=1:numel(degrees));
end

function [z,w] = independentGaussRule(domain,n)
% Golub-Welsch nodes provide a rule distinct from the production legpts call.
j = (1:n-1).';
offDiagonal = j./sqrt(4*j.^2-1);
[vectors,values] = eig(diag(offDiagonal,1)+diag(offDiagonal,-1));
[nodes,index] = sort(diag(values));
z = mean(domain)+diff(domain)*nodes/2;
w = diff(domain)*vectors(1,index).'.^2;
end

function [z,w] = endpointRule(domain,n)
[z,w] = independentGaussRule(domain,n);
z = [domain(1);z;domain(2)];
w = [0;w;0];
end
