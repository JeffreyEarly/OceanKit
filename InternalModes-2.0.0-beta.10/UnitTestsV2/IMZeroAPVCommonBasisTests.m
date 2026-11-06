classdef IMZeroAPVCommonBasisTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end

    methods (Test)
        function physicalFormsAndResponsesMatchIndependentQuadrature(testCase)
            for surface = ["freeSurface","rigidLid"]
                for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=64)}
                    boundary = makeEndpointBasis(solver{1},surface,1.2);
                    for current = {boundary,boundary.orthogonalizeGeneralizedEnergy(g0=-3,gd=.4).withNormalization("maxAbsF").selectColumns([2 1])}
                        b = current{1};
                        [H,B0,Bd,response] = independentForms(b);
                        testCase.verifyEqual(b.energyMatrix,H,RelTol=3e-10,AbsTol=3e-12);
                        testCase.verifyEqual(b.surfaceBuoyancyMatrix,B0,RelTol=3e-10,AbsTol=3e-12);
                        testCase.verifyEqual(b.bottomBuoyancyMatrix,Bd,RelTol=3e-10,AbsTol=3e-12);
                        testCase.verifyEqual(b.boundaryResponseMatrix,response,AbsTol=3e-11);
                        testCase.verifyEqual(b.generalizedEnergyMatrix(g0=-3,gd=.4),H-3*B0+.4*Bd,RelTol=3e-10,AbsTol=3e-12);
                    end
                end
            end
        end

        function theManuscriptBasisHasSignedEnergyAndOrthonormalResponses(testCase)
            b = makeEndpointBasis(IMSolverAnalytical(),"freeSurface",[1.2 .7]);
            o = b.orthogonalizeGeneralizedEnergy(g0=-3,gd=.4);
            testCase.verifyTrue(any(o.boundaryDepth<0,"all"));
            for page = 1:2
                C = o.boundaryResponseMatrix(:,:,page);
                H = b.generalizedEnergyMatrix(g0=-3,gd=.4);
                testCase.verifyEqual(C.'*C,eye(2),AbsTol=3e-14);
                testCase.verifyEqual(C.'*H(:,:,page)*C,diag(o.boundaryDepth(:,page)/(2*b.k(page)^2)),AbsTol=3e-12);
            end
            testCase.verifyFalse(isprop(o,"h"));
            testCase.verifyFalse(isprop(o,"eigenvalues"));
        end

        function priorFullSpanScalingAndOrderingDoNotChangePhysicalModes(testCase)
            b = makeEndpointBasis(IMSolverAnalytical(),"freeSurface",1.2);
            original = b.orthogonalizeGeneralizedEnergy(g0=-3,gd=.4);
            changed = b.withNormalization("maxAbsG").selectColumns([2 1]).orthogonalizeGeneralizedEnergy(g0=-3,gd=.4);
            testCase.verifyEqual(changed.boundaryResponseMatrix,original.boundaryResponseMatrix,AbsTol=0);
            testCase.verifyEqual(changed.boundaryDepth,original.boundaryDepth,AbsTol=0);
            testCase.verifyEmpty(changed.constructionRecord.selection.initialColumnIndices);
            testCase.verifyEqual(changed.constructionRecord.selection.sourceDirections,changed.endpointUnscaledResponseMatrix(),AbsTol=0);
        end

        function profileNormalizationSelectionAndCustomRulesStayAligned(testCase)
            b = makeEndpointBasis(IMSolverAnalytical(),"freeSurface",[1.2 .7]).orthogonalizeGeneralizedEnergy(g0=-3,gd=.4);
            z = [-1;-.6;-.1;0];
            changed = b.withNormalization("rmsF").selectColumns([2 1]);
            M = changed.boundaryResponseMatrix;
            raw = b.endpointUnscaledResponseMatrix();
            for page = 1:2
                T = raw(:,:,page)\M(:,:,page);
                for variable = ["F","G"]
                    for order = 0:(1+(variable=="G"))
                        before = b.evaluate(z,variable=variable,derivativeOrder=order,pages=page);
                        after = changed.evaluate(z,variable=variable,derivativeOrder=order,pages=page);
                        testCase.verifyEqual(after,before*T,AbsTol=3e-12);
                    end
                end
                H = b.energyMatrix; newH = changed.energyMatrix;
                testCase.verifyEqual(newH(:,:,page),T.'*H(:,:,page)*T,AbsTol=3e-12);
            end
            testCase.verifyEqual(changed.boundaryDepth,b.boundaryDepth([2 1],:),AbsTol=0);
            calls = 0;
            custom = b.addNormalization("counted",@countRule).withNormalization("counted");
            F = custom.F(z); first = calls;
            custom.G(z); custom.Gzz(z); custom.F(z);
            testCase.verifyEqual(calls,first);
            testCase.verifyGreaterThan(first,0);
            replaced = custom.addNormalization("counted",@(~,~) 7);
            testCase.verifyEqual(custom.F(z),F,AbsTol=0);
            testCase.verifyNotEqual(replaced.F(z),F);
            function value = countRule(~,j)
                calls = calls+1; value = j+1;
            end
        end

        function selectedCustomRulesValidateOnlyRetainedColumns(testCase)
            z = [-1;-.4;0]; pages = [2 1 2];
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=48)}
                b = makeEndpointBasis(solver{1},"freeSurface",[1.2 .7]);
                custom = b.addNormalization("retainedOnly",@(source,j) (j-1)*source.k);
                testCase.verifyError(@() custom.normalizationFactors("retainedOnly"),string(class(b))+":NullNormalization");
                selected = custom.selectColumns(2).selectColumns(1).withNormalization("retainedOnly");
                testCase.verifyEqual(selected.normalizationFactors(),b.k,AbsTol=0);
                expected = b.G(z,pages=pages);
                expected = expected(:,2,:)./reshape(b.k(pages),1,1,[]);
                testCase.verifyEqual(selected.G(z,pages=pages),expected,AbsTol=2e-14);
                testCase.verifyEqual(selected.columnLabels,b.columnLabels(2));
                testCase.verifyEqual(selected.constructionRecord.selection.initialColumnIndices,2);
                empty = custom.selectColumns([]).withNormalization("retainedOnly");
                testCase.verifySize(empty.G(z,pages=pages),[3 0 3]);
            end
        end

        function scientificEndpointNormalizationsCannotBeReplaced(testCase)
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=48)}
                b = makeEndpointBasis(solver{1},"freeSurface",1.2);
                for name = ["endpointResponse","rmsF","rmsG","maxAbsF","maxAbsG","surfaceF","waveEnergy"]
                    testCase.verifyError(@() b.addNormalization(name,@(~,~) 2),string(class(b))+":ReservedNormalization");
                end
                testCase.verifyError(@() b.withNormalization("generalizedEnergy"),string(class(b))+":UnsupportedNormalization");
                custom = b.addNormalization("myScaling",@(~,~) 2);
                replaced = custom.addNormalization("myScaling",@(~,~) 3);
                testCase.verifyEqual(custom.normalizationFactors("myScaling"),2*ones(2,1),AbsTol=0);
                testCase.verifyEqual(replaced.normalizationFactors("myScaling"),3*ones(2,1),AbsTol=0);
                orthogonal = replaced.orthogonalizeGeneralizedEnergy(g0=-3,gd=.4);
                M = orthogonal.boundaryResponseMatrix;
                testCase.verifyEqual(M.'*M,eye(2),AbsTol=5e-15);
                testCase.verifyEqual(orthogonal.normalization,"endpointResponse");
                Hg = orthogonal.generalizedEnergyMatrix(g0=-3,gd=.4);
                testCase.verifyEqual(Hg,diag(diag(Hg)),AbsTol=3e-14);
            end
        end

        function restrictionsNeverRestoreDiscardedDirections(testCase)
            b = makeEndpointBasis(IMSolverAnalytical(),"freeSurface",1.2);
            for columns = {1,[]}
                selected = b.selectColumns(columns{1});
                testCase.verifyError(@() selected.orthogonalizeGeneralizedEnergy(g0=-3,gd=.4),"IMAnalyticalGeostrophicZeroAPVModesBasis:IncompleteEndpointBasis");
            end
            o = b.orthogonalizeGeneralizedEnergy(g0=-3,gd=.4);
            selected = o.selectColumns(2).withNormalization("maxAbsG");
            testCase.verifySize(selected.boundaryResponseMatrix,[2 1]);
            testCase.verifyEqual(selected.canonicalEndpointLabels,["surface","bottom"]);
            testCase.verifyError(@() selected.orthogonalizeGeneralizedEnergy(g0=-3,gd=.4),"IMAnalyticalGeostrophicZeroAPVModesBasis:IncompleteEndpointBasis");
            directions=o.endpointUnscaledResponseMatrix();
            testCase.verifyEqual(selected.constructionRecord.selection.sourceDirections,directions(:,2),AbsTol=0);
        end

        function physicalDualAndSharedRecoveryCovaryForComplexOffSpanData(testCase)
            energy = struct(g0=-3,gd=.4);
            original = makeEndpointBasis(IMSolverAnalytical(),"freeSurface",1.2).withNormalization("maxAbsF");
            changed = original.orthogonalizeGeneralizedEnergy(g0=energy.g0,gd=energy.gd).withNormalization("rmsF");
            T = original.boundaryResponseMatrix\changed.boundaryResponseMatrix;
            [z,w] = legpts(48,original.zDomain); w = w(:);
            [A,metric] = physicalObservation(original,z,w,energy);
            [newA,newMetric] = physicalObservation(changed,z,w,energy);
            Hplus = A'*(abs(metric).*A);
            U = diag(exp(1i*linspace(.1,1.8,size(A,1))));
            A = U*A; newA = U*newA;
            W = diag(metric);
            H = original.generalizedEnergyMatrix(g0=energy.g0,gd=energy.gd);
            newH = changed.generalizedEnergyMatrix(g0=energy.g0,gd=energy.gd);
            B = H\(A'*W); newB = newH\(newA'*W);
            C = chol(Hplus,"lower")'\eye(original.nColumns);
            D = diag(sqrt(abs(metric)));
            F = diag(linspace(.7,1.6,numel(metric)))*D;
            complement = null(A');
            source = A*[1+.3i;-.2+1i]+(.4-.7i)*complement(:,1)+(.3+.2i)*complement(:,2);
            testCase.verifyEqual(newMetric,metric,AbsTol=0);
            testCase.verifyEqual(newA,A*T,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(A'*W*A,H,RelTol=3e-11,AbsTol=3e-12);
            testCase.verifyEqual(newH,T'*H*T,RelTol=3e-11,AbsTol=3e-12);
            testCase.verifyTrue(any(eig(H)<0));
            testCase.verifyEqual(newB,T\B,RelTol=3e-11,AbsTol=3e-12);
            testCase.verifyGreaterThan(norm(source-A*(A\source)),.5);
            for policy = ["sampledLeastSquares","correctedPhysicalDual","extendedSpan"]
                oldOptions = struct(policy=policy,columnLabels=original.columnLabels,physicalDual=B,coefficientScale=C,observationFactor=D,fitFactor=F);
                newOptions = struct(policy=policy,columnLabels=changed.columnLabels,physicalDual=newB,coefficientScale=T\C,observationFactor=D,fitFactor=F);
                if policy=="extendedSpan"
                    oldOptions.extendedSynthesis=[A,complement(:,1)]; oldOptions.extendedColumnLabels=[original.columnLabels,"latent:1"];
                    oldOptions.extendedActiveColumnMask=true(1,3); oldOptions.extendedCoefficientScale=blkdiag(C,1); oldOptions.retainedColumnIndices=1:2;
                    newOptions.extendedSynthesis=[newA,complement(:,1)]; newOptions.extendedColumnLabels=[changed.columnLabels,"latent:1"];
                    newOptions.extendedActiveColumnMask=true(1,3); newOptions.extendedCoefficientScale=blkdiag(T\C,1); newOptions.retainedColumnIndices=1:2;
                end
                oldPairs=namedargs2cell(oldOptions); newPairs=namedargs2cell(newOptions);
                oldRecovery=IMDiscreteTransform.fromObservationMatrices(A,oldPairs{:});
                newRecovery=IMDiscreteTransform.fromObservationMatrices(newA,newPairs{:});
                testCase.verifyEqual(newRecovery.forwardMatrix,T\oldRecovery.forwardMatrix,RelTol=3e-11,AbsTol=3e-12);
                testCase.verifyEqual(newRecovery.physicalDualMatrix*source,T\(oldRecovery.physicalDualMatrix*source),RelTol=3e-11,AbsTol=3e-12);
                testCase.verifyEqual(newRecovery.transformBack(newRecovery.transformForward(source)),oldRecovery.transformBack(oldRecovery.transformForward(source)),RelTol=3e-11,AbsTol=3e-12);
                testCase.verifyEqual(newRecovery.recoverySpecification.fitFactor,oldRecovery.recoverySpecification.fitFactor,AbsTol=0);
                testCase.verifyEqual(newRecovery.recoverySpecification.observationFactor,oldRecovery.recoverySpecification.observationFactor,AbsTol=0);
            end
        end
    end
end

function basis = makeEndpointBasis(solver,surface,k)
background=IMBackground.constantStratification(N0=1,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,g=2,surfaceCondition=surface);
problem=IMGeostrophicZeroAPVModes.atWavenumber(background,k,endpoints=["surface","bottom"]);
basis=solver.solveGeostrophicZeroAPVModes(problem);
end

function response = endpointResponse(basis)
F=basis.F([basis.zDomain(2);basis.zDomain(1)]); G=basis.G([basis.zDomain(2);basis.zDomain(1)]);
response=G; if basis.surfaceBoundary=="freeSurface", response(1,:)=response(1,:)-F(1,:); end
end

function [H,B0,Bd,response] = independentForms(basis)
[z,w]=legpts(64,basis.zDomain); w=w(:);
F=basis.F(z); G=basis.G(z); f0=basis.f0; g=basis.g; k=basis.k;
H=(k^2*(F'*(w.*F))+(f0^2/g^2)*(G'*(w.*basis.N2(z).*G)))/(2*k^4);
if basis.surfaceBoundary=="freeSurface", Fs=basis.F(basis.zDomain(2)); H=H+(f0^2/(2*g*k^4))*(Fs'*Fs); end
response=endpointResponse(basis);
B0=(f0^2/(2*g^2*k^4))*(response(1,:)'*response(1,:));
Bd=(f0^2/(2*g^2*k^4))*(response(2,:)'*response(2,:));
end

function [A,metric] = physicalObservation(basis,z,w,energy)
F=basis.F(z); G=basis.G(z); Fs=basis.F(basis.zDomain(2)); response=endpointResponse(basis);
f0=basis.f0; g=basis.g; k=basis.k;
A=[F;G;Fs;response];
metric=[w/(2*k^2);w.*basis.N2(z)*f0^2/(2*g^2*k^4);f0^2/(2*g*k^4);energy.g0*f0^2/(2*g^2*k^4);energy.gd*f0^2/(2*g^2*k^4)];
end
