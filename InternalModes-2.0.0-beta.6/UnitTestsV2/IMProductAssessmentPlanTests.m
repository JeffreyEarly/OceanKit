classdef IMProductAssessmentPlanTests < matlab.unittest.TestCase
    properties (Access = private)
        originalPath
    end
    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            testCase.originalPath=path;
            addpath(fileparts(fileparts(mfilename("fullpath"))));
        end
    end
    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end
    methods (Test)
        function budgetPreflightCallsNoEvaluatorAndIncludesZeros(testCase)
            [inventory,grids,readCounts]=constantInventory();
            testCase.verifyError(@() IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,retainedCounts=1,productBudget=1),"IMProductAssessmentPlan:ProductBudgetExceeded")
            testCase.verifyEqual(readCounts(),[0 0])
            plan=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,retainedCounts=1,productBudget=3);
            testCase.verifyEqual(readCounts(),[0 0])
            result=plan.assess(grids);
            testCase.verifyEqual(plan.reservedProducts,3)
            testCase.verifyEqual(result.costs.structuralZeros,1)
            testCase.verifyEqual(result.costs.nonzeroProducts,2)
            testCase.verifyEqual(result.costs.evaluatedProducts,3)
        end

        function bothSelectionRulesRejectBudgetBeforePreparation(testCase)
            [definitions,~,readCounts]=constantInventory();
            for selection=["fixedSparse","allProducts"]
                testCase.verifyError(@() IMProductAssessmentPlan(definitions.factors,definitions.products,definitions.outputs,retainedCounts=1,productBudget=2,selection=selection),"IMProductAssessmentPlan:ProductBudgetExceeded")
                testCase.verifyEqual(readCounts(),[0 0])
            end
        end

        function retainedOrdinalsPreserveSignMultiplicityAndFixedCounts(testCase)
            [definitions,~]=constantInventory();
            wave=definitions.factors{1}; wave.id="wave"; wave.family="wave";
            wave.labels=["1+","1-","2+","2-"]; wave.ordinals=[1 1 2 2];
            wave.frequencySigns=[1 -1 1 -1]; wave.countRole="retained";
            boundary=definitions.factors{2}; boundary.labels=["surface","bottom"];
            boundary.ordinals=[1 2]; boundary.frequencySigns=[0 0]; boundary.columnKind="endpoint";
            products=definitions.products(1,:); products.factorB=2;
            plan=IMProductAssessmentPlan({wave,boundary},products,definitions.outputs,retainedCounts=[1 2],productBudget=8);
            testCase.verifyEqual(plan.retainedCounts,[1 2])
            testCase.verifyEqual(plan.reservedProducts,8)
            testCase.verifyEqual(nnz(plan.firstCounts{1}<=1),4)
            testCase.verifyEqual(unique(wave.frequencySigns(plan.pairs{1}(1,:))),[-1 1])
            testCase.verifyEqual(unique(plan.pairs{1}(2,:)),[1 2])
            testCase.verifyEqual(plan.outputs{1}.labels,"mean")
        end

        function selectionIsDeterministicCumulativeAndPreservesFixedCoordinates(testCase)
            [inventory,~]=cosineInventory(16);
            fixed=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,retainedCounts=1:16,productBudget=100000);
            repeated=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,interactionIds="cosine",retainedCounts=1:16,productBudget=100000);
            testCase.verifyEqual(fixed.pairs,repeated.pairs)
            a=inventory.factors{1}; pairs=fixed.pairs{1};
            testCase.verifyTrue(any(a.ordinals(pairs(1,:))==1 & a.ordinals(pairs(2,:))==16))
            testCase.verifyTrue(any(a.ordinals(pairs(1,:))==16 & a.ordinals(pairs(2,:))==16))
            testCase.verifyLessThan(fixed.reservedProducts,16^2)
            for count=1:16
                selected=fixed.firstCounts{1}<=count;
                testCase.verifyFalse(any(a.ordinals(pairs(:,selected))>count,"all"))
            end
            [boundaryInventory,~]=constantInventory();
            factor=boundaryInventory.factors{1}; factor.labels=["surface","bottom"]; factor.ordinals=[1 2]; factor.frequencySigns=[0 0]; factor.countRole="fixed"; factor.columnKind="endpoint";
            table=boundaryInventory.products(1,:); table.factorB=1;
            boundaryInventory=productDefinitions({factor},table,boundaryInventory.outputs);
            plan=IMProductAssessmentPlan(boundaryInventory.factors,boundaryInventory.products,boundaryInventory.outputs,retainedCounts=1,productBudget=4);
            testCase.verifyEqual(plan.reservedProducts,4)
        end

        function executionReusesSelectedFactorsAndSeparateRunsDoNotCache(testCase)
            [inventory,grids,readCounts]=constantInventory();
            plan=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,interactionIds="good",retainedCounts=1,productBudget=2);
            first=plan.assess(grids,chunkSize=1); firstCounts=readCounts();
            testCase.verifyEqual(firstCounts,[2 1])
            testCase.verifyEqual(first.costs.evaluationCalls,2)
            testCase.verifyEqual(first.costs.evaluatedProducts,2)
            plan.assess(grids,chunkSize=2);
            testCase.verifyEqual(readCounts(),2*firstCounts)
        end

        function strictFailurePreservesCountsAndUnqualifiedReferenceIsInconclusive(testCase)
            [inventory,grids]=constantInventory();
            result=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,selection="allProducts",retainedCounts=1,productBudget=3).assess(grids);
            decision=checkProductAssessment(result,quadraticTolerance=0.1,referenceTolerance=1e-8,requestedCount=1);
            testCase.verifyEqual(decision.status,"rejected")
            testCase.verifyEqual(decision.requestedCount,1)
            testCase.verifyEqual(result.plan.outputs{1}.labels,"mean")
            [inventory,grids]=constantInventory(referenceStatus="unverified");
            result=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,retainedCounts=1,productBudget=3).assess(grids);
            decision=checkProductAssessment(result,quadraticTolerance=100,referenceTolerance=1e-8);
            testCase.verifyEqual(decision.status,"inconclusive")
            testCase.verifyFalse(decision.requestedCountAccepted)
        end

        function independentControlExposesAnUntestedFailure(testCase)
            [inventory,grids]=constantInventory();
            sparse=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,interactionIds="good",retainedCounts=1,productBudget=2).assess(grids);
            control=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,selection="allProducts",retainedCounts=1,productBudget=3).assess(grids);
            comparison=compareProductAssessmentCoverage(sparse,control,quadraticTolerance=0.1,referenceTolerance=1e-8);
            testCase.verifyTrue(comparison.falseAcceptance)
            testCase.verifyEqual(comparison.observedMissedFailures,1)
            testCase.verifyGreaterThan(comparison.worstMissedError,1)
            testCase.verifyEqual(comparison.limitingMissedProduct.interactionId,"hidden")
            testCase.verifyFalse(sparse.coverage.exhaustiveGuarantee)
            testCase.verifyFalse(sparse.coverage.superpositionGuarantee)
            testCase.verifyEqual(sparse.coverage.prefixCoverage.availableDeclaredProducts,3)
            testCase.verifyEqual(sparse.coverage.prefixCoverage.examinedProducts,2)
            testCase.verifyEqual(sparse.coverage.prefixCoverage.omittedProducts,1)
        end

        function cosineControlAndChunkSizesAgree(testCase)
            [inventory,grids]=cosineInventory(8);
            plan=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,retainedCounts=1:8,productBudget=1000);
            first=plan.assess(grids,chunkSize=3); second=plan.assess(grids,chunkSize=100);
            testCase.verifyEqual(first.measurements.value,second.measurements.value,AbsTol=2e-14)
            testCase.verifyLessThan(max(first.measurements.value(1:6)),1e-12)
            testCase.verifyGreaterThan(first.measurements.value(end),0.5)
            testCase.verifyEqual(first.evidence{1}.frequencySigns,zeros(size(first.evidence{1}.inputLabels)))
        end

        function collectionFactorsKeepPageDerivativeCoefficientAndLabels(testCase)
            solver=IMSolverSpectral(nEVP=32);
            evp=IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            basis=solver.solveEVP(evp,nModes=3);
            smaller=solver.solveEVP(evp,nModes=2);
            collection=IMBasisCollection({basis,smaller},basisIndex=[2 1]);
            factor=collection.productFactor(id="derivative",page=2,variable="u",derivativeOrder=1,coefficient=@(z) 1i*(1+z));
            z=linspace(-1,0,17).';
            expected=collection.evaluate(z,pages=2,variable="u",derivativeOrder=1,columns=[3 1]).*(1i*(1+z));
            testCase.verifyEqual(factor.evaluate(z,[3 1],"reference"),expected,AbsTol=0)
            testCase.verifyEqual(factor.provenance.page,2)
            testCase.verifyEqual(factor.provenance.derivativeOrder,1)
            testCase.verifyEqual(factor.labels,collection.metadata(1).columnLabels)
            smallFactor=collection.productFactor(id="smaller",page=1,variable="u");
            testCase.verifyNumElements(smallFactor.labels,2)
            testCase.verifyNumElements(factor.labels,3)
        end

        function invalidReferencesAndExplicitConditionGuardsAreActionable(testCase)
            [inventory,grids]=constantInventory(negativeNorm=true);
            testCase.verifyError(@() IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,retainedCounts=1,productBudget=3).assess(grids),"IMProductAssessmentPlan:InvalidReferenceMetric")
            [inventory,grids]=constantInventory();
            result=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,interactionIds="good",retainedCounts=1,productBudget=2).assess(grids,minimumReciprocalCondition=1.1);
            testCase.verifyTrue(all(isinf(result.measurements.value)))
            testCase.verifyEqual(result.costs.minimumReciprocalCondition,1.1)
        end

        function mixedSignedEndpointProductsUsePhysicalPairingsAndPositiveNorms(testCase)
            z=[0;0.5;1]; grids=struct(id={"sample","reference"},z={z,z});
            a=struct(id="page-a",family="boundary",labels="surface",ordinals=1,frequencySigns=0,countRole="fixed",columnKind="endpoint",evaluate=@(z,columns,id) [1;3;2]+double(id=="sample")*[0;-1;1],provenance=struct(page=1,variable="F",derivativeOrder=0));
            b=struct(id="page-b-derivative",family="wave",labels="2-",ordinals=1,frequencySigns=-1,countRole="fixed",evaluate=@(z,columns,id) 1i*ones(size(z)),provenance=struct(page=2,variable="G",derivativeOrder=1,coefficient="imaginary"));
            projection=IMProjection.fromPrescribedDual([-2 1 4],3,3,majorantGramMatrix=7,columnLabels="target",provenance=struct(source="signed endpoint dual"));
            reference=struct(gridId="reference",pairingMatrix=[-2 1 4],normMatrix=diag([2 1 4]),targetGramMatrix=3,majorantGramMatrix=7,role="primary",status="qualified",provenance=struct(source="exact endpoint reference"));
            context=struct(projection=projection,references={{reference}});
            output=struct(id="target",family="apv",labels="target",ordinals=1,countRole="fixed",prepare=@(grids) context,provenance=struct(source="exact endpoint dual"));
            products=table("mixed","derivative",1,2,1,VariableNames=["interactionId","channel","factorA","factorB","output"]);
            result=IMProductAssessmentPlan({a,b},products,{output},retainedCounts=1,productBudget=1).assess(grids);
            testCase.verifyEqual(result.measurements.value,sqrt(7/27),AbsTol=2e-14)
            testCase.verifyEqual(result.evidence{1}.frequencySigns,[0;-1])
            testCase.verifyEqual(result.evidence{1}.inputFamilies,["boundary","wave"])
            testCase.verifyEqual(result.evidence{1}.outputFamily,"apv")
        end

        function largeUnrequestedBandsDoNotAllocateTheirCartesianProducts(testCase)
            [inventory,~]=constantInventory(); factor=inventory.factors{1};
            factor.labels=string(1:100000); factor.ordinals=1:100000; factor.frequencySigns=zeros(1,100000); factor.countRole="retained";
            inventory=productDefinitions({factor},inventory.products(1,:),inventory.outputs);
            plan=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,retainedCounts=1,productBudget=1);
            testCase.verifyEqual(plan.pairs{1},[1;1])
            testCase.verifyEqual(plan.reservedProducts,1)
        end

        function coverageComparisonChecksScientificAndNumericalIdentity(testCase)
            [inventory,grids]=constantInventory();
            sparse=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,interactionIds="good",retainedCounts=1,productBudget=2).assess(grids);
            changed=grids; changed(1).z=2*changed(1).z;
            control=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,selection="allProducts",retainedCounts=1,productBudget=3).assess(changed);
            testCase.verifyError(@() compareProductAssessmentCoverage(sparse,control,quadraticTolerance=0.1,referenceTolerance=1e-8),"compareProductAssessmentCoverage:IncompatibleControl")
            factors=inventory.factors; factors{1}.provenance.source="different scientific basis";
            other=productDefinitions(factors,inventory.products,inventory.outputs);
            control=IMProductAssessmentPlan(other.factors,other.products,other.outputs,selection="allProducts",retainedCounts=1,productBudget=3).assess(grids);
            testCase.verifyError(@() compareProductAssessmentCoverage(sparse,control,quadraticTolerance=0.1,referenceTolerance=1e-8),"compareProductAssessmentCoverage:IncompatibleControl")
        end

        function coverageDetectsFailureAtAnEarlierUntestedPrefix(testCase)
            [inventory,grids]=cosineInventory(16);
            sparse=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,retainedCounts=1:9,productBudget=1000).assess(grids);
            control=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,selection="allProducts",retainedCounts=1:9,productBudget=1000).assess(grids);
            se=sparse.evidence; ce=control.evidence;
            se{1}.error(:)=0; ce{1}.error(:)=0;
            ordinal=ce{1}.inputOrdinals; index=find(ordinal(1,:)==5 & ordinal(2,:)==8);
            ce{1}.error(8,index)=1;
            sparse.evidence=se;
            control.evidence=ce;
            comparison=compareProductAssessmentCoverage(sparse,control,quadraticTolerance=0.1,referenceTolerance=1e-8);
            testCase.verifyTrue(comparison.falseAcceptance)
            testCase.verifyEqual(comparison.observedMissedFailures,1)
            testCase.verifyEqual(comparison.limitingMissedProduct.retainedCount,8)
        end

        function qualifiedFailureOverridesAnUnrelatedInconclusiveFamily(testCase)
            [inventory,grids]=constantInventory();
            result=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,selection="allProducts",retainedCounts=1,productBudget=3).assess(grids);
            evidence=result.evidence; evidence{1}.referenceStatus="inconclusive";
            result.evidence=evidence;
            decision=checkProductAssessment(result,quadraticTolerance=0.1,referenceTolerance=1e-8);
            testCase.verifyEqual(decision.status,"rejected")
        end

        function sampledConditionFailureAndReferenceFailureHaveDifferentStatuses(testCase)
            [inventory,grids]=constantInventory(sampleGram=0);
            result=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,interactionIds="good",retainedCounts=1,productBudget=2).assess(grids,minimumReciprocalCondition=1e-13);
            decision=checkProductAssessment(result,quadraticTolerance=0.1,referenceTolerance=1e-8);
            testCase.verifyEqual(decision.status,"rejected")
            testCase.verifyEqual(result.measurements.referenceError,0)
            [inventory,grids]=constantInventory(referenceGram=0);
            result=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,interactionIds="good",retainedCounts=1,productBudget=2).assess(grids,minimumReciprocalCondition=1e-13);
            decision=checkProductAssessment(result,quadraticTolerance=0.1,referenceTolerance=1e-8);
            testCase.verifyEqual(decision.status,"inconclusive")
        end

        function sparseCountRequestsStillRequireTheCompleteOutputBand(testCase)
            [inventory,~]=cosineInventory(8); outputs=inventory.outputs;
            outputs{1}.labels="mode8"; outputs{1}.ordinals=8;
            inventory=productDefinitions(inventory.factors,inventory.products,outputs);
            testCase.verifyError(@() IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,retainedCounts=8,productBudget=1000),"IMProductAssessmentPlan:InvalidOutputCounts")
        end

        function coverageComparisonRequiresIdenticalNumericalGuards(testCase)
            [inventory,grids]=constantInventory(); plan=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,retainedCounts=1,productBudget=3);
            sampled=plan.assess(grids); control=plan.assess(grids,minimumReciprocalCondition=1.1);
            testCase.verifyError(@() compareProductAssessmentCoverage(sampled,control,quadraticTolerance=0.1,referenceTolerance=1e-8),"compareProductAssessmentCoverage:IncompatibleControl")
        end

        function independentReferenceDifferenceRemainsSeparate(testCase)
            [inventory,grids]=constantInventory(independentScale=1.1);
            result=IMProductAssessmentPlan(inventory.factors,inventory.products,inventory.outputs,interactionIds="good",retainedCounts=1,productBudget=2).assess(grids);
            testCase.verifyLessThan(result.measurements.value,1e-13)
            testCase.verifyGreaterThan(result.measurements.independentSolveError,0.09)
            testCase.verifyEqual(result.measurements.quadratureError,0)
            decision=checkProductAssessment(result,quadraticTolerance=0.1,referenceTolerance=1e-4);
            testCase.verifyEqual(decision.status,"inconclusive")
        end
    end
end

function [inventory,grids,readCounts] = constantInventory(options)
arguments
    options.referenceStatus (1,1) string = "qualified"
    options.negativeNorm (1,1) logical = false
    options.independentScale (1,1) double = 1
    options.sampleGram (1,1) double = 1
    options.referenceGram (1,1) double = 1
end
calls=0; preparations=0;
z=linspace(0,1,5).'; zr=linspace(0,1,9).';
grids=struct(id={"sample","reference"},z={z,zr});
if options.independentScale~=1, grids(3)=struct(id="independent",z=zr); end
base=struct(id="good",family="scalar",labels="mean",ordinals=1,frequencySigns=0,countRole="fixed",evaluate=@good,provenance=struct(source="constant-exact"));
bad=base; bad.id="bad"; bad.evaluate=@badValues;
output=struct(id="output",family="scalar",labels="mean",ordinals=1,countRole="fixed",prepare=@prepare,provenance=struct(source="constant-exact"));
products=table(["good";"good";"hidden"],["nonzero";"zero";"hidden"],[1;1;2],[1;1;1],ones(3,1),[1;0;1],VariableNames=["interactionId","channel","factorA","factorB","output","coefficient"]);
inventory=productDefinitions({base,bad},products,{output}); readCounts=@getCounts;
    function counts=getCounts()
        counts=[calls preparations];
    end
    function v=good(zz,columns,id)
        calls=calls+1; v=ones(numel(zz),numel(columns));
        if id=="independent", v=v*options.independentScale; end
    end
    function v=badValues(zz,columns,id)
        calls=calls+1; v=ones(numel(zz),numel(columns));
        if id=="sample", v=3*v; end
    end
    function context=prepare(grid)
        preparations=preparations+1;
        w=trapezoidWeights(grid(1).z);
        projection=IMProjection.fromPrescribedDual(w.',options.sampleGram,1,columnLabels="mean",provenance=struct(source="constant dual"));
        references=cell(1,numel(grid)-1);
        for j=2:numel(grid)
            w=trapezoidWeights(grid(j).z); metric=diag(w); if options.negativeNorm, metric=-metric; end
            role="primary"; if j>2, role="independent"; end
            references{j-1}=struct(gridId=grid(j).id,pairingMatrix=w.',normMatrix=metric,targetGramMatrix=options.referenceGram,majorantGramMatrix=1,role=role,status=options.referenceStatus,provenance=struct(source="explicit-exact-reference",independentSolve=j>2));
        end
        context=struct(projection=projection,references={references});
    end
end

function [inventory,grids] = cosineInventory(n)
z=linspace(0,pi,10).'; zr=linspace(0,pi,257).';
grids=struct(id={"sample","reference"},z={z,zr});
labels=string(0:n-1);
factor=struct(id="cosine",family="scalar",labels=labels,ordinals=1:n,frequencySigns=zeros(1,n),countRole="retained",evaluate=@(z,columns,id) cos(z*(columns-1)),provenance=struct(source="analytical-cosines"));
output=struct(id="cosine",family="scalar",labels=labels,ordinals=1:n,countRole="retained",prepare=@prepare,provenance=struct(source="analytical-cosines"));
products=table("cosine","cos*cos",1,1,1,VariableNames=["interactionId","channel","factorA","factorB","output"]);
inventory=productDefinitions({factor},products,{output});
    function context=prepare(grid)
        target=diag([pi repmat(pi/2,1,n-1)]);
        w=trapezoidWeights(grid(1).z); projection=IMProjection(cos(grid(1).z*(0:n-1)),diag(w),target,columnLabels=labels);
        w=trapezoidWeights(grid(2).z);
        reference=struct(gridId="reference",pairingMatrix=cos(grid(2).z*(0:n-1)).'*diag(w),normMatrix=diag(w),targetGramMatrix=target,majorantGramMatrix=target,role="primary",status="qualified",provenance=struct(source="resolved-trigonometric-quadrature"));
        context=struct(projection=projection,references={{reference}});
    end
end

function weights=trapezoidWeights(z)
weights=[diff(z(1:2));z(3:end)-z(1:end-2);diff(z(end-1:end))]/2;
end

function definitions = productDefinitions(factors,products,outputs)
definitions = struct(factors={factors},products=products,outputs={outputs});
end
