% Assess a bounded cosine-product inventory on a caller-owned sample grid.
root=fileparts(fileparts(mfilename("fullpath"))); addpath(root);
n=8; sampleZ=linspace(0,pi,10).'; referenceZ=linspace(0,pi,129).';
grids=struct(id={"sample","reference"},z={sampleZ,referenceZ});
labels=string(0:n-1);
factor=struct(id="cosine",family="scalar",labels=labels,ordinals=1:n,frequencySigns=zeros(1,n),countRole="retained",evaluate=@(z,columns,id) cos(z*(columns-1)),provenance=struct(source="analytical cosine basis"));
output=struct(id="cosine",family="scalar",labels=labels,ordinals=1:n,countRole="retained",prepare=@(grids) prepareCosines(grids,n,labels),provenance=struct(source="positive L2 projection"));
products=table("supplied-interaction","cosine product",1,1,1,VariableNames=["interactionId","channel","factorA","factorB","output"]);
plan=IMProductAssessmentPlan({factor},products,{output},retainedCounts=1:n,productBudget=1000);
assessment=plan.assess(grids,chunkSize=16);
decision=checkProductAssessment(assessment,quadraticTolerance=0.1,referenceTolerance=1e-12,requestedCount=n);
disp(assessment.measurements)
disp(decision)

function context=prepareCosines(grids,n,labels)
target=diag([pi repmat(pi/2,1,n-1)]);
z=grids(1).z; w=trapWeights(z);
projection=IMProjection(cos(z*(0:n-1)),diag(w),target,columnLabels=labels);
z=grids(2).z; w=trapWeights(z);
reference=struct(gridId="reference",pairingMatrix=cos(z*(0:n-1)).'*diag(w),normMatrix=diag(w),targetGramMatrix=target,majorantGramMatrix=target,role="primary",status="qualified",provenance=struct(source="trigonometric products resolved by reference rule"));
context=struct(projection=projection,references={{reference}});
end

function w=trapWeights(z)
w=[diff(z(1:2));z(3:end)-z(1:end-2);diff(z(end-1:end))]/2;
end
