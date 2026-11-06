function result = assess(self,grids,options)
% Evaluate only reserved products and reuse preparation across their prefixes.
%
% `grids` is a struct row with distinct id and finite coordinate column z;
% its first id must be sample. Output prepare(grids) returns projection and
% a cell row of references. Each reference supplies gridId, pairingMatrix,
% normMatrix, targetGramMatrix, majorantGramMatrix, role (primary,
% quadrature, independent), status, and explicit provenance. Exactly one
% primary reference defines errors; other references measure coefficient and
% normalization changes without replacing the assessed scientific modes.
%
% - Topic: Execute assessments
% - Parameter grids: named sample and explicit reference grids
% - Parameter options.chunkSize: maximum factor/product columns per temporary batch
% - Parameter options.minimumReciprocalCondition: explicit coefficient-system guard, default zero
% - Returns result: measurements, individual evidence, coverage, costs, and separate policies
arguments (Input)
    self (1,1) IMProductAssessmentPlan
    grids (1,:) struct
    options.chunkSize (1,1) double {mustBeInteger,mustBePositive} = 128
    options.minimumReciprocalCondition (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 0
end
if isempty(grids) || ~all(isfield(grids,["id","z"]))
    error("IMProductAssessmentPlan:InvalidGrids","Supply named sample and reference coordinate grids.");
end
ids=string({grids.id});
if ids(1)~="sample" || numel(unique(ids))~=numel(ids) || any(ismissing(ids) | strlength(ids)==0)
    error("IMProductAssessmentPlan:InvalidGrids","Grid IDs must be distinct and begin with sample.");
end
for j=1:numel(grids)
    if ~isnumeric(grids(j).z) || ~isreal(grids(j).z) || size(grids(j).z,2)~=1 || isempty(grids(j).z) || any(~isfinite(grids(j).z))
        error("IMProductAssessmentPlan:InvalidGrids","Each grid needs a finite real nonempty coordinate column.");
    end
end
% The plan's constructor already reserved products, before this first callback.
started=tic; inventory=self; counts=self.retainedCounts;
outputIndices=unique(inventory.products.output(self.rows),"stable").';
outputData=cell(size(inventory.outputs)); preparationSeconds=0;
for index=outputIndices
    entry=inventory.outputs{index}; timer=tic;
    context=entry.prepare(grids);
    context=validateContext(context,entry,grids,ids);
    context.primary=find(cellfun(@(r) string(r.role)=="primary",context.references));
    if numel(context.primary)~=1
        error("IMProductAssessmentPlan:InvalidReferences","Each output requires exactly one primary reference.");
    end
    context.prefix=cell(size(counts)); context.referenceSystems=cell(numel(counts),numel(context.references));
    context.conditionAccepted=true(size(counts)); context.sampleConditionAccepted=true(size(counts)); context.referenceConditionAccepted=true(numel(counts),numel(context.references));
    for q=1:numel(counts)
        if entry.countRole=="fixed", n=numel(entry.labels); else, n=nnz(entry.ordinals<=counts(q)); end
        context.prefix{q}=context.projection.prefix(n);
        active=context.prefix{q}.activeColumnMask;
        sample=context.prefix{q}.gramMatrix(active,active);
        context.conditionAccepted(q)=isempty(sample) || (rank(sample)==nnz(active) && rcond(sample)>=options.minimumReciprocalCondition);
        context.sampleConditionAccepted(q)=context.conditionAccepted(q);
        for r=1:numel(context.references)
            reference=context.references{r}; target=reference.targetGramMatrix(1:n,1:n); target=target(active,active);
            good=isempty(target) || (rank(target)==nnz(active) && rcond(target)>=options.minimumReciprocalCondition);
            context.referenceSystems{q,r}=struct(active=active,accepted=good,solver=[]);
            context.referenceConditionAccepted(q,r)=good;
            if good && ~isempty(target), context.referenceSystems{q,r}.solver=decomposition(target); end
            context.conditionAccepted(q)=context.conditionAccepted(q) && good;
        end
    end
    preparationSeconds=preparationSeconds+toc(timer); outputData{index}=context;
end
% Cache keys are the fixed inventory factor identity, selected columns, and
% exact named grid within this execution. Nothing persists between calls.
selectedColumns=cell(numel(inventory.factors),numel(grids));
for j=1:numel(self.rows)
    row=inventory.products(self.rows(j),:); context=outputData{row.output};
    needed=unique([1 cellfun(@(r) find(ids==string(r.gridId),1),context.references)]);
    for side=1:2
        if side==1, factor=row.factorA; else, factor=row.factorB; end
        for g=needed
            selectedColumns{factor,g}=reshape(union(selectedColumns{factor,g},self.pairs{j}(side,:)),1,[]);
        end
    end
end
values=cell(size(selectedColumns)); evaluationCalls=0; evaluatedFactorColumns=0;
timer=tic;
for f=1:size(values,1)
    for g=1:size(values,2)
        columns=selectedColumns{f,g}; if isempty(columns), continue; end
        values{f,g}=zeros(numel(grids(g).z),numel(columns));
        for first=1:options.chunkSize:numel(columns)
            batch=first:min(first+options.chunkSize-1,numel(columns));
            sampled=inventory.factors{f}.evaluate(grids(g).z,columns(batch),ids(g));
            if ~isnumeric(sampled) || ~isequal(size(sampled),[numel(grids(g).z) numel(batch)]) || any(~isfinite(sampled),"all")
                error("IMProductAssessmentPlan:InvalidFactorValues","Factor %s must return finite nZ-by-selected-column values.",inventory.factors{f}.id);
            end
            values{f,g}(:,batch)=sampled;
            evaluationCalls=evaluationCalls+1; evaluatedFactorColumns=evaluatedFactorColumns+numel(batch);
        end
    end
end
evaluationSeconds=toc(timer); info=whos("values","outputData"); cacheBytes=sum([info.bytes]); peakTemporaryBytes=0;
evidence=cell(size(self.rows)); measured=0; structuralZeros=0; timer=tic;
for j=1:numel(self.rows)
    row=inventory.products(self.rows(j),:); context=outputData{row.output};
    a=inventory.factors{row.factorA}; b=inventory.factors{row.factorB}; pairs=self.pairs{j}; nPairs=size(pairs,2);
    errors=nan(numel(counts),nPairs); referenceErrors=zeros(numel(counts),nPairs);
    quadratureErrors=zeros(numel(counts),nPairs); independentErrors=zeros(numel(counts),nPairs); isZero=false(1,nPairs);
    for first=1:options.chunkSize:nPairs
        batch=first:min(first+options.chunkSize-1,nPairs); sampled=productsOnGrid(1,batch);
        refs=cell(size(context.references));
        for r=1:numel(refs)
            reference=context.references{r}; g=find(ids==string(reference.gridId),1);
            products=productsOnGrid(g,batch);
            pairings=reference.pairingMatrix*products;
            normSquared=real(sum(conj(products).*(reference.normMatrix*products),1));
            zero=all(products==0,1);
            if any(normSquared(~zero)<=0) || any(zero & any(sampled~=0,1))
                error("IMProductAssessmentPlan:InvalidReferenceNorm","A nonzero reference product requires a positive norm; an exact zero reference must also be zero on the sampled grid.");
            end
            refs{r}=struct(pairings=pairings,normSquared=normSquared,isZero=zero,coefficients={cell(size(counts))});
        end
        primary=refs{context.primary}; isZero(batch)=primary.isZero;
        for q=1:numel(counts)
            p=context.prefix{q}; n=p.columnCount; active=p.activeColumnMask;
            for r=1:numel(refs)
                system=context.referenceSystems{q,r}; coefficient=zeros(n,numel(batch));
                pairing=refs{r}.pairings(1:n,:);
                if any(pairing(~active,:)~=0,"all")
                    error("IMProductAssessmentPlan:InactiveReferencePairing","Inactive coordinates require zero signed reference pairings.");
                end
                if system.accepted && any(active), coefficient(active,:)=system.solver\pairing(active,:); end
                refs{r}.coefficients{q}=coefficient;
            end
            if ~context.referenceConditionAccepted(q,context.primary)
                errors(q,batch)=Inf; referenceErrors(q,batch)=Inf;
                continue
            end
            if context.sampleConditionAccepted(q)
                errors(q,batch)=p.productError(sampled,refs{context.primary}.coefficients{q},primary.normSquared);
            else
                errors(q,batch)=Inf;
            end
            for r=1:numel(refs)
                if r==context.primary, continue; end
                if ~context.referenceConditionAccepted(q,r)
                    referenceErrors(q,batch)=Inf;
                    if string(context.references{r}.role)=="quadrature", quadratureErrors(q,batch)=Inf; end
                    if string(context.references{r}.role)=="independent", independentErrors(q,batch)=Inf; end
                    continue
                end
                delta=refs{r}.coefficients{q}-refs{context.primary}.coefficients{q};
                discrepancy=IMProjection.relativeCoefficientNorm(delta,p.majorantGramMatrix,primary.normSquared);
                normalization=abs(sqrt(refs{r}.normSquared./primary.normSquared)-1);
                discrepancy=max(discrepancy,normalization);
                bothZero=primary.isZero & refs{r}.isZero;
                discrepancy(bothZero)=0;
                discrepancy(primary.isZero & ~refs{r}.isZero)=Inf;
                if any(~isfinite(discrepancy) & ~isinf(discrepancy)), discrepancy(isnan(discrepancy))=Inf; end
                referenceErrors(q,batch)=max(referenceErrors(q,batch),discrepancy);
                if string(context.references{r}.role)=="quadrature", quadratureErrors(q,batch)=max(quadratureErrors(q,batch),discrepancy); end
                if string(context.references{r}.role)=="independent", independentErrors(q,batch)=max(independentErrors(q,batch),discrepancy); end
            end
        end
        info=whos("sampled","refs","products","pairings"); peakTemporaryBytes=max(peakTemporaryBytes,sum([info.bytes]));
        measured=measured+numel(batch); structuralZeros=structuralZeros+nnz(primary.isZero);
    end
    statuses=cellfun(@(r) string(r.status),context.references);
    referenceStatus="qualified"; if any(statuses~="qualified"), referenceStatus="inconclusive"; end
    evidence{j}=struct(inventoryRow=self.rows(j),interactionId=row.interactionId,channel=row.channel,inputFamilies=[a.family b.family],outputFamily=inventory.outputs{row.output}.family,outputId=inventory.outputs{row.output}.id,inputFactorIds=[a.id b.id],inputLabels=[a.labels(pairs(1,:));b.labels(pairs(2,:))],frequencySigns=[a.frequencySigns(pairs(1,:));b.frequencySigns(pairs(2,:))],inputOrdinals=[a.ordinals(pairs(1,:));b.ordinals(pairs(2,:))],firstCounts=self.firstCounts{j},error=errors,referenceError=referenceErrors,quadratureError=quadratureErrors,independentSolveError=independentErrors,isZero=isZero,referenceStatus=referenceStatus,referenceProvenance={cellfun(@(r) r.provenance,context.references,UniformOutput=false)},conditionAccepted=context.conditionAccepted,sampleConditionAccepted=context.sampleConditionAccepted,referenceConditionAccepted=context.referenceConditionAccepted,outputColumnLabels=inventory.outputs{row.output}.labels);
end
assert(measured==self.reservedProducts,"Executed products differ from the reserved plan.");
assessmentSeconds=toc(timer);
costs=struct(planningSeconds=self.planningSeconds,preparationSeconds=preparationSeconds,evaluationSeconds=evaluationSeconds,assessmentSeconds=assessmentSeconds,totalSeconds=toc(started),productBudget=self.productBudget,reservedProducts=self.reservedProducts,evaluatedProducts=measured,nonzeroProducts=measured-structuralZeros,structuralZeros=structuralZeros,evaluationCalls=evaluationCalls,evaluatedFactorColumns=evaluatedFactorColumns,cachedArrayBytes=cacheBytes,peakTemporaryArrayBytes=peakTemporaryBytes,memoryDescription="MATLAB whos array payloads; excludes allocator/process overhead",minimumReciprocalCondition=options.minimumReciprocalCondition);
outputIdentities=cell(size(outputData));
for index=outputIndices
    outputIdentities{index}=struct(projection=outputData{index}.projection,references={outputData{index}.references});
end
identity=struct(grids=grids,outputs={outputIdentities});
result=summarizeAssessment(self,evidence,costs,identity);

    function products = productsOnGrid(grid,batch)
        [~,ia]=ismember(pairs(1,batch),selectedColumns{row.factorA,grid});
        [~,ib]=ismember(pairs(2,batch),selectedColumns{row.factorB,grid});
        products=row.coefficient.*values{row.factorA,grid}(:,ia).*values{row.factorB,grid}(:,ib);
    end
end

function context = validateContext(context,entry,grids,ids)
if ~isstruct(context) || ~isscalar(context) || ~all(isfield(context,["projection","references"])) || ~isa(context.projection,"IMProjection") || ~iscell(context.references) || isempty(context.references)
    error("IMProductAssessmentPlan:InvalidOutput","Output preparation must return projection and a nonempty reference cell row.");
end
p=context.projection; n=numel(entry.labels);
if p.columnCount~=n || p.sampleCount~=numel(grids(1).z) || ~isequal(p.columnLabels,entry.labels)
    error("IMProductAssessmentPlan:InvalidOutput","Output projection must preserve declared sample count and exact scientific column labels.");
end
for j=1:numel(context.references)
    r=context.references{j}; fields=["gridId","pairingMatrix","normMatrix","targetGramMatrix","majorantGramMatrix","role","status","provenance"];
    if ~isstruct(r) || ~isscalar(r) || ~all(isfield(r,fields))
        error("IMProductAssessmentPlan:InvalidReferences","References require explicit grid, signed pairings, positive metric, coefficient systems, role, status, and provenance.");
    end
    g=find(ids==string(r.gridId),1);
    if isempty(g) || g==1 || ~ismember(string(r.role),["primary","quadrature","independent"]) || ~ismember(string(r.status),["qualified","unverified","inconclusive"]) || ~isstruct(r.provenance) || isempty(fieldnames(r.provenance))
        error("IMProductAssessmentPlan:InvalidReferences","References must identify a separate declared grid, valid role/status, and explicit provenance.");
    end
    nz=numel(grids(g).z);
    if ~isnumeric(r.pairingMatrix) || ~isequal(size(r.pairingMatrix),[n nz]) || any(~isfinite(r.pairingMatrix),"all") || ~isequal(size(r.normMatrix),[nz nz]) || ~isreal(r.normMatrix) || any(~isfinite(r.normMatrix),"all") || norm(r.normMatrix-r.normMatrix.',Inf)>100*eps(max(1,norm(r.normMatrix,Inf)))
        error("IMProductAssessmentPlan:InvalidReferences","Reference pairings and positive metric must match the output columns and named grid.");
    end
    % Reference norms may contain unweighted observation coordinates; positive
    % semidefiniteness is required, then each nonzero product is checked.
    if isdiag(r.normMatrix)
        positive=all(diag(r.normMatrix)>=0);
    else
        positive=min(eig(full((r.normMatrix+r.normMatrix.')/2)))>=-100*eps(max(1,norm(r.normMatrix,2)));
    end
    if ~positive, error("IMProductAssessmentPlan:InvalidReferenceMetric","Reference error metrics must be positive semidefinite, never an absolute signed quadratic form."); end
    % Reuse the projection boundary validation for signed reference systems.
    IMProjection.fromPrescribedDual(zeros(n,0),r.targetGramMatrix,r.targetGramMatrix,majorantGramMatrix=r.majorantGramMatrix,activeColumnMask=p.activeColumnMask,columnLabels=p.columnLabels,provenance=r.provenance);
end
end

function result = summarizeAssessment(plan,evidence,costs,identity)
result = struct();
counts=plan.retainedCounts; n=numel(counts);
if numel(evidence)~=numel(plan.rows)
    error("IMProductAssessmentPlan:InvalidEvidence","Evidence must include every reserved product-family row.");
end
value=zeros(n,1); referenceError=zeros(n,1); quadratureError=zeros(n,1); independentSolveError=zeros(n,1);
status=repmat("measured",n,1); referenceStatus=repmat("qualified",n,1);
examinedCount=zeros(n,1); zeroCount=zeros(n,1); limitingInputI=strings(n,1); limitingInputJ=strings(n,1);
limitingInteraction=strings(n,1); limitingChannel=strings(n,1); limitingOutput=strings(n,1);
for q=1:n
    hasValue=false;
    for j=1:numel(evidence)
        e=evidence{j}; selected=e.firstCounts<=counts(q); indices=find(selected);
        if isempty(indices), continue; end
        examinedCount(q)=examinedCount(q)+nnz(selected); zeroCount(q)=zeroCount(q)+nnz(e.isZero(selected));
        [worst,i]=max(e.error(q,selected));
        if ~hasValue || worst>value(q)
            value(q)=worst; index=indices(i); limitingInputI(q)=e.inputLabels(1,index); limitingInputJ(q)=e.inputLabels(2,index);
            limitingInteraction(q)=e.interactionId; limitingChannel(q)=e.channel; limitingOutput(q)=e.outputId; hasValue=true;
        end
        referenceError(q)=max(referenceError(q),max(e.referenceError(q,selected)));
        quadratureError(q)=max(quadratureError(q),max(e.quadratureError(q,selected)));
        independentSolveError(q)=max(independentSolveError(q),max(e.independentSolveError(q,selected)));
        if e.referenceStatus~="qualified", referenceStatus(q)="inconclusive"; status(q)="inconclusive"; end
    end
    if ~hasValue, value(q)=NaN; status(q)="inconclusive"; end
end
retainedCount=counts(:); quantity=repmat("quadraticAliasing",n,1);
result.measurements=table(retainedCount,quantity,value,status,referenceStatus,referenceError,quadratureError,independentSolveError,examinedCount,zeroCount,limitingInputI,limitingInputJ,limitingInteraction,limitingChannel,limitingOutput);
products=plan.products; families=strings(numel(evidence),3); productCounts=zeros(numel(evidence),1);
for j=1:numel(evidence)
    e=evidence{j}; families(j,:)=[e.inputFamilies e.outputFamily]; productCounts(j)=numel(e.isZero);
end
[uniqueFamilies,~,group]=unique(families,"rows","stable");
examinedFamilies=table(uniqueFamilies(:,1),uniqueFamilies(:,2),uniqueFamilies(:,3),accumarray(group,productCounts),VariableNames=["inputA","inputB","output","productCount"]);
availableDeclaredProducts=zeros(n,1); selectedFamilyAvailableProducts=zeros(n,1);
for q=1:n
    for row=1:height(products)
        a=plan.factors{products.factorA(row)}; b=plan.factors{products.factorB(row)};
        na=numel(a.labels); nb=numel(b.labels);
        if a.countRole=="retained", na=nnz(a.ordinals<=counts(q)); end
        if b.countRole=="retained", nb=nnz(b.ordinals<=counts(q)); end
        availableDeclaredProducts(q)=availableDeclaredProducts(q)+na*nb;
        if ismember(row,plan.rows), selectedFamilyAvailableProducts(q)=selectedFamilyAvailableProducts(q)+na*nb; end
    end
end
prefixCoverage=table(counts(:),availableDeclaredProducts,selectedFamilyAvailableProducts,examinedCount,selectedFamilyAvailableProducts-examinedCount,availableDeclaredProducts-examinedCount,VariableNames=["retainedCount","availableDeclaredProducts","selectedFamilyAvailableProducts","examinedProducts","omittedWithinSelectedFamilies","omittedProducts"]);
result.coverage=struct(selection=plan.selection,selectedRows=plan.rows,omittedRows=setdiff(1:height(products),plan.rows),selectedInteractions=unique(products.interactionId(plan.rows),"stable"),omittedInteractions=setdiff(unique(products.interactionId,"stable"),unique(products.interactionId(plan.rows),"stable"),"stable"),examinedFamilies=examinedFamilies,prefixCoverage=prefixCoverage,selectedFamilyRowCount=numel(plan.rows),totalFamilyRowCount=height(products),exhaustiveGuarantee=false,superpositionGuarantee=false,scope="Only supplied valid interactions and explicit physical source recipes; no model qualification");
result.plan=plan; result.evidence=evidence; result.costs=costs; result.assessmentIdentity=identity;
end
