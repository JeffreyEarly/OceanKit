function comparison = compareProductAssessmentCoverage(assessment,control,options)
% Detect failures missed by this plan using a separate bounded control.
%
% A control must use the same declared inventory and counts, with
% qualified references. This reports observed missed products;
% no failure can be detected outside both examined inventories.
%
% - Topic: Inspect product evidence
% - Parameter control: independently selected or all-products assessment
% - Parameter assessment: product assessment structure returned by plan.assess
% - Parameter options.quadraticTolerance: explicit failure threshold
% - Parameter options.referenceTolerance: required reference qualification
arguments (Input)
    assessment (1,1) struct
    control (1,1) struct
    options.quadraticTolerance (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
    options.referenceTolerance (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
end
if ~isequal(assessment.plan.retainedCounts,control.plan.retainedCounts) || ~isequal(assessment.plan.products,control.plan.products) || ~isequaln(assessment.plan.factors,control.plan.factors) || ~isequaln(assessment.plan.outputs,control.plan.outputs) || ~isequaln(assessment.assessmentIdentity.grids,control.assessmentIdentity.grids) || assessment.costs.minimumReciprocalCondition~=control.costs.minimumReciprocalCondition
    error("compareProductAssessmentCoverage:IncompatibleControl","Coverage comparison requires the same scientific inventory, grids, and examined counts.");
end
shared=intersect(assessment.plan.products.output(assessment.plan.rows),control.plan.products.output(control.plan.rows));
for index=shared.'
    left=assessment.assessmentIdentity.outputs{index}; right=control.assessmentIdentity.outputs{index};
    if ~isequaln(left.projection,right.projection) || ~isequaln(left.references,right.references)
        error("compareProductAssessmentCoverage:IncompatibleControl","Coverage comparison requires identical prepared projections and reference recipes on shared outputs.");
    end
end
sampledDecision=checkProductAssessment(assessment,quadraticTolerance=options.quadraticTolerance,referenceTolerance=options.referenceTolerance);
controlDecision=checkProductAssessment(control,quadraticTolerance=options.quadraticTolerance,referenceTolerance=options.referenceTolerance);
worst=0; limiting=struct(); missed=0;
for j=1:numel(control.evidence)
    e=control.evidence{j}; own=find(assessment.plan.rows==e.inventoryRow,1);
    ownKeys=strings(1,0); ownFirst=zeros(1,0);
    if ~isempty(own), ownKeys=pairKeys(assessment.evidence{own}); ownFirst=assessment.evidence{own}.firstCounts; end
    [present,location]=ismember(pairKeys(e),ownKeys);
    for k=1:numel(e.isZero)
        checked=e.firstCounts(k)<=assessment.plan.retainedCounts;
        if present(k), checked=checked & assessment.plan.retainedCounts<ownFirst(location(k)); end
        q=find(checked); if isempty(q), continue; end
        [candidate,i]=max(e.error(q,k));
        if candidate>options.quadraticTolerance, missed=missed+1; end
        if candidate>worst
            worst=candidate; limiting=struct(interactionId=e.interactionId,channel=e.channel,inputLabels=e.inputLabels(:,k),frequencySigns=e.frequencySigns(:,k),outputId=e.outputId,retainedCount=assessment.plan.retainedCounts(q(i)));
        end
    end
end
qualified=all(control.measurements.referenceStatus=="qualified" & control.measurements.referenceError<=options.referenceTolerance) && all(assessment.measurements.referenceStatus=="qualified" & assessment.measurements.referenceError<=options.referenceTolerance);
comparison=struct(referenceQualified=qualified,observedMissedFailures=missed,worstMissedError=worst,limitingMissedProduct=limiting,falseAcceptance=qualified && sampledDecision.requestedCountAccepted && ~controlDecision.requestedCountAccepted,sampledDecision=sampledDecision,controlDecision=controlDecision,scope="Only failures exposed by this separate supplied control; unseen products remain unqualified");
if ~qualified, comparison.falseAcceptance=NaN; comparison.observedMissedFailures=NaN; end
end

function keys = pairKeys(evidence)
keys=strings(1,size(evidence.inputLabels,2));
for j=1:numel(keys)
    keys(j)=jsonencode(struct(labels=evidence.inputLabels(:,j),signs=evidence.frequencySigns(:,j)));
end
end
