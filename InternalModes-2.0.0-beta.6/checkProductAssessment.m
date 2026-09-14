function decision = checkProductAssessment(assessment,options)
% Apply explicit product/reference tolerances without reducing counts.
%
% Reference failures remain inconclusive. Qualified product errors
% exceeding tolerance reject. Strict requests refer to an examined
% count and retain all scientific columns supplied by the caller.
%
% - Topic: Apply product policies
% - Parameter assessment: product assessment structure returned by plan.assess
% - Parameter options.quadraticTolerance: nonnegative product-error tolerance
% - Parameter options.referenceTolerance: maximum qualified reference discrepancy
% - Parameter options.requestedCount: optional strict examined retained count
arguments (Input)
    assessment (1,1) struct
    options.quadraticTolerance (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
    options.referenceTolerance (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
    options.requestedCount (1,:) double {mustBeInteger,mustBePositive} = assessment.plan.retainedCounts(end)
end
if numel(options.requestedCount)~=1 || ~ismember(options.requestedCount,assessment.plan.retainedCounts)
    error("checkProductAssessment:InvalidRequestedCount","A strict requestedCount must be one explicitly examined count.");
end
m=assessment.measurements; statuses=repmat("accepted",height(m),1);
for q=1:height(m)
    inconclusive=false; trustedFailure=false; examined=false;
    for j=1:numel(assessment.evidence)
        e=assessment.evidence{j}; selected=e.firstCounts<=m.retainedCount(q);
        if ~any(selected), continue; end
        examined=true;
        qualified=e.referenceStatus=="qualified" & e.referenceError(q,selected)<=options.referenceTolerance;
        errors=e.error(q,selected);
        trustedFailure=trustedFailure || any(qualified & errors>options.quadraticTolerance);
        inconclusive=inconclusive || any(~qualified | isnan(errors));
    end
    if inconclusive || ~examined, statuses(q)="inconclusive"; end
    if trustedFailure, statuses(q)="rejected"; end
end
consecutive=cumprod(double(statuses=="accepted"))>0;
selected=find(m.retainedCount<=options.requestedCount);
status="accepted";
if any(statuses(selected)=="inconclusive"), status="inconclusive"; end
if any(statuses(selected)=="rejected"), status="rejected"; end
passing=m.retainedCount(consecutive & m.retainedCount<=options.requestedCount);
largest=0; if ~isempty(passing), largest=passing(end); end
decision=struct(status=status,requestedCount=options.requestedCount,requestedCountAccepted=status=="accepted",largestExaminedAcceptedPrefix=largest,prefixDecisions=table(m.retainedCount,statuses,consecutive,VariableNames=["retainedCount","status","consecutiveAccepted"]),tolerances=options,countDescription="Largest count passing the explicitly sampled products; scientific columns remain unchanged");
end
