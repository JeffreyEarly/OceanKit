function coefficients = applyAssessmentOperator(self,values)
% Apply the operator explicitly prepared for the existing assessment contract.
arguments (Input)
    self (1,1) IMProjection
    values (:,:) double {mustBeFinite}
end
arguments (Output)
    coefficients (:,:) double
end
if size(values,1) ~= self.sampleCount
    error("IMProjection:InvalidShape","values must have one row per sample.");
end
if self.projectionKind == "prescribedDual"
    coefficients = self.projectPhysical(values);
    return
end
if isempty(self.recovery) || ~self.recovery.success
    error("IMProjection:UnavailableRecovery","Sample recovery is unprepared or unresolved. Inspect recovery diagnostics, or use projectPhysical for the separately available continuous dual.");
end
coefficients = self.recovery.transformForward(values);
end
