classdef (Hidden, Sealed) IMDeclaredIntervals
    % Validate a caller's declaration of exact neutral or critical intervals.
    methods (Static)
        function validate(background,intervals,value,kind,errorPrefix)
            arguments
                background (1,1) IMBackground
                intervals (:,2) double {mustBeReal,mustBeFinite}
                value (1,1) double {mustBeReal,mustBeFinite}
                kind (1,1) string
                errorPrefix (1,1) string
            end
            if isempty(intervals) || any(intervals(:,1) >= intervals(:,2)) || any(intervals(:,1) < background.zDomain(1)) || any(intervals(:,2) > background.zDomain(2)) || any(intervals(2:end,1) <= intervals(1:end-1,2))
                error(errorPrefix+":InvalidIntervals","Declare ordered, separated nonempty %s intervals inside zDomain.",kind);
            end
            % Sampling rejects contradictions. Exact constancy remains a
            % scientific declaration, never an inference from sampled data.
            for j = 1:size(intervals,1)
                z = linspace(intervals(j,1),intervals(j,2),33).';
                values = background.N2(z);
                if ~isequal(size(values),size(z)) || ~isreal(values) || any(~isfinite(values)) || any(values ~= value)
                    error(errorPrefix+":InconsistentInterval","The declared %s interval contradicts the background N2 profile. Exact N2=%g is required throughout each interval.",kind,value);
                end
            end
        end
    end
end
