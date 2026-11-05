classdef IMUnqualifiedSecondModeTestProblem < IMEigenvalueProblem
    % Inject a qualification failure without removing its physical position.

    methods
        function self = IMUnqualifiedSecondModeTestProblem()
            self@IMEigenvalueProblem(name="qualificationTest",zDomain=[-1 0],p=1,q=0,r=1);
        end
    end

    methods (Hidden)
        function [A,B,map] = prepareSolvePencil(self,solver,A,B)
            [A,B,map] = prepareSolvePencil@IMEigenvalueProblem(self,solver,A,B);
            map.publicEigenvalues = @(~,lambda) self.qualifyEigenvalues(lambda);
        end
    end

    methods (Static,Access = private)
        function mapped = qualifyEigenvalues(lambda)
            qualified = true(size(lambda));
            [~,indices] = sort(lambda);
            qualified(indices(2)) = false;
            mapped = struct("values",lambda,"qualified",qualified);
        end
    end
end
