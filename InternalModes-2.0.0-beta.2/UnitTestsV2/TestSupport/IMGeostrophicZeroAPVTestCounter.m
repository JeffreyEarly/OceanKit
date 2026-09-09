classdef IMGeostrophicZeroAPVTestCounter < handle
    properties
        callCount = 0
        rightHandSideCounts = zeros(1,0)
    end

    methods
        function record(self,nRightHandSides)
            self.callCount = self.callCount+1;
            self.rightHandSideCounts(end+1) = nRightHandSides;
        end
    end
end
