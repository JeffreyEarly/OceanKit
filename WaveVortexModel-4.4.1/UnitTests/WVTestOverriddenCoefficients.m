classdef WVTestOverriddenCoefficients < WVCoefficients
    properties
        fluxCalls = 0
    end
    methods
        function self = WVTestOverriddenCoefficients(model)
            self@WVCoefficients(model);
        end
        function values = fluxAtTime(self,t,state)
            self.fluxCalls = self.fluxCalls + 1;
            values = fluxAtTime@WVCoefficients(self,t,state);
            values = cellfun(@(value)2*value,values,UniformOutput=false);
        end
    end
end
