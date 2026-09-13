classdef WVTestInheritedCoefficients < WVCoefficients
    methods
        function self = WVTestInheritedCoefficients(model)
            self@WVCoefficients(model);
        end
    end
end
