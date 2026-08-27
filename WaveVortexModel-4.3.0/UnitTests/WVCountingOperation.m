classdef WVCountingOperation < WVOperation
    properties (SetAccess=private)
        callCount (1,1) double = 0
        valuesFunction function_handle
    end

    methods
        function self = WVCountingOperation(name,outputVariables,valuesFunction)
            self@WVOperation(name,outputVariables,@disp);
            self.valuesFunction = valuesFunction;
        end

        function varargout = compute(self,wvt,varargin)
            self.callCount = self.callCount+1;
            values = self.valuesFunction(wvt,varargin{:});
            if length(values) ~= self.nVarOut
                error("The test operation returned an unexpected number of values.")
            end
            varargout = values;
        end
    end
end
