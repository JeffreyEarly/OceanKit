classdef EtaTrueOperationToolboxUnavailable < EtaTrueOperation

    methods
        function self = EtaTrueOperationToolboxUnavailable(wvt)
            self@EtaTrueOperation(wvt);
        end
    end

    methods (Access=protected)
        function tf = hasOptimizationToolboxSupport(self)
            tf = false;
        end
    end
end
