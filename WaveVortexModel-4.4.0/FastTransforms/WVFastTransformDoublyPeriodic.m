classdef WVFastTransformDoublyPeriodic < handle
    properties (SetAccess=protected)
        fourierStorageLayout
    end

    methods (Abstract)
        u_bar = transformFromSpatialDomainWithFourier(self,u)
        u = transformToSpatialDomainWithFourier(self,u_bar)
        du = diffX(wvg,u,options)
        du = diffY(wvg,u,options)
    end
end
