classdef WVFastTransformDoublyPeriodicMatlab < WVFastTransformDoublyPeriodic
    properties (Dependent)
        complexBuffer
    end

    properties
        wvg
        Nz
    end

    properties (Access=private)
        complexBufferRows
    end

    methods
        function self = WVFastTransformDoublyPeriodicMatlab(wvg,Nz)
            self.wvg = wvg;
            self.Nz=Nz;
            self.fourierStorageLayout = WVFourierStorageLayout(wvg,"full-complex");
            self.complexBufferRows = self.fourierStorageLayout.allocateFourierStorage(Nz);
        end

        function value = get.complexBuffer(self)
            value = self.fourierStorageLayout.reshapeFourierRowsToStorage(self.complexBufferRows);
        end

        function set.complexBuffer(self,value)
            self.complexBufferRows = self.fourierStorageLayout.reshapeFourierStorageToRows(value);
        end
    end

    methods
        u_bar = transformFromSpatialDomainWithFourier(self,u)
        u = transformToSpatialDomainWithFourier(self,u_bar)
        du = diffX(wvg,u,options)
        du = diffY(wvg,u,options)
    end
end
