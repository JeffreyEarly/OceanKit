classdef WVTestPortableLinearCoefficientForcing < WVForcing
    properties (SetAccess=immutable)
        rate (1,1) double
    end

    methods
        function self = WVTestPortableLinearCoefficientForcing(wvt,rate,options)
            arguments
                wvt WVTransform
                rate (1,1) double {mustBeFinite}
                options.name (1,1) string = "test portable linear coefficient forcing"
            end
            self@WVForcing(wvt,options.name,WVForcingType("Spectral"));
            self.rate = rate;
        end

        function [Fp,Fm,F0] = addSpectralForcing(self,wvt,Fp,Fm,F0)
            Fp = Fp + self.rate*wvt.Ap;
            Fm = Fm + self.rate*wvt.Am;
            F0 = F0 + self.rate*wvt.A0;
        end

        function contract = portableImplementationContract(self)
            payload = struct("name",string(self.name),"forcingTypes",string(self.forcingType),"priority",self.priority,"rate",self.rate);
            contract = self.supportedPortableImplementationContract("WVTestPortableLinearCoefficientForcing",payload);
        end
    end
end
