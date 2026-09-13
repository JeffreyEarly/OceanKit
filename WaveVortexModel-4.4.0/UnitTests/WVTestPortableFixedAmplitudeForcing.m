classdef WVTestPortableFixedAmplitudeForcing < WVFixedAmplitudeForcing
    methods
        function self = WVTestPortableFixedAmplitudeForcing(wvt,options)
            arguments
                wvt WVTransform
                options.name (1,1) string = "test portable fixed amplitude"
            end
            self@WVFixedAmplitudeForcing(wvt,name=options.name);
        end

        function contract = portableImplementationContract(self)
            payload = struct("name",string(self.name),"forcingTypes",string(self.forcingType),"priority",self.priority,"ApIndices",self.Ap_indices,"ApValues",self.Apbar,"AmIndices",self.Am_indices,"AmValues",self.Ambar,"A0Indices",self.A0_indices,"A0Values",self.A0bar);
            contract = self.supportedPortableImplementationContract("WVTestPortableFixedAmplitudeForcing",payload);
        end
    end
end
