classdef WVTestUnregisteredFixedAmplitudeForcing < WVFixedAmplitudeForcing
    methods
        function self = WVTestUnregisteredFixedAmplitudeForcing(wvt)
            self@WVFixedAmplitudeForcing(wvt,name="test unregistered fixed amplitude");
        end
    end
end
