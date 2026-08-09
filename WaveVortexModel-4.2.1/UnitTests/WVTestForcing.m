classdef WVTestForcing < WVForcing
    properties (GetAccess=public, SetAccess=protected)
        value
    end

    methods
        function self = WVTestForcing(wvt,name,forcingType,priority,value)
            self@WVForcing(wvt,name,forcingType);
            self.priority = priority;
            self.value = value;
        end

        function [Fu,Fv,Feta] = addHydrostaticSpatialForcing(self,~,Fu,Fv,Feta)
            Fu = Fu+self.value;
            Fv = Fv+2*self.value;
            Feta = Feta+3*self.value;
        end

        function [Fu,Fv,Fw,Feta] = addNonhydrostaticSpatialForcing(self,~,Fu,Fv,Fw,Feta)
            Fu = Fu+self.value;
            Fv = Fv+2*self.value;
            Fw = Fw+3*self.value;
            Feta = Feta+4*self.value;
        end

        function [Fp,Fm,F0] = addSpectralForcing(self,~,Fp,Fm,F0)
            Fp = Fp+self.value;
            Fm = Fm+2*self.value;
            F0 = F0+3*self.value;
        end

        function [Ap,Am,A0] = setSpectralAmplitude(self,~,Ap,Am,A0)
            Ap(:) = self.value;
            Am(:) = 2*self.value;
            A0(:) = 3*self.value;
        end

        function F0 = addPotentialVorticitySpatialForcing(self,~,F0)
            F0 = F0+self.value;
        end

        function F0 = addPotentialVorticitySpectralForcing(self,~,F0)
            F0 = F0+self.value;
        end

        function A0 = setPotentialVorticitySpectralAmplitude(self,~,A0)
            A0(:) = self.value;
        end
    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {};
        end

        function annotations = classDefinedPropertyAnnotations()
            annotations = CAPropertyAnnotation.empty(0,0);
        end
    end
end
