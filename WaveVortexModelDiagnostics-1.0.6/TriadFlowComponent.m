classdef TriadFlowComponent
   properties
      hasMDA logical
      hasGeostrophic logical
      hasInertial logical
      hasInternalGravityWave logical
      name
      fancyName
   end
   properties (Dependent)
       vectorContents
   end
   methods
       function er = TriadFlowComponent(hasMDA, hasGeostrophic, hasInertial, hasInternalGravityWave,name,fancyName)
         er.hasMDA = hasMDA;
         er.hasGeostrophic = hasGeostrophic;
         er.hasInertial = hasInertial;
         er.hasInternalGravityWave = hasInternalGravityWave;
         er.name = name;
         er.fancyName = fancyName;
       end
       function v = get.vectorContents(self)
           v = [self.hasMDA,self.hasGeostrophic,self.hasInertial,self.hasInternalGravityWave].';
       end
       function fc = flowComponent(self,wvt)
           switch self
               case TriadFlowComponent.geostrophic
                   fc = wvt.geostrophicComponent;
               case TriadFlowComponent.mda
                   fc = wvt.mdaComponent;
               case TriadFlowComponent.geostrophic_mda
                   fc = wvt.geostrophicComponent + wvt.mdaComponent;
               case TriadFlowComponent.igw
                   fc = wvt.waveComponent;
               case TriadFlowComponent.io
                   fc = wvt.inertialComponent;
               case TriadFlowComponent.wave
                   fc = wvt.waveComponent + wvt.inertialComponent;
           end
       end
   end
   enumeration
      mda                       (true,false,false,false,"mda","mean density anomaly")
      geostrophic               (false,true,false,false,"g","geostrophic")
      geostrophic_mda           (true,true,false,false,"gmda","geostrophic + mda")
      igw                       (false,false,false,true,"igw","internal gravity wave")
      io                        (false,false,true,false,"io","inertial")
      wave                      (false,false,true,true,"wave","wave")
   end
end