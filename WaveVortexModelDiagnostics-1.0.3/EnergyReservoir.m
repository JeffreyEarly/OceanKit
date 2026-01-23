classdef EnergyReservoir
   properties
      hasMDA logical
      hasGeostrophicKinetic logical
      hasGeostrophicPotential logical
      hasInertial logical
      hasInternalGravityWave logical
      name
      fancyName
   end
   properties (Dependent)
        vectorContents
   end
   methods
       function er = EnergyReservoir(hasMDA, hasGeostrophicKinetic, hasGeostrophicPotential, hasInertial, hasInternalGravityWave,name,fancyName)
         er.hasMDA = hasMDA;
         er.hasGeostrophicKinetic = hasGeostrophicKinetic;
         er.hasGeostrophicPotential = hasGeostrophicPotential;
         er.hasInertial = hasInertial;
         er.hasInternalGravityWave = hasInternalGravityWave;
         er.name = name;
         er.fancyName = fancyName;
       end

       function v = get.vectorContents(self)
            v = [self.hasMDA,self.hasGeostrophicKinetic,self.hasGeostrophicPotential,self.hasInertial,self.hasInternalGravityWave].';
       end

       function k = kFromKRadial(self,k)
           switch self
               case {EnergyReservoir.geostrophic_kinetic, EnergyReservoir.geostrophic_potential,EnergyReservoir.geostrophic,EnergyReservoir.igw}
                   k = k(2:end);
               case {EnergyReservoir.mda,EnergyReservoir.io}
                   k = k(1);
           end
       end
   end
   enumeration
      mda                       (true,false,false,false,false,"te_mda","mean density anomaly")
      geostrophic_kinetic       (false,true,false,false,false,"ke_g","geostrophic kinetic")
      geostrophic_kinetic_mda   (false,true,false,false,false,"ke_g","geostrophic kinetic + mda")
      geostrophic_potential     (false,false,true,false,false,"pe_g","geostrophic potential")
      geostrophic_potential_mda (true,false,true,false,false,"pe_g","geostrophic potential + mda")
      geostrophic               (false,true,true,false,false,"te_g","geostrophic")
      geostrophic_mda           (true,true,true,false,false,"te_gmda","geostrophic + mda")
      igw                       (false,false,false,false,true,"te_igw","internal gravity wave")
      io                        (false,false,false,true,false,"te_io","inertial")
      wave                      (false,false,false,true,true,"te_wave","wave")
      total                     (true,true,true,true,true,"te_quadratic","total quadratic")
      damp                      (true,true,true,true,true,"te_damp","closure region")
   end

   methods (Static)
        function eFlux = energyFluxForReservoirFromStructure(Ejk,reservoirNames)
            eFlux = cell(length(reservoirNames),1);
            for iReservoir = 1:length(reservoirNames)
                switch reservoirNames(iReservoir)
                    case EnergyReservoir.geostrophic_kinetic
                        eFlux{iReservoir} = Ejk.KE0(:,2:end,:);
                    case EnergyReservoir.geostrophic_kinetic_mda
                        eFlux{iReservoir} = Ejk.KE0;
                    case EnergyReservoir.geostrophic_potential
                        eFlux{iReservoir} = Ejk.PE0(:,2:end,:);
                    case EnergyReservoir.geostrophic_potential_mda
                        eFlux{iReservoir} = Ejk.PE0;
                    case EnergyReservoir.geostrophic
                        eFlux{iReservoir} = Ejk.KE0(:,2:end,:)+Ejk.PE0(:,2:end,:);
                    case EnergyReservoir.mda
                        eFlux{iReservoir} = Ejk.PE0(:,1,:);
                    case EnergyReservoir.geostrophic_mda
                        eFlux{iReservoir} = Ejk.KE0 + Ejk.PE0;
                    case EnergyReservoir.igw
                        eFlux{iReservoir} = Ejk.Ep(:,2:end,:) + Ejk.Em(:,2:end,:);
                    case EnergyReservoir.io
                        eFlux{iReservoir} = Ejk.Ep(:,1,:) + Ejk.Em(:,1,:);
                    case EnergyReservoir.wave
                        eFlux{iReservoir} = Ejk.Ep+Ejk.Em;
                    case EnergyReservoir.total
                        eFlux{iReservoir} = Ejk.Ep+Ejk.Em+Ejk.KE0+Ejk.PE0;
                    otherwise
                        error("unknown energy reservoir");
                end
            end
        end
   end
end