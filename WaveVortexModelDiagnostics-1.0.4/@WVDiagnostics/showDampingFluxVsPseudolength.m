function showDampingFluxVsPseudolength(self,options)
% Show Damping Flux Vs Pseudolength.
%
% A plot focused on the damping flux as a function of pseudolength scale, for both exact and quadratic fluxes.
%
% - Topic: Figures — Energy
% - Declaration: showDampingFluxVsPseudolength(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
arguments
    self WVDiagnostics
    options.timeIndices = Inf
end
    options.energyReservoirs = [EnergyReservoir.geostrophic_mda, EnergyReservoir.wave];

    kp = self.kPseudoRadial;
    forcing_fluxes_kp = self.quadraticEnergyFluxesTemporalAverage(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices);
    forcing_exact_kp = self.exactEnergyFluxesTemporalAverage(timeIndices=options.timeIndices);

    % First convert the forcing fluxes to kPseudoRadial
    for iForce=1:length(forcing_fluxes_kp)
        forcing_fluxes_kp(iForce).te_gmda = self.transformToPseudoRadialWavenumberA0(forcing_fluxes_kp(iForce).te_gmda);
        forcing_fluxes_kp(iForce).te_wave = self.transformToPseudoRadialWavenumberApm(forcing_fluxes_kp(iForce).te_wave);
    end
    for iForce=1:length(forcing_exact_kp)
        forcing_exact_kp(iForce).te = self.transformToPseudoRadialWavenumberA0(forcing_exact_kp(iForce).te);
    end

    damping_index = [forcing_fluxes_kp.name] == "adaptive_damping";
    damping_kp = forcing_fluxes_kp(damping_index).te_gmda + forcing_fluxes_kp(damping_index).te_wave;
    
filter = @(v) abs(cumsum(v))/self.flux_scale;

figure
tiledlayout(2,1)
nexttile
plot(2*pi./kp,filter(damping_kp)), yscale('log')
xlabel('meters')
ylabel('abs(cumsum(damping))')
nexttile
plot(2*pi./kp,filter(forcing_exact_kp(damping_index).te)), yscale('log')
xlabel('meters')
ylabel('abs(cumsum(damping))')

end