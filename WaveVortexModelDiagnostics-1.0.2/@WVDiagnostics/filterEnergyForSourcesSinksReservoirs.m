function [sources, sinks, inertial_tx, inertial_cascade, ddt, energy] = filterEnergyForSourcesSinksReservoirs(self,options)
% This function returns values assuming three reservoirs: geo, wave, and.
%
% This function can probably be deprecated after we modernize to the  fluxesForReservoirGroup
%
% This function returns values assuming three reservoirs: geo, wave, and
% damping. The damping resevoir is just scales below a threshold, wave or
% geostrophic. It also returns the exact and exact-damp resevoirs.
% The forcing struct has the the forcing on the three/two different
% reservoirs
% The inertial struct has the flux from the two reservoirs (wave
% geostrophic) to each other and to the damping region.
% The forcing struct also include the nonlinear advection, which has the
% flux to the damping region.
% The ddt struct contains the change in total energy, closing the energy
% budget.
%
% - Topic: Internal — Support functions for createReservoirGroup
% - Declaration: [sources, sinks, inertial_tx, inertial_cascade, ddt, energy] = filterEnergyForSourcesSinksReservoirs(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter customNames: (optional) input argument `customNames` (default: configureDictionary("string","string"))
% - Parameter fluxTolerance: (optional) input argument `fluxTolerance` (default: 1e-2)
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
% - Parameter shouldShowReservoirEnergy: (optional) energy or enstrophy reservoir selection (default: true)
% - Parameter shouldShowExactValues: (optional) input argument `shouldShowExactValues` (default: true)
% - Parameter shouldSeparateClosureRegion: (optional) input argument `shouldSeparateClosureRegion` (default: true)
% - Returns sources: output value `sources`
% - Returns sinks: output value `sinks`
% - Returns inertial_tx: output value `inertial_tx`
% - Returns inertial_cascade: output value `inertial_cascade`
% - Returns ddt: output value `ddt`
% - Returns energy: diagnosed energy as a function of time and/or scale
arguments
    self WVDiagnostics
    options.customNames = configureDictionary("string","string")
    options.fluxTolerance = 1e-2;
    options.timeIndices = Inf;
    options.shouldShowReservoirEnergy = true
    options.shouldShowExactValues = true
    options.shouldSeparateClosureRegion = true
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

    % Now find which indices are outside the damping region
    damping_index = [forcing_fluxes_kp.name] == "adaptive_damping";
    if ~any(damping_index)
        error("Unable to find adaptive damping and thus do not know how to proceed. I guess we could look for another closure.");
    end

    % Method 1
    damping_kp = forcing_fluxes_kp(damping_index).te_gmda + forcing_fluxes_kp(damping_index).te_wave;
    NoDamp = abs(cumsum(damping_kp/self.flux_scale)) < options.fluxTolerance;
    
    % Method 2
    svv = self.wvt.forcingWithName("adaptive damping");
    NoDamp = kp < (svv.k_damp+svv.k_no_damp)/2;

    % Sort the forcing (both quadratic and exact) into damped and undamped
    forcing = self.filterFluxesForReservoir(forcing_fluxes_kp,filter=@(v) sum(sum(v(NoDamp))));
    forcing_fluxes_damp = self.filterFluxesForReservoir(forcing_fluxes_kp,filter=@(v) sum(sum(v(~NoDamp))));
    forcing_exact = self.filterFluxesForReservoir(forcing_exact_kp,filter=@(v) sum(sum(v(NoDamp))));
    forcing_exact_damp = self.filterFluxesForReservoir(forcing_exact_kp,filter=@(v) sum(sum(v(~NoDamp))));

    % create a new reservoir, te_damp, which contains the combined geo and
    % wave forcing in the damped region. Same for the exact flux.
    for iForce=1:length(forcing)
        forcing(iForce).te_damp = forcing_fluxes_damp(iForce).te_gmda + forcing_fluxes_damp(iForce).te_wave;
        forcing(iForce).te_exact = forcing_exact(iForce).te;
        forcing(iForce).te_exact_damp = forcing_exact_damp(iForce).te;
    end

    % inertial_triads = self.filterFluxesForReservoir(inertial_jk,filter=@(v) sum(sum(v(NoDamp))));
    
    % Now deal with the inertial triads. Here we use the *sparse* pseudo
    % radial axis, so we need to find the appropriate indices again.
    [inertial_fluxes_g_kps, inertial_fluxes_w_kps, kps] = self.quadraticEnergyPrimaryTriadFluxesTemporalAverage1D(timeIndices=options.timeIndices);
    NoDampKps = kps <= max(kp(NoDamp));

    % the total transfers are opposite and equal
    % but the transfers into/from the undamped region will not be
    gmda_tx_wave_no_damp = sum(inertial_fluxes_g_kps([inertial_fluxes_g_kps.name] == "tx-wwg").flux(NoDampKps) + inertial_fluxes_g_kps([inertial_fluxes_g_kps.name] == "tx-ggw").flux(NoDampKps));
    wave_tx_gmda_no_damp = sum(inertial_fluxes_w_kps([inertial_fluxes_w_kps.name] == "tx-wwg").flux(NoDampKps) + inertial_fluxes_w_kps([inertial_fluxes_w_kps.name] == "tx-ggw").flux(NoDampKps));
    gmda_tx_wave_damp = sum(inertial_fluxes_g_kps([inertial_fluxes_g_kps.name] == "tx-wwg").flux(~NoDampKps) + inertial_fluxes_g_kps([inertial_fluxes_g_kps.name] == "tx-ggw").flux(~NoDampKps));
    wave_tx_gmda_damp = sum(inertial_fluxes_w_kps([inertial_fluxes_w_kps.name] == "tx-wwg").flux(~NoDampKps) + inertial_fluxes_w_kps([inertial_fluxes_w_kps.name] == "tx-ggw").flux(~NoDampKps));

    % There is one route for te_gmda to transfer to te_wave, and that is
    % the direct transfer terms summed. Although, *some* of that energy
    % might land in the damping region of the waves, and thus actually land
    % in te_damp. The two cascade terms will also transfer to te_damp.
    inertial_tx(1).name = "te_gmda";
    inertial_tx(1).fancyName = WVDiagnostics.fancyNameForName(inertial_tx(1).name);
    inertial_tx(1).te_gmda = 0;
    inertial_tx(1).te_wave = wave_tx_gmda_no_damp;
    inertial_tx(1).te_damp = wave_tx_gmda_damp;

    inertial_tx(2).name = "te_wave";
    inertial_tx(2).fancyName = WVDiagnostics.fancyNameForName(inertial_tx(2).name);
    inertial_tx(2).te_gmda = gmda_tx_wave_no_damp;
    inertial_tx(2).te_wave = 0;
    inertial_tx(2).te_damp = gmda_tx_wave_damp;

    inertial_tx(3).name = "te_damp";
    inertial_tx(3).fancyName = WVDiagnostics.fancyNameForName(inertial_tx(3).name);
    inertial_tx(3).te_gmda = 0;
    inertial_tx(3).te_wave = 0;
    inertial_tx(3).te_damp = 0;
    
    inertial_tx(4).name = "te_exact";
    inertial_tx(4).fancyName = "exact undamped reservoir";
    inertial_tx(4).te_gmda = 0;
    inertial_tx(4).te_damp = forcing(1).te_exact;
    inertial_tx(4).te_wave = 0;

    % Categorization method 1: Only count transfer that both start and land
    % in the un-damped region. Then allocate the rest as transfers to the
    % damped region.
    % if abs(gmda_tx_wave_no_damp) < abs(wave_tx_gmda_no_damp)
    %     % Two possibilities
    %     % 1. gmda (< 0) loses energy to wave (> 0), wave *gain* from damped geo,
    %     % 2. gmda (> 0) gain energy from wave (< 0), but wave more negative than gmda positive, so wave lose to damped geo
    %     gmda_tx_wave = gmda_tx_wave_no_damp;
    %     wave_tx_gmda = -gmda_tx_wave_no_damp;
    %     gmda_tx_damp = 0;
    %     wave_tx_damp = wave_tx_gmda_no_damp + gmda_tx_wave_no_damp;
    % else
    %     gmda_tx_wave = -wave_tx_gmda_no_damp;
    %     wave_tx_gmda = wave_tx_gmda_no_damp;
    %     gmda_tx_damp = wave_tx_gmda_no_damp + gmda_tx_wave_no_damp;
    %     wave_tx_damp = 0;
    % end

    % Categorization method 2: Count transfers as in the undamped region
    % where they originate. The portion that lands in the damped region,
    % add that the forward flux
    % if abs(gmda_tx_wave_no_damp) > abs(wave_tx_gmda_no_damp)
    %     % Two possibilities
    %     % 1. gmda (< 0) loses energy to wave (> 0), wave *gain* from damped geo,
    %     % 2. gmda (> 0) gain energy from wave (< 0), but wave more negative than gmda positive, so wave lose to damped geo
    %     gmda_tx_wave = gmda_tx_wave_no_damp;
    %     wave_tx_gmda = -gmda_tx_wave_no_damp;
    %     gmda_tx_damp = 0;
    %     wave_tx_damp = wave_tx_gmda_no_damp + gmda_tx_wave_no_damp;
    % else
    %     gmda_tx_wave = -wave_tx_gmda_no_damp;
    %     wave_tx_gmda = wave_tx_gmda_no_damp;
    %     gmda_tx_damp = wave_tx_gmda_no_damp + gmda_tx_wave_no_damp;
    %     wave_tx_damp = 0;
    % end

    % but the transfers into/from the undamped region will not be
    g_cascade = sum(inertial_fluxes_g_kps([inertial_fluxes_g_kps.name] == "ggg").flux(NoDampKps) + inertial_fluxes_g_kps([inertial_fluxes_g_kps.name] == "ggw").flux(NoDampKps));
    w_cascade = sum(inertial_fluxes_w_kps([inertial_fluxes_w_kps.name] == "www").flux(NoDampKps) + inertial_fluxes_w_kps([inertial_fluxes_w_kps.name] == "wwg").flux(NoDampKps));
    g_cascade_damp = sum(inertial_fluxes_g_kps([inertial_fluxes_g_kps.name] == "ggg").flux(~NoDampKps) + inertial_fluxes_g_kps([inertial_fluxes_g_kps.name] == "ggw").flux(~NoDampKps));
    w_cascade_damp = sum(inertial_fluxes_w_kps([inertial_fluxes_w_kps.name] == "www").flux(~NoDampKps) + inertial_fluxes_w_kps([inertial_fluxes_w_kps.name] == "wwg").flux(~NoDampKps));

    inertial_cascade(1).name = "te_gmda";
    inertial_cascade(1).fancyName = WVDiagnostics.fancyNameForName(inertial_tx(1).name);
    inertial_cascade(1).te_gmda = 0;
    inertial_cascade(1).te_wave = 0;
    inertial_cascade(1).te_damp = -g_cascade;

    inertial_cascade(2).name = "te_wave";
    inertial_cascade(2).fancyName = WVDiagnostics.fancyNameForName(inertial_tx(2).name);
    inertial_cascade(2).te_gmda = 0;
    inertial_cascade(2).te_wave = 0;
    inertial_cascade(2).te_damp = -w_cascade;

    inertial_cascade(3).name = "te_damp";
    inertial_cascade(3).fancyName = WVDiagnostics.fancyNameForName(inertial_tx(3).name);
    inertial_cascade(3).te_gmda = g_cascade;
    inertial_cascade(3).te_wave = w_cascade;
    inertial_cascade(3).te_damp = 0;

    inertial_cascade(4).name = "te_exact";
    inertial_cascade(4).fancyName = "exact undamped reservoir";
    inertial_cascade(4).te_gmda = 0;
    inertial_cascade(4).te_wave = 0;
    inertial_cascade(4).te_damp = -forcing(1).te_exact;

    % remove nonlinear advection, now that we copied the values we needed
    forcing(1) = [];

    % divide into sources and sinks
    iSink = 1; iSource = 1;
    for iForce=1:length(forcing)
        if forcing(iForce).te_exact + forcing(iForce).te_exact_damp < 0
            sinks(iSink) = forcing(iForce); iSink = iSink + 1;
        else
            sources(iSource) = forcing(iForce); iSource = iSource + 1;
        end
    end

    [reservoirEnergy, t] = self.quadraticEnergyOverTime(energyReservoirs=options.energyReservoirs,timeIndices=options.timeIndices,shouldIncludeExactTotalEnergy=true);

    % Let's be precise here, and pull out the energy for the damped region
    % also
    if self.diagnosticsHasExplicitAntialiasing
        wvt = self.wvt_aa;
    else
        wvt = self.wvt;
    end
    self.iTime = options.timeIndices(end);
    E_pm_kp = self.transformToPseudoRadialWavenumberApm(wvt.transformToRadialWavenumber(wvt.Apm_TE_factor.*(abs(wvt.Ap).^2 + abs(wvt.Am).^2)));
    E_0_kp = self.transformToPseudoRadialWavenumberA0(wvt.transformToRadialWavenumber(wvt.A0_TE_factor.*(abs(wvt.A0).^2)));
    E_damp_pm_f = sum(E_pm_kp(:) .* ~NoDamp(:));
    E_damp_0_f = sum(E_0_kp(:) .* ~NoDamp(:));
    self.iTime = options.timeIndices(1);
    E_pm_kp = self.transformToPseudoRadialWavenumberApm(wvt.transformToRadialWavenumber(wvt.Apm_TE_factor.*(abs(wvt.Ap).^2 + abs(wvt.Am).^2)));
    E_0_kp = self.transformToPseudoRadialWavenumberA0(wvt.transformToRadialWavenumber(wvt.A0_TE_factor.*(abs(wvt.A0).^2)));
    E_damp_pm_i = sum(E_pm_kp(:) .* ~NoDamp(:));
    E_damp_0_i = sum(E_0_kp(:) .* ~NoDamp(:));

    energy.te_gmda = mean(reservoirEnergy(1).energy) - (E_damp_0_i + E_damp_0_f)/2;
    energy.te_wave = mean(reservoirEnergy(2).energy) - (E_damp_pm_i + E_damp_pm_f)/2;
    energy.te_damp = (E_damp_0_i + E_damp_0_f + E_damp_pm_i + E_damp_pm_f)/2;
    energy.te_exact = mean(reservoirEnergy(3).energy);
    
    ddt.te_gmda = (reservoirEnergy(1).energy(end) - reservoirEnergy(1).energy(1))/(t(end)-t(1));
    ddt.te_wave = (reservoirEnergy(2).energy(end) - reservoirEnergy(2).energy(1))/(t(end)-t(1));
    ddt.te_damp = (E_damp_pm_f + E_damp_0_f - E_damp_pm_i - E_damp_0_i)/(t(end)-t(1));
    ddt.te_exact = (reservoirEnergy(3).energy(end) - reservoirEnergy(3).energy(1))/(t(end)-t(1));
end