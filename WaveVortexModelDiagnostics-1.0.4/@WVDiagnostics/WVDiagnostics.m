classdef WVDiagnostics < handle & CAAnnotatedClass
    %WVDiagnostics Produces diagnostics and figures from WVModel output
    %   This is a collection of diagnostic tools for analyzing model output
    %
    % - Topic: Initialization
    % - Topic: Diagnostics Generation
    % - Topic: Figures — Model Snapshot
    % - Topic: Figures — Energy
    % - Topic: Figures — Potential Enstrophy
    % - Topic: Figures — Ancillary
    % - Topic: Summaries
    % - Topic: Diagnostics — Energy
    % - Topic: Diagnostics — Energy — Time series, [t 1]
    % - Topic: Diagnostics — Energy Fluxes
    % - Topic: Diagnostics — Energy Fluxes — General, [j kRadial t]
    % - Topic: Diagnostics — Energy Fluxes — Temporal averages, [j kRadial]
    % - Topic: Diagnostics — Energy Fluxes — Time series, [t 1]
    % - Topic: Diagnostics — Energy Fluxes — Spatial-temporal averages, [1 1]
    % - Topic: Diagnostics — Energy Fluxes — Temporal averages, 1D axes
    % - Topic: Diagnostics — Energy Fluxes — Temporal averages, 2D axes [sparseJWavenumberAxis sparseKRadialAxis]
    % - Topic: Diagnostics — Potential Enstrophy
    % - Topic: Diagnostics — Potential Enstrophy — Time series, [t 1]
    % - Topic: Diagnostics — Potential Enstrophy Fluxes
    % - Topic: Diagnostics — Potential Enstrophy Fluxes — General, [j kRadial t]
    % - Topic: Diagnostics — Potential Enstrophy Fluxes — Temporal averages, [j kRadial]
    % - Topic: Diagnostics — Potential Enstrophy Fluxes — Time series, [t 1]
    % - Topic: Diagnostics — Potential Enstrophy Fluxes — Spatial-temporal averages, [1 1]
    % - Topic: Computing Spectra
    % - Topic: Transformations — Cosine and Sine
    % - Topic: Transformations — Axes
    % - Topic: Utilities
    % - Topic: Utilities — Sparse matrices
    % - Topic: Internal — Support functions for createReservoirGroup
    properties
        wvpath      % path to the WaveVortexModel output
        diagpath    % path to the diagnostics file

        % wvaapath - path to a WVTransform with explicit anti-aliasing
        % The WVTransformBoussinesq can take a very long time to
        % initialize, so if explicity anti-aliasing is requested, we cache
        % a copy of the transform.
        wvaapath
        wvfile      % NetCDFFile instance for the WaveVortexModel output
        diagfile    % NetCDFFile instance for the diagnostics file
        wvt         % WVTransform instance, set to the current iTime

        % time axis scaling (numeric value)
        % The default is 86400 seconds, corresponding to 1 day.
        %
        % Set both tscale and tscale_units simultaneously to change the
        % axis in the figures.
        % - Topic: Units
        % - Declaration: tscale
        % - Returns
        tscale = 86400

        % time axis units (string value)
        % The default is `days'.
        %
        % Set both tscale and tscale_units simultaneously to change the
        % axis in the figures.
        % - Topic: Units
        tscale_units = "days"

        % energy scaling (numeric value)
        % This is a depth-integrated area-averaged energy.
        %
        % The default is 3.74 $m^3 s^{-2}$, corresponding to 1 GM.
        %
        % Set both escale and escale_units simultaneously to change the
        % axis in the figures.
        % - Topic: Units
        escale = 3.74

        % energy units (string value)
        % The default is `GM'.
        %
        % Set both escale and escale_units simultaneously to change the
        % axis in the figures.
        % - Topic: Units
        escale_units = "GM"

        % energy flux scaling (numeric value)
        % This is a depth-integrated area-averaged energy per unit time.
        %
        % The default is 3.74/(86400*365), which corresponds to 1 GM per
        % year.
        %
        % Set both flux_scale and flux_scale_units simultaneously to change the
        % axis in the figures.
        % - Topic: Units
        flux_scale = 3.74/(86400*365)

        % energy flux units (string value)
        % The default is `GM/yr'.
        %
        % Set both flux_scale and flux_scale_units simultaneously to change the
        % axis in the figures.
        % - Topic: Units
        flux_scale_units = "GM/yr"

        % potential enstrophy scaling (numeric value)
        % This is a depth-integrated area-averaged potential enstrophy.
        %
        % The default units are $m f^2$, where f gets set upon
        % initialization of the WaveVortexModel.
        %
        % Set both zscale and zscale_units simultaneously to change the
        % axis in the figures.
        % - Topic: Units
        zscale

        % potential enstrophy units (string value)
        %
        % The default units are $m f^2$, where f gets set upon
        % initialization of the WaveVortexModel.
        %
        % Set both zscale and zscale_units simultaneously to change the
        % axis in the figures.
        % - Topic: Units
        zscale_units = "m f^2";

        % potential enstrophy flux scaling (numeric value)
        % This is a depth-integrated area-averaged potential enstrophy per unit time..
        %
        % The default units are $m f^2/yr$, where f gets set upon
        % initialization of the WaveVortexModel.
        %
        % Set both z_flux_scale and z_flux_scale_units simultaneously to change the
        % axis in the figures.
        % - Topic: Units
        z_flux_scale

        % potential enstrophy flux units (string value)
        % This is a depth-integrated area-averaged potential enstrophy per unit time..
        %
        % The default units are $m f^2/yr$, where f gets set upon
        % initialization of the WaveVortexModel.
        %
        % Set both z_flux_scale and z_flux_scale_units simultaneously to change the
        % axis in the figures.
        % - Topic: Units
        z_flux_scale_units = "m f^2/yr"

        % choose an algorithm for pseudoRadialBinning
        pseudoRadialBinning (1,1) string {mustBeMember(pseudoRadialBinning,["k2+j2","adaptive"])} = "k2+j2"
    end

    properties (Access=private)
        wvt_aa_cache
    end

    properties (Dependent)
        % Get time vector from the diagnostics file
        %
        % Reads the 't' variable from the diagnostics file. This might not be the same as the time vector in the model output file.
        %
        % - Topic: Dependent property getter
        % - Declaration: t = self.t_diag
        % - Returns t: time vector from diagnostics file
        t_diag
        t_wv
        j
        jWavenumber
        kRadial
        kPseudoRadial
        omegaAxis
        kePeAxis
        Lr2
        Lr2_pm
        omega_jk
        geo_hke_jk
        geo_pe_jk
        forcingNames
        wvt_aa
        diagnosticsHasExplicitAntialiasing
    end

    properties (SetObservable)
        iTime = Inf
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Diagnostics Generation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        createDiagnosticsFile(self,options)
        create1DMirrorFluxes(self,options)
        create2DMirrorFluxes(self,options)
        createReservoirGroup(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Internal — Support functions for createReservoirGroup
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [triadVar, forcingVar, energyVar] = variablesForReservoirGroup(self,options)
        addTriadFluxesForReservoirGroupAtTime(self,options)
        [transferFlux, forcingFlux, ddt, energy] = fluxesForReservoirGroup(self,options)
        [sources, sinks, inertial_tx, inertial_cascade, ddt, energy] = filterEnergyForSourcesSinksReservoirs(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Figures — Model Snapshot
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = plotFluidStateMultipanel(self,options)
        fig = plotFluidDecompositionMultipanel(self,options)
        fig = plotEnergySpectrum(self,options)
        fig = plotEnstrophySpectrum(self,options)
        fig = plotPotentialEnergySpectrum(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Figures — Energy
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = plotEnergyOverTime(self,options)
        fig = plotEnergyFluxOverTime(self,options)
        fig = plotEnergyTriadFluxOverTime(self,options)
        fig = plotEnergyFluxTemporalAverage(self,options)

        fig = plotEnergyFluxes1D(self,options)
        [fig, boxDiagram] = plotSourcesSinksForReservoirGroup(self,options)
        fig = plotSourcesSinksReservoirsDiagram(self,options)
        showDampingFluxVsPseudolength(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Figures — Potential Enstrophy
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = plotEnstrophyOverTime(self,options)
        fig = plotEnstrophyFluxOverTime(self,options)
        fig = plotEnstrophyTriadFluxOverTime(self,options)
        fig = plotEnstrophyFluxTemporalAverage(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Figures — Ancillary
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = plotMooringRotarySpectrum(self)
        setLogWavelengthXAxis(self,options)
        [labels, ticks] = logWavelengthAxis(self,options)
        overlayFrequencyContours(self,options)
        overlayGeostrophicKineticPotentialRatioContours(self,options)
        overlayGeostrophicKineticPotentialFractionContours(self,options)
        showRossbyRadiusYAxis(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Summaries
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tableString = createEnstrophyFluxSummaryTable(self,options)
        summarizeSourcesSinksReservoirs(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Diagnostics — Energy
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % - Topic: Diagnostics — Energy — Time series, [t 1]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [energy, t] = exactEnergyOverTime(self, options)
        [reservoirs, t] = quadraticEnergyOverTime(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Diagnostics — Energy Fluxes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % - Topic: Diagnostics — Energy Fluxes — General, [j kRadial t]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        energy_fluxes = exactEnergyFluxes(self)
        forcing_fluxes = quadraticEnergyFluxes(self,options)
        inertial_fluxes = quadraticEnergyTriadFluxes(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % - Topic: Diagnostics — Energy Fluxes — Temporal averages, [j kRadial]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        energy_fluxes = quadraticEnergyFluxesTemporalAverage(self,options)
        inertial_fluxes = quadraticEnergyTriadFluxesTemporalAverage(self,options)
        energy_fluxes = exactEnergyFluxesTemporalAverage(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % - Topic: Diagnostics — Energy Fluxes — Time series, [t 1]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [energy_fluxes, t] = exactEnergyFluxesOverTime(self,options)
        [energy_fluxes,t] = quadraticEnergyFluxesOverTime(self,options)
        [energy_fluxes,t] = quadraticEnergyTriadFluxesOverTime(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % - Topic: Diagnostics — Energy Fluxes — Spatial-temporal averages, [1 1]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        energy_fluxes = exactEnergyFluxesSpatialTemporalAverage(self,options)
        forcing_fluxes = quadraticEnergyFluxesSpatialTemporalAverage(self,options)
        inertial_fluxes = quadraticEnergyTriadFluxesSpatialTemporalAverage(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % - Topic: Diagnostics — Energy Fluxes — Temporal averages, 1D axes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [M_wwg, M_ggw, kp] = quadraticEnergyMirrorTriadFluxes1D(self,options)
        [M_wwg, omegaAxis] = quadraticEnergyMirrorTriadFluxes1D_omega(self,options)
        [M_ggw, kePeAxis] = quadraticEnergyMirrorTriadFluxes1D_kepe(self,options)

        [inertial_fluxes_g, inertial_fluxes_w, kp] = quadraticEnergyPrimaryTriadFluxesTemporalAverage1D(self,options)
        [inertial_fluxes_w, omegaAxis] = quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_omega(self,options)
        [inertial_fluxes_g, kePeAxis] = quadraticEnergyPrimaryTriadFluxesTemporalAverage1D_kepe(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % - Topic: Diagnostics — Energy Fluxes — Temporal averages, 2D axes [sparseJWavenumberAxis sparseKRadialAxis]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [M_wwg, F_wwg, ks, js] = quadraticEnergyMirrorTriadFluxes2D(self,options)
        [inertial_fluxes_g, inertial_fluxes_w, ks, js] = quadraticEnergyPrimaryTriadFluxesTemporalAverage2D(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Diagnostics — Potential Enstrophy
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % - Topic: Diagnostics — Potential Enstrophy — Time series, [t 1]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [enstrophy, t] = exactEnstrophyOverTime(self, options)
        [enstrophy, t] = quadraticEnstrophyOverTime(self, options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Diagnostics — Potential Enstrophy Fluxes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % - Topic: Diagnostics — Potential Enstrophy Fluxes — General, [j kRadial t]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        enstrophy_fluxes = exactEnstrophyFluxes(self)
        enstrophy_fluxes = quadraticEnstrophyFluxes(self)
        inertial_fluxes = quadraticEnstrophyTriadFluxes(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % - Topic: Diagnostics — Potential Enstrophy Fluxes — Temporal averages, [j kRadial]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        enstrophy_fluxes = quadraticEnstrophyFluxesTemporalAverage(self,options)
        enstrophy_fluxes = exactEnstrophyFluxesTemporalAverage(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % - Topic: Diagnostics — Potential Enstrophy Fluxes — Time series, [t 1]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [enstrophy_fluxes, t] = exactEnstrophyFluxesOverTime(self,options)
        [enstrophy_fluxes,t] = quadraticEnstrophyFluxesOverTime(self,options)
        [enstrophy_fluxes,t] = quadraticEnstrophyTriadFluxesOverTime(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % - Topic: Diagnostics — Potential Enstrophy Fluxes — Spatial-temporal averages, [1 1]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        enstrophy_fluxes = exactEnstrophyFluxesSpatialTemporalAverage(self,options)
        enstrophy_fluxes = quadraticEnstrophyFluxesSpatialTemporalAverage(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Transformations — Axes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [varargout] = transformToPseudoRadialWavenumber(self,energyReservoir,varargin);
        [varargout] = transformToPseudoRadialWavenumberA0(self,varargin);
        [varargout] = transformToPseudoRadialWavenumberApm(self,varargin)
        [varargout] = transformToOmegaAxis(self,varargin)
        [varargout] = transformToKePeAxis(self,varargin)

        [kp,bins_0,bins_pm] = sparsePseudoRadialAxis(self)
        [omegaAxis,bins_omega] = sparseOmegaAxis(self)
        [kePeAxis,bins_kepe] = sparseKePeAxis(self)
        k = sparseKRadialAxis(self)
        j = sparseJWavenumberAxis(self)
        [S_0, S_pm, mask_0, mask_pm] = sparseJKAxisBinMatrices(self)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Computing Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        S_f = spectrumWithFgTransform(self,f)
        S_f = spectrumWithGgTransform(self,f)
        S_f = crossSpectrumWithFgTransform(self,phi,gamma)
        S_f = crossSpectrumWithGgTransform(self,phi,gamma)


        fig = plotPoissonFlowOverPcolor(wvd,options)
        fig = plotPoissonFlowOverContours(wvd,options)
        [logX,logY,Uprime,Vprime] = RescalePoissonFlowFluxForLogSpace(wvd,X,Y,U,V,options)
        [X,Y,U,V] = PoissonFlowFromFlux(wvd, flux)
        [X,Y,U,V] = PoissonFlowFromFluxType1(wvd, flux)


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Dependent property implementations
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function wvt_aa = get.wvt_aa(self)
            if isempty(self.wvt_aa_cache) || ~isvalid(self.wvt_aa_cache)
                if isa(self.wvt,"WVTransformBoussinesq")
                    if ~exist(self.wvaapath,"file")
                        fprintf("No existing Boussinesq transform with explicity antialiasing found. Creating a new transform.\n");
                        wvt_ = self.wvt.waveVortexTransformWithExplicitAntialiasing();
                        ncfile_ = wvt_.writeToFile(self.wvaapath);
                        ncfile_.close;
                    end
                    self.wvt_aa_cache = WVTransform.waveVortexTransformFromFile(self.wvaapath,iTime=self.iTime,shouldReadOnly=true);
                else
                    self.wvt_aa_cache = self.wvt.waveVortexTransformWithExplicitAntialiasing();
                end
                self.wvt_aa_cache.addOperation(SpatialForcingOperation(self.wvt_aa_cache));
            end
            wvt_aa = self.wvt_aa_cache;
        end

        function t = get.t_diag(self)
            t = self.diagfile.readVariables('t');
        end

        function t = get.j(self)
            % j axis
            %
            % Reads the 'j' variable from the diagnostics file, if available, otherwise the model file.
            %
            % - Topic: Axes
            % - Declaration: j
            % - Returns j: array of values
            if ~isempty(self.diagfile)
                t = self.diagfile.readVariables('j');
            else
                t = self.wvt.j;
            end
        end

        function t = get.kRadial(self)
            % Get kRadial indices
            %
            % Reads the 'kRadial' variable from the diagnostics file, if available, otherwise the model file.
            %
            % - Topic: Axes
            % - Declaration: kRadial
            % - Returns kRadial: array of values
            if ~isempty(self.diagfile)
                t = self.diagfile.readVariables('kRadial');
            else
                t = self.wvt.kRadial;
            end
        end

        function t = get.Lr2(self)
            % Hydrostatic deformation radius squared
            %
            % Reads the 'Lr2' variable from the diagnostics file if it
            % exists, or from the wvt if not
            %
            % - Topic: Axes
            % - Declaration: Lr2
            % - Returns Lr2: array of values
            if ~isempty(self.diagfile)
                t = self.diagfile.readVariables('Lr2');
            else
                t = self.wvt.Lr2;
            end
        end

        function t = get.Lr2_pm(self)
            % Wave deformation radius squared
            %
            % Reads the 'Lr2_pm' variable from the diagnostics file if it
            % exists, or from the wvt if not
            %
             % - Topic: Axes
            % - Declaration: Lr2_pm
            % - Returns Lr2_pm: array of values
            if ~isempty(self.diagfile)
                t = self.diagfile.readVariables('Lr2_pm');
            else
                t = self.wvt.g*self.wvt.h_pm/self.wvt.f/self.wvt.f;
            end
        end

        function t = get.omega_jk(self)
            % Wave frequency omega in jk space
            %
            % Reads the 'omega_jk' variable from the diagnostics file if it
            % exists, or from the wvt if not
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.omega_jk(self)
            % - Returns t: omega_jk matrix from diagnostics file
            if ~isempty(self.diagfile)
                if self.diagfile.hasVariableWithName('omega_jk')
                    t = self.diagfile.readVariables('omega_jk');
                else
                    t = self.diagfile.readVariables('omega_axis');
                end
            else
                [omegaN,n] = self.wvt.transformToRadialWavenumber(abs(self.wvt.Omega),ones(size(self.wvt.Omega)));
                t = (omegaN./n);
            end
        end

        function t = get.geo_hke_jk(self)
            % Geostrophic kinetic energy in jk space
            %
            % Reads the 'geo_hke_jk' variable from the diagnostics file if it
            % exists, or from the wvt if not
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.geo_hke_jk(self)
            % - Returns t: geo_hke_jk matrix from diagnostics file
            if ~isempty(self.diagfile)
                if self.diagfile.hasVariableWithName('geo_hke_jk')
                    t = self.diagfile.readVariables('geo_hke_jk');
                else
                    t = self.diagfile.readVariables('geo_hke_axis');
                end
            else
                t = self.wvt.transformToRadialWavenumber(self.wvt.A0_KE_factor);
            end
        end

        function t = get.geo_pe_jk(self)
            % Geostrophic potential energy in jk space
            %
            % Reads the 'geo_pe_jk' variable from the diagnostics file if it
            % exists, or from the wvt if not
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.geo_pe_jk(self)
            % - Returns t: geo_pe_jk matrix from diagnostics file
            if ~isempty(self.diagfile)
                if self.diagfile.hasVariableWithName('geo_pe_jk')
                    t = self.diagfile.readVariables('geo_pe_jk');
                else
                    t = self.diagfile.readVariables('geo_pe_axis');
                end
            else
                t = self.wvt.transformToRadialWavenumber(self.wvt.A0_PE_factor);
            end
        end

        function jWavenumber = get.jWavenumber(self)
            jWavenumber = 1./sqrt(self.Lr2);
            jWavenumber(1) = 0; % barotropic mode is a mean?
        end

        function kPseudoRadial = get.kPseudoRadial(self)
            if self.pseudoRadialBinning == "k2+j2"
                [kj,kr] = ndgrid(self.jWavenumber,self.kRadial);
                Kh = sqrt(kj.^2 + kr.^2);
                allKs = unique(reshape(abs(Kh),[],1),'sorted');
                deltaK = max(diff(allKs));

                kAxis_ = 0:deltaK:(max(allKs)+deltaK/2);
            elseif self.pseudoRadialBinning == "adaptive"
                % n = find(self.kRadial < self.jWavenumber(2),1,"last");
                n = find(self.kRadial > self.jWavenumber(2),1,"first");
                kAxis_ = cat(1,self.kRadial(1:n),self.jWavenumber(3:end));
            end

            kPseudoRadial = reshape(kAxis_,[],1);

            % kPseudoRadial = reshape(kAxis_,[],1);
            %
            % n = find(kPseudoRadial < self.jWavenumber(2),1,"last");
            % kPseudoRadial = cat(1,kPseudoRadial(1:n),self.jWavenumber(2:end));

            % n = find(self.kRadial < self.jWavenumber(2),1,"last");
            % kPseudoRadial = cat(1,self.kRadial(1:n),self.jWavenumber(2:end));

        end

        function omega = get.omegaAxis(self)
            omega=self.wvt.Omega(2:end,:);
            omegaj1=omega(1,:);
            dOmega=max(diff(sort(omegaj1(:))));
            omega=min(omega(:)):dOmega:max(omega(:));
        end

        function kepe = get.kePeAxis(self)
            % fraction = self.geo_hke_jk./(self.geo_hke_jk+self.geo_pe_jk);
            a = cat(1,0,10.^linspace(log10(0.01),log10(0.5),10).');
            kepe = cat(1,a,1-flip(a(1:end-1)));
        end

        function t = get.t_wv(self)
            % Get time vector from the model output file
            %
            % Reads the 't' variable from the wave-vortex file.
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.t_wv(self)
            % - Returns t: time vector from wave-vortex file
            t = self.wvfile.readVariables('wave-vortex/t');
        end

        function forcingNames = get.forcingNames(self)
            % Get names of the forcings
            %
            % Reads the 'forcingNames' from the wave-vortex file and adds
            % in antialiasing if appropriate.
            %
            % - Topic: Dependent property getter
            % - Declaration: t = get.forcingNames(self)
            % - Returns forcingNames: strings
            forcingNames = self.wvt.forcingNames;
            if self.diagfile.hasVariableWithName("Z0_antialias_filter")
                forcingNames(end+1) = "antialias filter";
            end
        end

        function bool = get.diagnosticsHasExplicitAntialiasing(self)
            bool = false;
            if ~isempty(self.diagfile)
                if self.diagfile.hasVariableWithName("Z0_antialias_filter")
                    bool = true;
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % - Topic: Units
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        setEnergyUnits(self, units)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = WVDiagnostics(filename,options)
            % Initializes the WVDiagnostics object, loads the wave-vortex transform and diagnostics files.
            %
            % - Topic: Initialization
            % - Declaration: self = WVDiagnostics(filename,options)
            % - Parameter filename: path to the WVModel output file
            % - Parameter options.diagnosticsFilePath: (optional) path to the diagnostics file
            % - Returns self: the constructed WVDiagnostics object
            arguments
                filename
                options.diagnosticsFilePath
            end
            self.wvpath = filename;


            [self.wvt, self.wvfile] = WVTransform.waveVortexTransformFromFile(filename,iTime=Inf,shouldReadOnly=true);
            self.wvt.addOperation(SpatialForcingOperation(self.wvt));

            [fpath,fname,~] = fileparts(filename);
            if ~isfield(options,"diagnosticsFilePath")
                if ~isempty(fpath)
                    self.diagpath = fullfile(fpath,strcat(fname,"-diagnostics.nc"));
                else
                    self.diagpath= fullfile(pwd,strcat(fname,"-diagnostics.nc"));
                end
            else
                self.diagpath = options.diagnosticsFilePath;
            end
            if exist(self.diagpath,"file")
                self.diagfile = NetCDFFile(self.diagpath);
            else
                warning("No diagnostics file found. Some functionality will not be available.")
            end

            if ~isempty(fpath)
                self.wvaapath = fullfile(fpath,strcat(fname,"-wvt-aa.nc"));
            else
                self.wvaapath= fullfile(pwd,strcat(fname,"-wvt-aa.nc"));
            end

            self.zscale = self.wvt.f^2;
            self.z_flux_scale = self.zscale/(86400*365);

            addlistener(self,'iTime','PostSet',@WVDiagnostics.iTimeChanged);
        end
    end

    methods (Static)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);

            propertyAnnotations(end+1) = CANumericProperty('t_diag',{}, 's', 'time of observations');
            propertyAnnotations(end).attributes('standard_name') = 'time';
            propertyAnnotations(end).attributes('axis') = 'T';

            propertyAnnotations(end+1) = CANumericProperty('t_wv',{}, 's', 'time of observations');
            propertyAnnotations(end).attributes('standard_name') = 'time';
            propertyAnnotations(end).attributes('axis') = 'T';

            propertyAnnotations(end+1) = CADimensionProperty('j', 'mode number', 'vertical mode number');
            propertyAnnotations(end+1) = CADimensionProperty('kRadial', 'rad/m', 'isotropic wavenumber dimension');

            propertyAnnotations(end+1) = CANumericProperty('Lr2',{'j'},'m^2/rad^2', 'squared Rossby radius');
            propertyAnnotations(end+1) = CANumericProperty('Lr2_pm',{'j','kRadial'},'m^2/rad^2', 'squared deformation radius of each wave mode', detailedDescription='- topic: Domain Attributes — Stratification');
            propertyAnnotations(end+1) = CANumericProperty('omega_jk',{'j','kRadial'},'rad s^{-1}', 'intrinsic frequency wave mode', detailedDescription='- topic: Domain Attributes — Stratification');
            propertyAnnotations(end+1) = CANumericProperty('geo_hke_jk',{'j','kRadial'},'m^{-1}', 'multiplicative factor that multiplies $$A_0^2$$ to compute kinetic energy.');
            propertyAnnotations(end+1) = CANumericProperty('geo_pe_jk',{'j','kRadial'},'m^{-1}', 'multiplicative factor that multiplies $$A_0^2$$ to compute potential energy.');

        end
        
        cmap = cmocean(ColormapName,varargin)
        cmap = crameri(ColormapName,varargin)

        [x, t] = CosineTransformBack( f, xbar, varargin )
        [xbar, f] = CosineTransformForward( t, x, varargin )
        [x, t] = SineTransformBack( f, xbar, varargin )
        [xbar, f] = SineTransformForward( t, x, varargin )

        matrix = DCT1(N)
        matrix = iDCT1(N)
        matrix = DST1(N)
        matrix = iDST1(N)

        matrix = DCT2(N)
        matrix = iDCT2(N)
        matrix = DST2(N)
        matrix = iDST2(N)

        Epm = geostrophicGeostrophicWaveEnergy(wvt,mask)
        E0 = waveWaveGeostrophicEnergy(wvt,mask)
        E0 = waveWaveGeostrophicEnergyForMode(wvt,maskKU,maskKUx,Nj)

        cmap = symmetricTintMap(c,options)

        function iTimeChanged(~,eventData)
            self = eventData.AffectedObject;
            self.wvt.initFromNetCDFFile(self.wvfile,iTime=self.iTime);
            if ~isempty(self.wvt_aa_cache)
                self.wvt_aa_cache.t = self.wvt.t;
                [self.wvt_aa_cache.A0,self.wvt_aa_cache.Ap,self.wvt_aa_cache.Am] = self.wvt.spectralVariableWithResolution(self.wvt_aa_cache,self.wvt.A0,self.wvt.Ap,self.wvt.Am);
            end
        end

        function [X,Y,U,V] = PoissonFlowFromFluxWithAxes(x, y, flux)
            % We will treat the first dimension as `x' and the second
            % dimension as `y'. This means that the flux in the usual form,
            % which is j by kRadial, might need to be transposed to get
            % what you want.
            %
            % [X,Y,U,V] = WVDiagnostics.PoissonFlowFromFlux(wvt.kRadial,jWavenumber,flux.');
            % quiver(X,Y,10*U,10*V,'off',Color=0*[1 1 1])
            % For the DCT2/DST2 we use a half-shift grid
            N = length(x);
            M = length(y);
            dk = 1/(2*N*(x(2)-x(1)));
            dl = 1/(2*M*(y(2)-y(1)));

            k=dk*(0:(N-1)).';
            l=dl*(0:(M-1)).';

            DCTx = WVDiagnostics.DCT2(N);
            DCTy = WVDiagnostics.DCT2(M);

            flux_ky = DCTx*flux;
            flux_kl = shiftdim(DCTy*shiftdim(flux_ky,1),1);

            [K,L] = ndgrid(k,l);

            D = -((2*pi*K).^2 + (2*pi*L).^2);
            D(1,1) = Inf;

            UFactor = 2*pi*K./D;
            VFactor = 2*pi*L./D;

            iDCTx = WVDiagnostics.iDCT2(N);
            iDSTy = WVDiagnostics.iDST2(M);
            iDSTy = circshift(iDSTy,1,2);
            iDSTy(:,1) = 0;

            V_xl = iDCTx*(VFactor.*flux_kl);
            V = shiftdim(iDSTy*shiftdim(V_xl,1),1);

            iDCTy = WVDiagnostics.iDCT2(M);
            iDSTx = WVDiagnostics.iDST2(N);
            iDSTx = circshift(iDSTx,1,2);
            iDSTx(:,1) = 0;

            U_xl = iDSTx*(UFactor.*flux_kl);
            U = shiftdim(iDCTy*shiftdim(U_xl,1),1);

            [X,Y] = ndgrid(x,y);
        end

        function [X,Y,U,V] = PoissonFlowFromFluxDCTI(x,y,flux)
            % We will treat the first dimension as `x' and the second
            % dimension as `y'. This means that the flux in the usual form,
            % which is j by kRadial, might need to be transposed to get
            % what you want.
            %
            % [X,Y,U,V] = WVDiagnostics.PoissonFlowFromFlux(wvt.kRadial,jWavenumber,flux.');
            % quiver(X,Y,10*U,10*V,'off',Color=0*[1 1 1])

            [X,Y] = ndgrid(x,y);
            [flux_bar, f_alpha] = WVDiagnostics.CosineTransformForward( x, flux, 1 );
            [flux_bar2, f_beta] = WVDiagnostics.CosineTransformForward( y, flux_bar, 2 );
            [ALPHA,BETA] = ndgrid(f_alpha,f_beta);
            D = -((2*pi*ALPHA).^2 + (2*pi*BETA).^2);
            D(1,1) = Inf;
            UFactor = 2*pi*ALPHA./D;
            VFactor = 2*pi*BETA./D;
            tmp = WVDiagnostics.CosineTransformBack(f_beta,UFactor.*flux_bar2,2);
            U = WVDiagnostics.SineTransformBack(f_alpha(2:end-1,:),tmp(2:end-1,:),1);
            V = WVDiagnostics.CosineTransformBack(f_alpha,WVDiagnostics.SineTransformBack(f_beta(2:end-1),VFactor(:,2:end-1).*flux_bar2(:,2:end-1),2),1);
        end

        function bool = areEnergyReservoirsComplete(reservoirs)
            bool = all(sum([reservoirs.vectorContents],2)==1);
            if ~bool
                warning('The collection of energy reservoirs is either not complete or over complete')
            end
        end

        function bool = areTriadComponentsComplete(reservoirs)
            bool = all(sum([reservoirs.vectorContents],2)==1);
            if ~bool
                warning('The collection of triad components is either not complete or over complete')
            end
        end

        function fancyName = fancyNameForName(name)
            switch name
                case "ke_g"
                    fancyName = "geostrophic kinetic";
                case "pe_g"
                    fancyName = "geostrophic potential";
                case "te_g"
                    fancyName = "geostrophic";
                case "te_mda"
                    fancyName = "mean density anomaly";
                case "te_gmda"
                    fancyName = "geostrophic + mda";
                case "te_igw"
                    fancyName = "internal gravity wave";
                case "te_io"
                    fancyName = "inertial";
                case "te_wave"
                    fancyName = "wave";
                case "te"
                    fancyName = "total";
                case "te_quadratic"
                    fancyName = "total quadratic";
                case "te_damp"
                    fancyName = "closure region";
                otherwise
                    error("unknown energy reservoir");
            end
        end

        function filtered_fluxes = filterFluxesForReservoir(fluxes,options)
            arguments
                fluxes
                options.filter = @(v) v;
            end

            reservoirNames = string(fields(fluxes));
            reservoirNames(reservoirNames=="name") = [];
            reservoirNames(reservoirNames=="fancyName") = [];
            for iForce=1:length(fluxes)
                for iReservoir = 1:length(reservoirNames)
                    fluxes(iForce).(reservoirNames(iReservoir)) = options.filter(fluxes(iForce).(reservoirNames(iReservoir)));
                end
            end
            filtered_fluxes = fluxes;
        end

        function v = version()
            % Locate folder containing this class file
            classFolder = fileparts(mfilename('fullpath'));
            packageRoot = fileparts(classFolder);

            % Path to mpackage.json
            jsonFile = fullfile(packageRoot, 'resources', 'mpackage.json');

            % Read + decode JSON
            info = jsondecode(fileread(jsonFile));

            % Extract version field
            v = info.version;
        end

    end
end