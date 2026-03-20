function createDiagnosticsFile(self,options)
% Create a new diagnostics file and compute diagnostics from WVModel output.
%
% Create a new diagnostics file and compute diagnostics from WVModel output.
% Initializes a diagnostics NetCDF file, computes and stores energy, enstrophy, APE, APV, and flux diagnostics for each time step.
%
% - Topic: Diagnostics Generation
% - Declaration: createDiagnosticsFile(self,options)
% - Parameter self: WVDiagnostics object.
% - Parameter stride: (optional) Stride for time sampling (default: 1).
% - Parameter timeIndices: Indices of time steps to process.
% - Parameter filename: Output path for the diagnostics file.
% - Parameter shouldMeasureAntialiasingFlux: (optional) (optional, logical) If true, computes antialiasing flux diagnostics (default: false).
% - Parameter shouldUseHigherOrderFlux: (optional) (optional, logical) If true, uses higher-order flux calculation (default: false).
arguments
    self WVDiagnostics
    options.stride = 1
    options.timeIndices
    options.filename
    options.shouldMeasureAntialiasingFlux logical = false
    options.shouldUseHigherOrderFlux logical = false
end

if ~(isa(self.wvt,"WVTransformHydrostatic") || isa(self.wvt,"WVTransformBoussinesq") || isa(self.wvt,"WVTransformConstantStratification"))
    error("Your transformation type is not supported.");
end

if exist(self.diagpath,"file") && ~isfield(options,"filename")
    [found, idx] = ismember(self.t_diag, self.t_wv);
    if ~all(found)
        error('Some entries of t_diag are not exactly in t_wv.');
    end
    didx = diff(idx);
    stride = didx(1);
    timeIndices = (idx(end)+stride):stride:length(self.t_wv);
    wvt = WVTransform.waveVortexTransformFromFile(self.wvfile.path,iTime=Inf);
    diagfile = self.diagfile;
    ncfile = self.wvfile;

    outputIndexOffset = length(self.t_diag);

    wvt.addOperation(SpatialForcingOperation(wvt));
    int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);
    
    % 1. Measures of energy, APV and enstrophy
    EnergyByComponent = configureDictionary("string","cell");
    flowComponentNames = wvt.flowComponentNames;
    for i=1:length(flowComponentNames)
        comp = wvt.flowComponentWithName(flowComponentNames(i));
        name = "E_" + comp.abbreviatedName;
        EnergyByComponent{name} = diagfile.variableWithName(name);
    end
    EnergyByComponent{"KE_g"} = diagfile.variableWithName("KE_g");
    EnergyByComponent{"PE_g"} = diagfile.variableWithName("PE_g");

    variable_ke = diagfile.variableWithName("ke");
    variable_pe = diagfile.variableWithName("pe_quadratic");
    variable_ape = diagfile.variableWithName("ape");
    variable_ape_bg = diagfile.variableWithName("ape_bg");
    variable_z = diagfile.variableWithName("enstrophy_quadratic");
    variable_apv = diagfile.variableWithName("enstrophy_apv");

    % 2. Forcing fluxes
    EnergyFlux = configureDictionary("string","cell");
    EnstrophyFlux = configureDictionary("string","cell");
    EnergyFluxTrue = configureDictionary("string","cell");
    EnstrophyFluxTrue = configureDictionary("string","cell");
    forcingNames = wvt.forcingNames;
    for i=1:length(forcingNames)
        name = replace(forcingNames(i),"-","_");
        name = replace(name," ","_");
        EnergyFlux{forcingNames(i)}.Ep = diagfile.variableWithName("Ep_" + name);
        EnergyFlux{forcingNames(i)}.Em = diagfile.variableWithName("Em_" + name);
        EnergyFlux{forcingNames(i)}.KE0 = diagfile.variableWithName("KE0_" + name);
        EnergyFlux{forcingNames(i)}.PE0 = diagfile.variableWithName("PE0_" + name);
        EnstrophyFlux{forcingNames(i)} = diagfile.variableWithName("Z0_" + name);
        EnergyFluxTrue{forcingNames(i)} = diagfile.variableWithName("E_" + name);
        EnstrophyFluxTrue{forcingNames(i)} = diagfile.variableWithName("Z_" + name);
    end

    % 3. Triads
    triadFlowComponents = [wvt.flowComponentWithName('wave'); wvt.flowComponentWithName('inertial'); wvt.flowComponentWithName('geostrophic'); wvt.flowComponentWithName('mda')];
    EnergyTriads = configureDictionary("string","cell");
    EnstrophyTriads = configureDictionary("string","cell");
    for i=1:length(triadFlowComponents)
        for j=1:length(triadFlowComponents)
            name = triadFlowComponents(i).abbreviatedName + "_" + triadFlowComponents(j).abbreviatedName;
            EnergyTriads{name}.Ep = diagfile.variableWithName("Ep_" + name);
            EnergyTriads{name}.Em = diagfile.variableWithName("Em_" + name);
            EnergyTriads{name}.KE0 = diagfile.variableWithName("KE0_" + name);
            EnergyTriads{name}.PE0 = diagfile.variableWithName("PE0_" + name);
            EnstrophyTriads{name} = diagfile.variableWithName("Z0_" + name);
        end
    end
else
    outputIndexOffset = 0;
    if options.shouldMeasureAntialiasingFlux
        [wvt_lowres, ncfile] = WVTransform.waveVortexTransformFromFile(self.wvfile.path,iTime=Inf,shouldReadOnly=true);
        if isa(wvt_lowres,"WVTransformBoussinesq")
            if exist(self.wvaapath,"file")
                wvt = WVTransform.waveVortexTransformFromFile(self.wvaapath,iTime=Inf,shouldReadOnly=true);
            else
                wvt = wvt_lowres.waveVortexTransformWithExplicitAntialiasing();
                wvt.writeToFile(self.wvaapath);
            end
        else
            wvt = wvt_lowres.waveVortexTransformWithExplicitAntialiasing();
        end
    else
        [wvt, ncfile] = WVTransform.waveVortexTransformFromFile(self.wvfile.path,iTime=Inf,shouldReadOnly=true);
    end
    if ncfile.hasVariableWithName('wave-vortex/t')
        tDim = ncfile.readVariables('wave-vortex/t');
    else
        tDim = ncfile.readVariables('t');
    end

    if ~isfield(options,"timeIndices")
        timeIndices = 1:options.stride:length(tDim);
    else
        timeIndices = options.timeIndices;
    end

    wvt.addOperation(SpatialForcingOperation(wvt));
    int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);

    %% setup diagnostic output file
    if isfield(options,"filename")
        fpath = fileparts(self.diagpath);
        if ~isempty(fpath)
            filepath = fullfile(fpath,options.filename);
        else
            filepath= fullfile(pwd,options.filename);
        end
        disp("writing diagnostics file: " + filepath);
        diagfile = NetCDFFile(filepath);
    else
        diagfile = NetCDFFile(self.diagpath);
    end
    self.diagfile = diagfile;

    dimensionNames = ["j","kRadial"];
    for iDim=1:length(dimensionNames)
        dimAnnotation = wvt.dimensionAnnotationWithName(dimensionNames(iDim));
        dimAnnotation.attributes('units') = dimAnnotation.units;
        dimAnnotation.attributes('long_name') = dimAnnotation.description;
        diagfile.addDimension(dimAnnotation.name,wvt.(dimAnnotation.name),attributes=dimAnnotation.attributes);
    end

    varAnnotation = wvt.propertyAnnotationWithName('Lr2');
    diagfile.addVariable(varAnnotation.name,varAnnotation.dimensions,wvt.Lr2,isComplex=varAnnotation.isComplex,attributes=varAnnotation.attributes );

    if isa(wvt,"WVTransformHydrostatic")
        h_pm = repmat(wvt.h_pm,[1 wvt.Nkl]);
    else
        h_pm = wvt.h_pm;
    end

    [omegaN,n,hke_jk,pe_jk,hN] = wvt.transformToRadialWavenumber(abs(wvt.Omega),ones(size(wvt.Omega)),wvt.A0_KE_factor,wvt.A0_PE_factor,h_pm);
    omegaJK = (omegaN./n);
    h_kj = (hN./n);
    Lr2_pm = wvt.g * h_kj / (wvt.f*wvt.f);
    diagfile.addVariable("omega_jk",dimensionNames,omegaJK,isComplex=false);
    diagfile.addVariable("geo_hke_jk",dimensionNames,hke_jk,isComplex=false);
    diagfile.addVariable("geo_pe_jk",dimensionNames,pe_jk,isComplex=false);
    diagfile.addVariable("Lr2_pm",dimensionNames,Lr2_pm,isComplex=false);

    varAnnotation = wvt.propertyAnnotationWithName('t');
    varAnnotation.attributes('units') = varAnnotation.units;
    varAnnotation.attributes('long_name') = varAnnotation.description;
    varAnnotation.attributes('standard_name') = 'time';
    varAnnotation.attributes('long_name') = 'time';
    varAnnotation.attributes('units') = 'seconds since 1970-01-01 00:00:00';
    varAnnotation.attributes('axis') = 'T';
    varAnnotation.attributes('calendar') = 'standard';
    diagfile.addDimension(varAnnotation.name,length=Inf,type="double",attributes=varAnnotation.attributes);

    % 1. Measures of energy, APV and enstrophy
    EnergyByComponent = configureDictionary("string","cell");
    flowComponentNames = wvt.flowComponentNames;
    for i=1:length(flowComponentNames)
        comp = wvt.flowComponentWithName(flowComponentNames(i));
        name = "E_" + comp.abbreviatedName;
        EnergyByComponent{name} = diagfile.addVariable(name,"t",type="double",isComplex=false);
    end
    EnergyByComponent{"KE_g"} = diagfile.addVariable("KE_g","t",type="double",isComplex=false);
    EnergyByComponent{"PE_g"} = diagfile.addVariable("PE_g","t",type="double",isComplex=false);
    variable_ke = diagfile.addVariable("ke","t",type="double",isComplex=false);
    variable_pe = diagfile.addVariable("pe_quadratic","t",type="double",isComplex=false);
    variable_ape = diagfile.addVariable("ape","t",type="double",isComplex=false);
    variable_ape_bg = diagfile.addVariable("ape_bg","t",type="double",isComplex=false);
    variable_z = diagfile.addVariable("enstrophy_quadratic","t",type="double",isComplex=false);
    variable_apv = diagfile.addVariable("enstrophy_apv","t",type="double",isComplex=false);

    % 2. Forcing fluxes
    EnergyFlux = configureDictionary("string","cell");
    EnstrophyFlux = configureDictionary("string","cell");
    EnergyFluxTrue = configureDictionary("string","cell");
    EnstrophyFluxTrue = configureDictionary("string","cell");
    forcingNames = wvt.forcingNames;
    dimensionNames = ["j", "kRadial", "t"];
    for i=1:length(forcingNames)
        name = replace(forcingNames(i),"-","_");
        name = replace(name," ","_");
        EnergyFlux{forcingNames(i)}.Ep = diagfile.addVariable("Ep_" + name,dimensionNames,type="double",isComplex=false);
        EnergyFlux{forcingNames(i)}.Em = diagfile.addVariable("Em_" + name,dimensionNames,type="double",isComplex=false);
        EnergyFlux{forcingNames(i)}.KE0 = diagfile.addVariable("KE0_" + name,dimensionNames,type="double",isComplex=false);
        EnergyFlux{forcingNames(i)}.PE0 = diagfile.addVariable("PE0_" + name,dimensionNames,type="double",isComplex=false);
        EnstrophyFlux{forcingNames(i)} = diagfile.addVariable("Z0_" + name,dimensionNames,type="double",isComplex=false);
        EnergyFluxTrue{forcingNames(i)} = diagfile.addVariable("E_" + name,dimensionNames,type="double",isComplex=false);
        EnstrophyFluxTrue{forcingNames(i)} = diagfile.addVariable("Z_" + name,dimensionNames,type="double",isComplex=false);
    end

    % 3. Triads
    triadFlowComponents = [wvt.flowComponentWithName('wave'); wvt.flowComponentWithName('inertial'); wvt.flowComponentWithName('geostrophic'); wvt.flowComponentWithName('mda')];
    EnergyTriads = configureDictionary("string","cell");
    EnstrophyTriads = configureDictionary("string","cell");
    dimensionNames = ["j", "kRadial", "t"];
    for i=1:length(triadFlowComponents)
        for j=1:length(triadFlowComponents)
            name = triadFlowComponents(i).abbreviatedName + "_" + triadFlowComponents(j).abbreviatedName;
            EnergyTriads{name}.Ep = diagfile.addVariable("Ep_" + name,dimensionNames,type="double",isComplex=false);
            EnergyTriads{name}.Em = diagfile.addVariable("Em_" + name,dimensionNames,type="double",isComplex=false);
            EnergyTriads{name}.KE0 = diagfile.addVariable("KE0_" + name,dimensionNames,type="double",isComplex=false);
            EnergyTriads{name}.PE0 = diagfile.addVariable("PE0_" + name,dimensionNames,type="double",isComplex=false);
            EnstrophyTriads{name} = diagfile.addVariable("Z0_" + name,dimensionNames,type="double",isComplex=false);
        end
    end
end
%% loop through time computing diagnostics
integrationLastInformWallTime = datetime('now');
integrationLastInformLoopNumber = 1;
integrationInformTime = 10;
disp("Starting loop");
for timeIndex = 1:length(timeIndices)
    deltaWallTime = datetime('now')-integrationLastInformWallTime;
    if ( seconds(deltaWallTime) > integrationInformTime)
        wallTimePerLoopTime = deltaWallTime / (timeIndex - integrationLastInformLoopNumber);
        wallTimeRemaining = wallTimePerLoopTime*(length(timeIndices) - timeIndex);
        fprintf('Time index %d of %d. Estimated time to finish is %s (%s)\n', timeIndex, length(timeIndices), wallTimeRemaining, datetime(datetime('now')+wallTimeRemaining,TimeZone='local',Format='d-MMM-y HH:mm:ss Z')) ;
        integrationLastInformWallTime = datetime('now');
        integrationLastInformLoopNumber = timeIndex;
    end

    iTime = timeIndices(timeIndex);
    outputIndex = timeIndex + outputIndexOffset;
    if options.shouldMeasureAntialiasingFlux
        wvt_lowres.initFromNetCDFFile(ncfile,iTime=iTime);
        wvt.t = wvt_lowres.t;
        [wvt.A0,wvt.Ap,wvt.Am] = wvt_lowres.spectralVariableWithResolution(wvt,wvt_lowres.A0,wvt_lowres.Ap,wvt_lowres.Am);
    else
        wvt.initFromNetCDFFile(ncfile,iTime=iTime);
    end

    diagfile.variableWithName('t').setValueAlongDimensionAtIndex(wvt.t,'t',outputIndex);

    rho_nm = wvt.chebfunForZArray(wvt.rho_nm);
    N2 = (-wvt.g/wvt.rho0)*diff(rho_nm);
    p_nm = (-wvt.g/wvt.rho0)*cumsum(rho_nm-wvt.rho0);
    p_nm = p_nm - p_nm(0);

    % 1. Measures of energy, APV and enstrophy
    for i=1:length(flowComponentNames)
        comp = wvt.flowComponentWithName(flowComponentNames(i));
        name = "E_" + comp.abbreviatedName;
        E = wvt.totalEnergyOfFlowComponent(comp);
        EnergyByComponent{name}.setValueAlongDimensionAtIndex(E,'t',outputIndex);
    end
    EnergyByComponent{"KE_g"}.setValueAlongDimensionAtIndex(wvt.geostrophicKineticEnergy,'t',outputIndex);
    EnergyByComponent{"PE_g"}.setValueAlongDimensionAtIndex(wvt.geostrophicPotentialEnergy,'t',outputIndex);
    if wvt.isHydrostatic
        [u,v,ape,apv] = wvt.variableWithName('u','v','ape','apv');
        ke = (u.^2 + v.^2)/2;
    else
        [u,v,w,ape,apv] = wvt.variableWithName('u','v','w','ape','apv');
        ke = (u.^2 + v.^2 + w.^2)/2;
    end
    variable_ke.setValueAlongDimensionAtIndex(int_vol(ke),'t',outputIndex);
    variable_pe.setValueAlongDimensionAtIndex(int_vol(shiftdim(wvt.N2,-2).*wvt.eta.*wvt.eta/2),'t',outputIndex);
    variable_ape.setValueAlongDimensionAtIndex(int_vol(ape),'t',outputIndex);
    variable_ape_bg.setValueAlongDimensionAtIndex(int_vol(p_nm(wvt.Z-wvt.eta_true)),'t',outputIndex);
    variable_z.setValueAlongDimensionAtIndex(wvt.totalEnstrophy,'t',outputIndex);
    variable_apv.setValueAlongDimensionAtIndex(0.5*int_vol(apv.*apv),'t',outputIndex);

    % 2. Forcing fluxes

    F = wvt.fluxForForcing();
    for i=1:length(forcingNames)
        [Ep,Em,KE0,PE0] = wvt.energyFluxFromNonlinearFlux(F{forcingNames(i)}.Fp,F{forcingNames(i)}.Fm,F{forcingNames(i)}.F0);
        [Ep_jk,Em_jk,KE0_jk,PE0_jk] = wvt.transformToRadialWavenumber(Ep,Em,KE0,PE0);
        if wvt.isHydrostatic
            [Fu,Fv,Feta] = wvt.spatialFluxForForcingWithName(forcingNames(i));
            % F_density = wvt.u .* Fu + wvt.v .* Fv+ wvt.eta_true .* shiftdim(wvt.N2,-2) .* Feta;

            S_energy = wvt.crossSpectrumWithFgTransform(wvt.u,Fu);
            S_energy = S_energy + wvt.crossSpectrumWithFgTransform(wvt.v,Fv);
            S_energy = S_energy + wvt.crossSpectrumWithGgTransform(wvt.eta_true,Feta);

            if forcingNames(i) == "nonlinear advection"
                Fu = Fu + wvt.f*wvt.v;
                Fv = Fv - wvt.f*wvt.u;
            end
            DF_x =  - wvt.diffZF(Fv); % w_y - v_z
            DF_y = wvt.diffZF(Fu);  % u_z - w_x
            DF_z = wvt.diffX(Fv) - wvt.diffY(Fu);  % v_x - u_y
        else
            [Fu,Fv,Fw,Feta] = wvt.spatialFluxForForcingWithName(forcingNames(i));
            % F_density = wvt.u .* Fu + wvt.v .* Fv +  wvt.w .* Fw + wvt.eta_true .* shiftdim(wvt.N2,-2) .* Feta;

            S_energy = wvt.crossSpectrumWithFgTransform(wvt.u,Fu);
            S_energy = S_energy + wvt.crossSpectrumWithFgTransform(wvt.v,Fv);
            S_energy = S_energy + wvt.crossSpectrumWithGgTransform(wvt.w./wvt.N2Function(wvt.Z),Fw);
            S_energy = S_energy + wvt.crossSpectrumWithGgTransform(wvt.eta_true,Feta);

            if forcingNames(i) == "nonlinear advection"
                Fu = Fu + wvt.f*wvt.v;
                Fv = Fv - wvt.f*wvt.u;
            end
            DF_x = wvt.diffY(Fw) - wvt.diffZF(Fv); % w_y - v_z
            DF_y = wvt.diffZF(Fu) - wvt.diffX(Fw);  % u_z - w_x
            DF_z = wvt.diffX(Fv) - wvt.diffY(Fu);  % v_x - u_y
        end

        if forcingNames(i) == "nonlinear advection"
            % F_density = F_density + wvt.w .* shiftdim(wvt.N2,-2) .* (wvt.eta_true-wvt.eta);
            S_energy = S_energy + wvt.crossSpectrumWithGgTransform(wvt.w,wvt.eta_true-wvt.eta);
            G_eta = (wvt.N2Function(wvt.Z)./N2(wvt.Z - wvt.eta_true)).*(Feta + wvt.w);
        else
            G_eta = (wvt.N2Function(wvt.Z)./N2(wvt.Z - wvt.eta_true)).*Feta;
        end

        FZ_L = wvt.diffX(Fv) - wvt.diffY(Fu) - wvt.f*wvt.diffZG(G_eta);
        FZ_NL = - wvt.zeta_x .* wvt.diffX(G_eta) - wvt.zeta_y .* wvt.diffY(G_eta)- wvt.zeta_z .* wvt.diffZG(G_eta);
        FZ_NL = FZ_NL - wvt.diffX(wvt.eta_true) .* DF_x - wvt.diffY(wvt.eta_true) .* DF_y - wvt.diffZG(wvt.eta_true) .* DF_z;
        FZ = FZ_L + FZ_NL;
        Z_density = wvt.crossSpectrumWithFgTransform(wvt.apv, FZ);
        Z_jk = wvt.transformToRadialWavenumber( Z_density);

        S_jk = wvt.transformToRadialWavenumber( S_energy);
        EnergyFluxTrue{forcingNames(i)}.setValueAlongDimensionAtIndex(S_jk,'t',outputIndex);
        EnstrophyFluxTrue{forcingNames(i)}.setValueAlongDimensionAtIndex(Z_jk,'t',outputIndex);

        EnergyFlux{forcingNames(i)}.Ep.setValueAlongDimensionAtIndex(Ep_jk,'t',outputIndex);
        EnergyFlux{forcingNames(i)}.Em.setValueAlongDimensionAtIndex(Em_jk,'t',outputIndex);
        EnergyFlux{forcingNames(i)}.KE0.setValueAlongDimensionAtIndex(KE0_jk,'t',outputIndex);
        EnergyFlux{forcingNames(i)}.PE0.setValueAlongDimensionAtIndex(PE0_jk,'t',outputIndex);

        Z = 2*wvt.A0_TZ_factor.*real( F{forcingNames(i)}.F0 .* conj(wvt.A0) );
        EnstrophyFlux{forcingNames(i)}.setValueAlongDimensionAtIndex(wvt.transformToRadialWavenumber(Z),'t',outputIndex);
    end

    % 3. Triads
    for i=1:length(triadFlowComponents)
        for j=1:length(triadFlowComponents)
            if options.shouldUseHigherOrderFlux
                [Fp,Fm,F0] = wvt.rk4NonlinearFluxForFlowComponents(triadFlowComponents(i),triadFlowComponents(j));
            else
                [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(triadFlowComponents(i),triadFlowComponents(j));
            end
            [Ep,Em,KE0,PE0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
            Z0 = wvt.enstrophyFluxFromNonlinearFlux(F0);
            [Ep_jk,Em_jk,KE0_jk,PE0_jk,Z0_jk] = wvt.transformToRadialWavenumber(Ep,Em,KE0,PE0,Z0);

            name = triadFlowComponents(i).abbreviatedName + "_" + triadFlowComponents(j).abbreviatedName;
            EnergyTriads{name}.Ep.setValueAlongDimensionAtIndex(Ep_jk,'t',outputIndex);
            EnergyTriads{name}.Em.setValueAlongDimensionAtIndex(Em_jk,'t',outputIndex);
            EnergyTriads{name}.KE0.setValueAlongDimensionAtIndex(KE0_jk,'t',outputIndex);
            EnergyTriads{name}.PE0.setValueAlongDimensionAtIndex(PE0_jk,'t',outputIndex);
            EnstrophyTriads{name}.setValueAlongDimensionAtIndex(Z0_jk,'t',outputIndex);
        end
    end
end
fprintf("\n");
end