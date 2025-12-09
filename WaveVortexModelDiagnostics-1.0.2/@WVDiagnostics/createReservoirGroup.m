function createReservoirGroup(self,options)
% Create a reservoir group in the diagnostics NetCDF and populate it.
%
% Create a reservoir group in the diagnostics NetCDF and populate it.
% Create (or open) a diagnostics group containing reservoir definitions as defined by specified flow components.
% The function will create variables and dimensions as needed, then loop over
% the requested time indices to compute and write reservoir diagnostics.
%
% - Topic: Diagnostics Generation
% - Declaration: createReservoirGroup(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter outputfile: NetCDFGroup to write into (default: self.diagfile).
% - Parameter name: (optional) Name of the reservoir group (default: "reservoir-damped-wave-geo").
% - Parameter flowComponents: WVFlowComponent array specifying reservoirs to create.
% - Parameter timeIndices: Time indices to process (default: all times in diagnostics file).
arguments
    self WVDiagnostics
    options.outputfile NetCDFGroup
    options.name string = "reservoir-damped-wave-geo"
    options.flowComponents WVFlowComponent
    options.timeIndices
end

if self.diagnosticsHasExplicitAntialiasing
    wvt = self.wvt_aa;
else
    wvt = self.wvt;
end

if ~isfield(options,"outputfile")
    options.outputfile = self.diagfile;
end

if ~isfield(options,"timeIndices")
    t = self.diagfile.readVariables("t");
    timeIndices = 1:length(t);
else
    timeIndices = options.timeIndices;
end

if isfield(options,"flowComponents")
    flowComponents=options.flowComponents;
else
    svv = wvt.forcingWithName("adaptive damping");
    NoDampMask = (wvt.Kh < (svv.k_damp+svv.k_no_damp)/2) & (wvt.J < (svv.j_damp + svv.j_no_damp)/2);

    i = 1;
    flowComponents(i) = wvt.flowComponentWithName("wave") + wvt.flowComponentWithName("inertial");
    flowComponents(i).name = "wave";
    flowComponents(i).maskAp = flowComponents(i).maskAp & NoDampMask;
    flowComponents(i).maskAm = flowComponents(i).maskAm & NoDampMask;

    i = i+1;
    flowComponents(i) = wvt.flowComponentWithName("geostrophic") + wvt.flowComponentWithName("mda");
    flowComponents(i).name = "geostrophic";
    flowComponents(i).maskA0 = flowComponents(i).maskA0 & NoDampMask;

    i = i+1;
    flowComponents(i) = wvt.flowComponentWithName("wave") + wvt.flowComponentWithName("inertial");
    flowComponents(i).name = "damped wave";
    flowComponents(i).maskAp = flowComponents(i).maskAp & ~NoDampMask;
    flowComponents(i).maskAm = flowComponents(i).maskAm & ~NoDampMask;

    i = i+1;
    flowComponents(i) = wvt.flowComponentWithName("geostrophic") + wvt.flowComponentWithName("mda");
    flowComponents(i).name = "damped geostrophic";
    flowComponents(i).maskA0 = flowComponents(i).maskA0 & ~NoDampMask;
end

[triadVar, forcingVar, energyVar] = self.variablesForReservoirGroup(outputfile=options.outputfile,name=options.name,flowComponents=flowComponents);
forcingNames = wvt.forcingNames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Loop over the the requested time indices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

integrationLastInformWallTime = datetime('now');
loopStartTime = integrationLastInformWallTime;
integrationLastInformLoopNumber = 1;
integrationInformTime = 10;
fprintf("Starting loop to compute reservoir fluxes for %d time indices.\n",length(timeIndices));
for timeIndex = 1:length(timeIndices)
    deltaWallTime = datetime('now')-integrationLastInformWallTime;
    if ( seconds(deltaWallTime) > integrationInformTime)
        wallTimePerLoopTime = deltaWallTime / (timeIndex - integrationLastInformLoopNumber);
        wallTimeRemaining = wallTimePerLoopTime*(length(timeIndices) - timeIndex + 1);
        fprintf('Time index %d of %d. Estimated time to finish is %s (%s)\n', timeIndex, length(timeIndices), wallTimeRemaining, datetime(datetime('now')+wallTimeRemaining,TimeZone='local',Format='d-MMM-y HH:mm:ss Z')) ;
        integrationLastInformWallTime = datetime('now');
        integrationLastInformLoopNumber = timeIndex;
    end

    outputIndex = timeIndices(timeIndex);
    self.iTime = timeIndices(timeIndex);

    self.addTriadFluxesForReservoirGroupAtTime(triadVar=triadVar,flowComponents=flowComponents,wvt=wvt,outputIndex=outputIndex);

    F = wvt.fluxForForcing();
    for i=1:length(forcingNames)
        [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(F{forcingNames(i)}.Fp,F{forcingNames(i)}.Fm,F{forcingNames(i)}.F0);
        for k=1:length(flowComponents)
            dE = sum(flowComponents(k).maskAp(:).*Ep(:) + flowComponents(k).maskAm(:).*Em(:) + flowComponents(k).maskA0(:).*E0(:));
            forcingVar{forcingNames(i)+"_"+k}.setValueAlongDimensionAtIndex(dE,'t',outputIndex);
        end
    end

    for k=1:length(flowComponents)
        dE = wvt.totalEnergyOfFlowComponent(flowComponents(k));
        energyVar{"E_"+k}.setValueAlongDimensionAtIndex(dE,'t',outputIndex);
    end
end
deltaLoopTime = datetime('now')-loopStartTime;
fprintf("Total loop time %s, which is %s per time index.\n",deltaLoopTime,deltaLoopTime/length(timeIndices));

end