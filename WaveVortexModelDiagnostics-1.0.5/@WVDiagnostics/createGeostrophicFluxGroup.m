function createGeostrophicFluxGroup(self,options)
% Create a geostrophic flux group in the diagnostics NetCDF and populate it.
%
% Computes the geostrophic fluxes from geostrophicFlux, valid only for
% constant stratification..
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

if ~options.outputfile.hasVariableWithName("geostrophic-flux/ggg")
    group = options.outputfile.addGroup("geostrophic-flux");
    dimensionNames = ["j", "kRadial", "t"];

    ggg = group.addVariable("ggg",dimensionNames,type="double",isComplex=false);
    ggw = group.addVariable("ggw",dimensionNames,type="double",isComplex=false);
    ggw_tx = group.addVariable("ggw_tx",dimensionNames,type="double",isComplex=false);
    wwg_tx = group.addVariable("wwg_tx",dimensionNames,type="double",isComplex=false);
else
    ggg = options.outputfile.variableWithName("geostrophic-flux/ggg");
    ggw = options.outputfile.variableWithName("geostrophic-flux/ggw");
    ggw_tx = options.outputfile.variableWithName("geostrophic-flux/ggw_tx");
    wwg_tx = options.outputfile.variableWithName("geostrophic-flux/wwg_tx");
end

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

    [~, spectralFlux] = self.geostrophicFlux();
    ggg.setValueAlongDimensionAtIndex(wvt.transformToRadialWavenumber(spectralFlux.ggg),'t',outputIndex);
    ggw.setValueAlongDimensionAtIndex(wvt.transformToRadialWavenumber(spectralFlux.ggw),'t',outputIndex);
    ggw_tx.setValueAlongDimensionAtIndex(wvt.transformToRadialWavenumber(spectralFlux.ggw_tx),'t',outputIndex);
    wwg_tx.setValueAlongDimensionAtIndex(wvt.transformToRadialWavenumber(spectralFlux.wwg_tx),'t',outputIndex);
end
deltaLoopTime = datetime('now')-loopStartTime;
fprintf("Total loop time %s, which is %s per time index.\n",deltaLoopTime,deltaLoopTime/length(timeIndices));

end