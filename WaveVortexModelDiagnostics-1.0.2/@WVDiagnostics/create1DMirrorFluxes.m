function create1DMirrorFluxes(self,options)
% Create 1D mirror flux diagnostics and write them to the diagnostics NetCDF.
%
% Create 1D mirror flux diagnostics and write them to the diagnostics NetCDF.
% Adds 1D mirror-flux variables (kp, omegaAxis, kePeAxis) to an existing diagnostics
% NetCDF file and computes mirror flux summaries (e.g. F_wwg_kp, pi_w_wwg_kp,
% F_ggw_kp, pi_g_ggw_kp, and the corresponding omega/kePe projections) from the
% current WVTransform for the requested time indices. The function will create
% required dimensions/variables if they do not already exist and then populate
% the variables for each requested time index.
%
% - Topic: Diagnostics Generation
% - Declaration: create1DMirrorFluxes(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter stride: (optional) Stride for time sampling (default: 1).
% - Parameter timeIndices: Indices of time steps to process (default: all times in diagnostics file).
arguments
    self WVDiagnostics
    options.stride = 1
    options.timeIndices
end

if ~exist(self.diagpath,"file")
    error("No existing diagnostics file found.");
end
diagfile = self.diagfile;
wvt = self.wvt;

[kp,bins_0,bins_pm] = self.sparsePseudoRadialAxis;
[omegaAxis,bins_omega] = self.sparseOmegaAxis;
[kePeAxis,bins_kepe] = self.sparseKePeAxis;

if ~isfield(options,"timeIndices")
    t = diagfile.readVariables("t");
    timeIndices = 1:length(t);
else
    timeIndices = options.timeIndices;
end

if ~diagfile.hasVariableWithName("kp")
    diagfile.addDimension("kp",kp);
    dimensionNames = ["kp", "t"];
    pi_w_wwg_kp = diagfile.addVariable("pi_w_wwg_kp",dimensionNames,type="double",isComplex=false);
    F_wwg_kp = diagfile.addVariable("F_wwg_kp",dimensionNames,type="double",isComplex=false);
    pi_g_ggw_kp = diagfile.addVariable("pi_g_ggw_kp",dimensionNames,type="double",isComplex=false);
    F_ggw_kp = diagfile.addVariable("F_ggw_kp",dimensionNames,type="double",isComplex=false);

    diagfile.addDimension("omegaAxis",omegaAxis);
    dimensionNames = ["omegaAxis", "t"];
    pi_w_wwg_omega = diagfile.addVariable("pi_w_wwg_omega",dimensionNames,type="double",isComplex=false);

    diagfile.addDimension("kePeAxis",kePeAxis);
    dimensionNames = ["kePeAxis", "t"];
    pi_g_ggw_kepe = diagfile.addVariable("pi_g_ggw_kepe",dimensionNames,type="double",isComplex=false);
else
    pi_w_wwg_kp = diagfile.variableWithName("pi_w_wwg_kp");
    F_wwg_kp = diagfile.variableWithName("F_wwg_kp");
    pi_g_ggw_kp = diagfile.variableWithName("pi_g_ggw_kp");
    F_ggw_kp = diagfile.variableWithName("F_ggw_kp");
    pi_w_wwg_omega = diagfile.variableWithName("pi_w_wwg_omega");
    pi_g_ggw_kepe = diagfile.variableWithName("pi_g_ggw_kepe");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Build the wavenumber/mode mask
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

valid = ~isnan(bins_0);
S_0 = sparse(find(valid), bins_0(valid), 1, numel(wvt.Ap), numel(kp), nnz(valid));

valid = ~isnan(bins_pm);
S_pm = sparse(find(valid), bins_pm(valid), 1, numel(wvt.Ap), numel(kp), nnz(valid));

mask_0 = false(wvt.Nj,wvt.Nkl,length(kp));
mask_pm = false(wvt.Nj,wvt.Nkl,length(kp));
for iK = 1:1:length(kp)
    mask_0(:,:,iK) = (bins_0 <= iK);
    mask_pm(:,:,iK) = (bins_pm <= iK);
end

mask_omega = false(wvt.Nj,wvt.Nkl,length(omegaAxis));
for iK = 1:1:length(omegaAxis)
    mask_omega(:,:,iK) = (bins_omega <= iK);
end

mask_kepe = false(wvt.Nj,wvt.Nkl,length(kePeAxis));
for iK = 1:1:length(kePeAxis)
    mask_kepe(:,:,iK) = (bins_kepe <= iK);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Loop over the the requested time indices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

integrationLastInformWallTime = datetime('now');
loopStartTime = integrationLastInformWallTime;
integrationLastInformLoopNumber = 1;
integrationInformTime = 10;
fprintf("Starting loop to compute the 1d mirror fluxes for %d time indices, over Nk=%d, Nomega=%d, Nkepe=%d.\n",length(timeIndices),length(kp),length(omegaAxis),length(kePeAxis));
for timeIndex = 1:length(timeIndices)
    deltaWallTime = datetime('now')-integrationLastInformWallTime;
    if ( seconds(deltaWallTime) > integrationInformTime)
        wallTimePerLoopTime = deltaWallTime / (timeIndex - integrationLastInformLoopNumber);
        wallTimeRemaining = wallTimePerLoopTime*(length(timeIndices) - timeIndex);
        fprintf('Time index %d of %d. Estimated time to finish is %s (%s)\n', timeIndex, length(timeIndices), wallTimeRemaining, datetime(datetime('now')+wallTimeRemaining,TimeZone='local',Format='d-MMM-y HH:mm:ss Z')) ;
        integrationLastInformWallTime = datetime('now');
        integrationLastInformLoopNumber = timeIndex;
    end

    outputIndex = timeIndices(timeIndex);
    self.iTime = timeIndices(timeIndex);

    E0 = WVDiagnostics.waveWaveGeostrophicEnergy(wvt,1);
    F_wwg_kp_val = reshape(E0(:).' * S_0,[],1);
    F_wwg_kp.setValueAlongDimensionAtIndex(F_wwg_kp_val,'t',outputIndex);

    Epm = WVDiagnostics.geostrophicGeostrophicWaveEnergy(wvt,1);
    F_ggw_kp_val = reshape(Epm(:).' * S_pm,[],1);
    F_ggw_kp.setValueAlongDimensionAtIndex(F_ggw_kp_val,'t',outputIndex);

    pi_w_wwg_kp_val = zeros(length(kp),1);
    pi_g_ggw_kp_val = zeros(length(kp),1);

    for i=1:length(kp)
        E0 = WVDiagnostics.waveWaveGeostrophicEnergy(wvt,mask_pm(:,:,i));
        pi_w_wwg_kp_val(i) = sum(E0(:));
        Epm = WVDiagnostics.geostrophicGeostrophicWaveEnergy(wvt,mask_0(:,:,i));
        pi_g_ggw_kp_val(i) = sum(Epm(:));
    end
    pi_w_wwg_kp.setValueAlongDimensionAtIndex(pi_w_wwg_kp_val,'t',outputIndex);
    pi_g_ggw_kp.setValueAlongDimensionAtIndex(pi_g_ggw_kp_val,'t',outputIndex);

    pi_w_wwg_omega_val = zeros(length(omegaAxis),1);
    for i=1:length(omegaAxis)
        E0 = WVDiagnostics.waveWaveGeostrophicEnergy(wvt,mask_omega(:,:,i));
        pi_w_wwg_omega_val(i) = sum(E0(:));
    end
    pi_w_wwg_omega.setValueAlongDimensionAtIndex(pi_w_wwg_omega_val,'t',outputIndex);

    pi_g_ggw_kepe_val = zeros(length(kePeAxis),1);
    for i=1:length(kePeAxis)
        Epm = WVDiagnostics.geostrophicGeostrophicWaveEnergy(wvt,mask_kepe(:,:,i));
        pi_g_ggw_kepe_val(i) = sum(Epm(:));
    end
    pi_g_ggw_kepe.setValueAlongDimensionAtIndex(pi_g_ggw_kepe_val,'t',outputIndex);
end
deltaLoopTime = datetime('now')-loopStartTime;
fprintf("Total loop time %s, which is %s per time index.\n",deltaLoopTime,deltaLoopTime/length(timeIndices));

end