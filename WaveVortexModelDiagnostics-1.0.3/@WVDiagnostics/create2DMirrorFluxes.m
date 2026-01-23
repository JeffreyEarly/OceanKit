function create2DMirrorFluxes(self,options)
% Compute 2D mirror-flux diagnostics and write them to the diagnostics NetCDF.
%
% Compute 2D mirror-flux diagnostics and write them to the diagnostics NetCDF.
% Generates 2D mirror-flux summaries (triad primary and mirror variables)
% binned on the js/ks grid for the specified triad type ("wwg" or "ggw").
% Creates a diagnostics group (mirror-flux-2d-<triad>) and variables if they
% do not exist, then computes and appends the per-time-index 2D flux fields.
%
% - Topic: Diagnostics Generation
% - Declaration: create2DMirrorFluxes(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter stride: (optional) Stride between model time steps when creating new output (default: 1)
% - Parameter timeIndices: Explicit list of model time indices to process (default: uses stride over all times)
% - Parameter mirrorTriad: (optional) Which mirror triad to compute, one of "wwg" or "ggw" (default: "wwg")
% - Parameter shouldOverwriteExisting: (optional) (optional, logical) If true, overwrite any existing group data rather than appending (default: false)
arguments
    self WVDiagnostics
    options.stride = 1
    options.timeIndices
    options.mirrorTriad {mustBeMember(options.mirrorTriad, ["wwg","ggw"])} = "wwg";
    options.shouldOverwriteExisting = false;
end

if ~exist(self.diagpath,"file")
    error("No existing diagnostics file found.");
end
diagfile = self.diagfile;
wvt = self.wvt;

ks = self.sparseKRadialAxis;
js = self.sparseJWavenumberAxis;
[S_0, S_pm, mask_0, mask_pm] = self.sparseJKAxisBinMatrices();
N = length(js)*length(ks);
matrixSize = [length(js) length(ks)];

groupName = "mirror-flux-2d-" + options.mirrorTriad;
tName = "t_" + options.mirrorTriad;
if options.mirrorTriad=="wwg"
    triadPrimaryName = "F_wwg_js_ks";
    triadMirrorName = "pi_w_wwg_js_ks";
else
    triadPrimaryName = "F_ggw_js_ks";
    triadMirrorName = "pi_g_ggw_js_ks";
end

if diagfile.hasGroupWithName(groupName)
    %%%%%%%%%%%%%%%%%%
    % If existing file
    %%%%%%%%%%%%%%%%%%
    group = diagfile.groupWithName(groupName);
    if options.shouldOverwriteExisting
        outputIndexOffset = 0;
        if ~isfield(options,"timeIndices")
            timeIndices = 1:options.stride:length(self.t_wv);
        else
            timeIndices = options.timeIndices;
        end
    else
        t = group.readVariables(tName);
        [found, idx] = ismember(t, self.t_wv);
        if ~all(found)
            error('Some entries of mirror-fluxes-2d/t are not exactly in t_wv.');
        end
        if length(t) > 1
            stride = idx(2)-idx(1);
        else
            stride = options.stride;
        end
        timeIndices = (idx(end)+stride):stride:length(self.t_wv);
        if isfield(options,"timeIndices")
            timeIndices = intersect(options.timeIndices,timeIndices);
        end
        outputIndexOffset = length(t);
    end

    Pi_var = group.variableWithName(triadMirrorName);
    F_var = group.variableWithName(triadPrimaryName);
else
    %%%%%%%%%%%%%%%%%%
    % If new file
    %%%%%%%%%%%%%%%%%%
    outputIndexOffset = 0;

    if ~isfield(options,"timeIndices")
        timeIndices = 1:options.stride:length(self.t_wv);
    else
        timeIndices = options.timeIndices;
    end

    group = diagfile.addGroup(groupName);
    group.addDimension("js",js);
    group.addDimension("ks",ks);

    varAnnotation = wvt.propertyAnnotationWithName("t");
    varAnnotation.attributes('units') = varAnnotation.units;
    varAnnotation.attributes('long_name') = varAnnotation.description;
    varAnnotation.attributes('standard_name') = 'time';
    varAnnotation.attributes('long_name') = 'time';
    varAnnotation.attributes('units') = 'seconds since 1970-01-01 00:00:00';
    varAnnotation.attributes('axis') = 'T';
    varAnnotation.attributes('calendar') = 'standard';
    group.addDimension(tName,length=Inf,type="double",attributes=varAnnotation.attributes);

    dimensionNames = ["js", "ks", tName];
    Pi_var = group.addVariable(triadMirrorName,dimensionNames,type="double",isComplex=false);
    F_var = group.addVariable(triadPrimaryName,dimensionNames,type="double",isComplex=false);
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
fprintf("Starting loop to compute the 2d mirror fluxes for %d time indices, over N=%d.\n",length(timeIndices),N);
for timeIndex = 1:length(timeIndices)
    deltaWallTime = datetime('now')-integrationLastInformWallTime;
    if ( seconds(deltaWallTime) > integrationInformTime)
        wallTimePerLoopTime = deltaWallTime / (timeIndex - integrationLastInformLoopNumber);
        wallTimeRemaining = wallTimePerLoopTime*(length(timeIndices) - timeIndex + 1);
        fprintf('Time index %d of %d. Estimated time to finish is %s (%s)\n', timeIndex, length(timeIndices), wallTimeRemaining, datetime(datetime('now')+wallTimeRemaining,TimeZone='local',Format='d-MMM-y HH:mm:ss Z')) ;
        integrationLastInformWallTime = datetime('now');
        integrationLastInformLoopNumber = timeIndex;
    end

    outputIndex = timeIndex + outputIndexOffset;
    self.iTime = timeIndices(timeIndex);
    group.variableWithName(tName).setValueAlongDimensionAtIndex(wvt.t,tName,outputIndex);

    Pi_val = zeros(matrixSize);
    if options.mirrorTriad=="wwg"
        E0 = WVDiagnostics.waveWaveGeostrophicEnergy(wvt,1);
        F_val = reshape(reshape(E0(:).' * S_0,[],1),matrixSize);
        F_var.setValueAlongDimensionAtIndex(F_val,tName,outputIndex);

        for i=1:N
            E0 = WVDiagnostics.waveWaveGeostrophicEnergy(wvt,mask_pm(:,:,i));
            Pi_val(i) = sum(E0(:));
        end
    else
        Epm = WVDiagnostics.geostrophicGeostrophicWaveEnergy(wvt,1);
        F_val = reshape(reshape(Epm(:).' * S_pm,[],1),matrixSize);
        F_var.setValueAlongDimensionAtIndex(F_val,tName,outputIndex);

        for i=1:N
            Epm = WVDiagnostics.geostrophicGeostrophicWaveEnergy(wvt,mask_0(:,:,i));
            Pi_val(i) = sum(Epm(:));
        end
    end
    Pi_var.setValueAlongDimensionAtIndex(Pi_val,tName,outputIndex);

end
deltaLoopTime = datetime('now')-loopStartTime;
fprintf("Total loop time %s, which is %s per time index.\n",deltaLoopTime,deltaLoopTime/length(timeIndices));

end