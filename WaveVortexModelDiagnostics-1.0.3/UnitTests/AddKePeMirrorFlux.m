% basedir = "/Users/Shared/CimRuns_June2025/output/";
% basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";
% basedir = "/Volumes/Samsung_T7/CimRuns_June2025/output/";
% basedir = '/Users/cwortham/Documents/research/Energy-Pathways-Group/garrett-munk-spin-up/CimRuns/output/';
basedir = '/Volumes/SanDiskExtremePro/research/Energy-Pathways-Group/garrett-munk-spin-up/CimRuns_June2025_v2/output/';

runs = 1:18;
% runs = [1,9,18]; % do these at 512 after.

for ii=1:length(runs)


    % runNumber=1;
    runNumber = runs(ii);
    wvd = WVDiagnostics(basedir + getRunParameters(runNumber) + ".nc");
    % wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

    fprintf('Making KE/PE mirror flux for runNumber %d\n',runNumber)

    %%
    
    diagfile = wvd.diagfile;
    wvt = wvd.wvt;
    
    t = diagfile.readVariables("t");
    timeIndices = 1:length(t);
    
    [kePeAxis,bins_kepe] = wvd.sparseKePeAxis;
    
    diagfile.addDimension("kePeAxis",kePeAxis);
    dimensionNames = ["kePeAxis", "t"];
    pi_g_ggw_kepe = diagfile.addVariable("pi_g_ggw_kepe",dimensionNames,type="double",isComplex=false);
    
    mask_kepe = false(wvt.Nj,wvt.Nkl,length(kePeAxis));
    for iK = 1:1:length(kePeAxis)
        mask_kepe(:,:,iK) = (bins_kepe <= iK);
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
        wvd.iTime = timeIndices(timeIndex);
    
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