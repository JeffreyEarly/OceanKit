% % basedir = "/Users/Shared/CimRuns_June2025/output/";
% % basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";
% basedir = '/Volumes/SanDiskExtremePro/research/Energy-Pathways-Group/garrett-munk-spin-up/CimRuns_June2025_v2/output/';
% 
% 
% runNumber=18;
% wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");
% % wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","256") + ".nc");

wvd = WVDiagnostics("/Users/jearly/Documents/ProjectRepositories/EddyTideInteraction/bottom-generated-tide-forced-const-N-5cms.nc");

%%

if wvd.diagnosticsHasExplicitAntialiasing
    wvt = wvd.wvt_aa;
else
    wvt = wvd.wvt;
end

if isa(wvt,"WVTransformHydrostatic")
    h_pm = repmat(wvt.h_pm,[1 wvt.Nkl]);
else
    h_pm = wvt.h_pm;
end

dimensionNames = ["j","kRadial"];
[n,hN] = wvt.transformToRadialWavenumber(ones(size(wvt.Omega)),h_pm);
h_kj = (hN./n);
Lr2_pm = wvt.g * h_kj / (wvt.f*wvt.f);

var = wvd.diagfile.variableWithName("Lr2_pm");
var.value = Lr2_pm;

% wvd.diagfile.addVariable("Lr2_pm",dimensionNames,Lr2_pm,isComplex=false);