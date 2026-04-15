basedir = "/Users/Shared/CimRuns_June2025/output/";
basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";

%%
runNumber=1; runName = "hydrostatic: geostrophic";
wvd1 = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

%%
runNumber=9; runName = "hydrostatic: geostrophic + waves";
wvd9 = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

%%
wvd.plotFluidStateMultipanel();

%%
wvd.plotFluidDecompositionMultipanel();

%%
wvd.plotEnstrophySpectrum();

%%
wvd.plotEnergySpectrum();

%%
wvd.plotMooringRotarySpectrum();

%%
wvd.plotEnergyFluxTemporalAverage(approximation="exact",overSaturationFactor=100);

%%
wvd9.plotEnergyFluxTemporalAverage(approximation="quadratic",energyReservoir = EnergyReservoir.geostrophic_mda,overSaturationFactor=20,timeIndices=51:251);

%%
wvd.plotEnstrophyOverTime();

%%
wvd.plotEnergyOverTime();

%%
wvd.plotEnergyFluxOverTime(approximation="exact",filter=@(v,t) movmean(v,21));

%%
wvd.plotEnergyFluxOverTime(approximation="quadratic",filter=@(v,t) movmean(v,21));

%%
wvd1.plotEnergyTriadFluxOverTime(filter=@(v,t) movmean(v,51),energyReservoirs = [EnergyReservoir.geostrophic_kinetic, EnergyReservoir.geostrophic_potential_mda]);

%%
wvd.plotEnstrophyFluxOverTime(approximation="exact",filter=@(v,t) movmean(v,21));
wvd.plotEnstrophyFluxOverTime(approximation="quadratic",filter=@(v,t) movmean(v,21));

%%
wvd.plotEnstrophyTriadFluxOverTime(filter=@(v,t) movmean(v,51))