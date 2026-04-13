function ensureBootstrapDiagnostics(self)
if numel(self.bootstrapKappaCache) == self.nBootstraps && numel(self.bootstrapCoherenceCache) == self.nBootstraps
    return
end

if self.nBootstraps == 0
    self.bootstrapKappaCache = zeros(1, 0);
    self.bootstrapCoherenceCache = zeros(1, 0);
    return
end

observedTotalDisplacementWeights = GriddedStreamfunctionBootstrap.totalDisplacementWeightsForSampleTimes( ...
    self.observedTrajectorySampleData.tCell);
self.bootstrapKappaCache = zeros(1, self.nBootstraps);
self.bootstrapCoherenceCache = zeros(1, self.nBootstraps);

for iBootstrap = 1:self.nBootstraps
    fit = self.fitForBootstrap(iBootstrap);
    diagnostics = GriddedStreamfunctionBootstrap.diagnosticsFromFit( ...
        fit, observedTotalDisplacementWeights(self.bootstrapIndices(iBootstrap, :)));
    self.bootstrapKappaCache(1, iBootstrap) = diagnostics.kappa;
    self.bootstrapCoherenceCache(1, iBootstrap) = diagnostics.coherence;
end
end
