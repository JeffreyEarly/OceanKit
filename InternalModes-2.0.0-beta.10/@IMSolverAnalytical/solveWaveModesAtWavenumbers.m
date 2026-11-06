function [collection,diagnostics] = solveWaveModesAtWavenumbers(self,background,kappa,options)
% Solve each distinct exact analytical wave page once.
%
% nModes counts all retained columns, including on kappa=0 pages. Supply a
% scalar or one count per requested page. At each distinct kappa, solve the
% largest count once and retain exact prefixes for smaller requests. Identical
% (kappa,count) pairs share one stored basis; no approximate sharing is used.
%
% - Topic: Solve analytical problems
% - Parameter background: authoritative source-specific physical column
% - Parameter kappa: nonempty requested nonnegative horizontal wavenumbers
% - Parameter options.nModes: required scalar or exact total counts aligned with every kappa entry
% - Parameter options.normalization: optional active profile normalization
% - Returns collection: capable analytical bases and requested-page mapping
% - Returns diagnostics: exact repeated-page reuse and retained counts
arguments
    self (1,1) IMSolverAnalytical
    background (1,1) IMBackground
    kappa (1,:) double {mustBeReal,mustBeFinite,mustBeNonnegative}
    options.nModes (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = []
    options.normalization = []
end
plan = IMModeSelection.wavePagePlan(kappa,options.nModes);
bases = cell(size(plan.solvedKappa));
for j = 1:numel(plan.solvedKappa)
    problem = IMInternalModes.waveModesAtWavenumber(background,k=plan.solvedKappa(j));
    bases{j} = self.solveEVP(problem,nModes=plan.solveCounts(j),normalization=options.normalization);
end
collection = IMModeSelection.wavePrefixCollection(bases,plan,self);
diagnostics = struct(requestedKappa=kappa,solvedKappa=plan.solvedKappa,basisIndex=plan.basisIndex,requestedPageCount=numel(kappa),solveCount=numel(plan.solvedKappa),reusedPageCount=numel(kappa)-numel(plan.solvedKappa),reuseKind="exactEqualKappa",solverClass="IMSolverAnalytical",nModes=options.nModes,requestedCounts=plan.pageCounts);
diagnostics.solves = table(plan.solvedKappa(:),plan.solveCounts(:),VariableNames=["kappa","requestedColumnCount"]);
end
