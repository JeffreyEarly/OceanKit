%% See how aligned F/G products alias between physical channels
% Quadratic products of paired hydrostatic variables do not all project
% back into the same variable. The aligned transform assesses the three
% physically consistent channels F_i F_j -> F, F_i G_j -> G, and
% G_i G_j -> F. Constant stratification on a DCT-I/DST-I grid makes both
% the channel routing and the first aliased product exact and visible.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Build the exact cosine/sine family and endpoint rule
D = 4000;
N0 = 5.2e-3;
zDomain = [-D 0];
N2 = @(z) N0*N0*ones(size(z));

evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=zDomain);
solver = IMSolverSpectral(nEVP=160);
basisSet = solver.solveEVP(evp,nModes=30);
basisSet.normalization = "geostrophic";

M = 12;
z = linspace(zDomain(1),zDomain(2),M+1).';
weights = [0.5;ones(M-1,1);0.5]*(D/M);
variables = ["F","G"];
quadraticTolerance = 1e-8;
[transform, assessment] = basisSet.discreteTransform(z=z,weights=weights,variables=variables,quadraticAliasingTolerance=quadraticTolerance);

%% Recover the exact two-thirds cutoff
% With M equal grid intervals, quadratic products remain unaliased through
%
%   K_max = floor((2 M - 1)/3).
%
% Physical labels begin at zero, so the accepted aligned family contains
% K_max+1 columns. The assessment computes this cutoff from independently
% integrated continuous product projections rather than from the formula.
Kmax = floor((2*M-1)/3);
expectedFamilyCount = Kmax+1;
diagnostics = assessment.prefixDiagnostics;
firstRejectedRow = expectedFamilyCount+1;

fprintf("\nDCT-I/DST-I intervals M: %d\n",M);
fprintf("Exact K_max: %d; expected family count: %d; assessed family count: %d.\n",Kmax,expectedFamilyCount,assessment.retainedModeCount);
fprintf("Last accepted quadratic error: %.3e\n",diagnostics.quadraticAliasingError(expectedFamilyCount));
fprintf("First rejected quadratic error: %.3e\n",diagnostics.quadraticAliasingError(firstRejectedRow));

cutoffRows = diagnostics(max(1,expectedFamilyCount-1):min(height(diagnostics),firstRejectedRow+1), ...
    ["modeCount" "lastModeNumber" "quadraticAliasingError" "quadraticLimitingChannel" ...
    "quadraticLimitingModeNumberI" "quadraticLimitingModeNumberJ" "quadraticAccepted"]);
disp(cutoffRows)

%% Follow one resolved product through every target channel
% Modes 2 and 3 generate sum and difference wavenumbers 5 and 1, both
% inside the accepted family. Each product is projected into its physical
% target variable and reconstructed on the same sample points.
iColumn = find(transform.modeNumber == 2,1);
jColumn = find(transform.modeNumber == 3,1);
F = transform.inverseMatrix(variable="F");
G = transform.inverseMatrix(variable="G");

FF = F(:,iColumn).*F(:,jColumn);
FG = F(:,iColumn).*G(:,jColumn);
GG = G(:,iColumn).*G(:,jColumn);
FFReconstruction = transform.transformBack(transform.transformForward(FF,variable="F"),variable="F");
FGReconstruction = transform.transformBack(transform.transformForward(FG,variable="G"),variable="G");
GGReconstruction = transform.transformBack(transform.transformForward(GG,variable="F"),variable="F");

channel = ["FF->F";"FG->G";"GG->F"];
relativeResidual = [norm(FF-FFReconstruction)/norm(FF);norm(FG-FGReconstruction)/norm(FG);norm(GG-GGReconstruction)/norm(GG)];
resolvedProducts = table(channel,relativeResidual);
disp(resolvedProducts)

figure(Name="V2 resolved coupled F-G products",Color="w");
tiledlayout(1,3,TileSpacing="compact",Padding="compact");

nexttile
plot(FF,z,"k-",FFReconstruction,z,"--",LineWidth=1.3)
grid on
xlabel("F_2 F_3")
ylabel("z (m)")
title("FF projects to F")

nexttile
plot(FG,z,"k-",FGReconstruction,z,"--",LineWidth=1.3)
grid on
xlabel("F_2 G_3")
ylabel("z (m)")
title("FG projects to G")

nexttile
plot(GG,z,"k-",GGReconstruction,z,"--",LineWidth=1.3)
grid on
xlabel("G_2 G_3")
ylabel("z (m)")
title("GG projects to F")
legend(["sampled product" "target-channel reconstruction"],Location="best")

%% Inspect the first product rejected by the policy
% The combined diagnostic identifies both source labels and the physical
% product channel. Rebuild that prefix without the quadratic policy so its
% sampled alias can be inspected directly rather than intentionally
% throwing the strict-policy exception.
limitingChannel = diagnostics.quadraticLimitingChannel(firstRejectedRow);
modeNumberI = diagnostics.quadraticLimitingModeNumberI(firstRejectedRow);
modeNumberJ = diagnostics.quadraticLimitingModeNumberJ(firstRejectedRow);
firstRejectedTransform = basisSet.discreteTransform(z=z,weights=weights,nModes=firstRejectedRow,variables=variables,gramTolerance=1);
iColumn = find(firstRejectedTransform.modeNumber == modeNumberI,1);
jColumn = find(firstRejectedTransform.modeNumber == modeNumberJ,1);
F = firstRejectedTransform.inverseMatrix(variable="F");
G = firstRejectedTransform.inverseMatrix(variable="G");

if limitingChannel == "FF->F"
    aliasedProduct = F(:,iColumn).*F(:,jColumn);
    targetVariable = "F";
elseif limitingChannel == "FG->G"
    aliasedProduct = F(:,iColumn).*G(:,jColumn);
    targetVariable = "G";
else
    aliasedProduct = G(:,iColumn).*G(:,jColumn);
    targetVariable = "F";
end
aliasedCoefficients = firstRejectedTransform.transformForward(aliasedProduct,variable=targetVariable);
[dominantAliasMagnitude,dominantAliasColumn] = max(abs(aliasedCoefficients));
dominantAliasModeNumber = firstRejectedTransform.modeNumber(dominantAliasColumn);

fprintf("First rejected product: %s from physical modes %g and %g.\n",limitingChannel,modeNumberI,modeNumberJ);
fprintf("Its sampled projection is dominated by mode %g with magnitude %.3e.\n",dominantAliasModeNumber,dominantAliasMagnitude);

figure(Name="V2 first coupled-product alias",Color="w");
stem(firstRejectedTransform.modeNumber,abs(aliasedCoefficients),"filled",LineWidth=1.2)
grid on
xlabel("target physical mode number")
ylabel("absolute sampled coefficient")
title(sprintf("First rejected %s product aliases into the retained band",limitingChannel))

%% Relate the exact example to general stratification
% The exact cutoff belongs to this uniform constant-stratification rule.
% For exponential stratification, `DiscreteTransformQualityPolicies.m`
% applies the same three product channels using refined continuous
% projections and reports the worst channel at every candidate prefix.
