function jointLikelihood = EstimateSolutionLikelihoodFromBootstraps(bootstraps,parameters,stride)
% EstimateSolutionLikelihoodFromBootstraps Computes the log-likelihood of
% each bootstrap solution, by assuming the bootstraps themselves create a
% PDF.
%
% Inputs are,
%   bootstraps - struct of [u0,v0,u1,v1,sigma_n,sigma_s,zeta,delta] where
%                each variable is size [nT nBootstraps].
%   parameters - array of ModelParameter objects indicating which
%                parameters should be used to compute the likelihood.
%   stride - integer indicating how many time points to skip
%
% Output is,
%   jointLikelihood - log10(likelihood) of each bootstrap
%                     size is [1 nBootstraps]
%
% This is the implementation of equation 36 in Oscroft, Sykulski and Early.
%
% Gaussian bootstrap densities are fit with `KernelDensityEstimate`.

shouldEstimateU0V0 = 0;
shouldEstimateU1V1 = 0;
shouldEstimateStrain = 0;
shouldEstimateVorticity = 0;
shouldEstimateDivergence = 0;

nT = 0;
for i=1:length(parameters)
    if  parameters(i) == ModelParameter.u0v0
        shouldEstimateU0V0 = 1;
        if nT > 0 && size(bootstraps.u0,1) ~= nT
            error('inconsistent time series sizes');
        else
            nT = size(bootstraps.u0,1);
            nBootstraps = size(bootstraps.u0,2);
        end
    elseif  parameters(i) == ModelParameter.u1v1
        shouldEstimateU1V1 = 1;
        if nT > 0 && size(bootstraps.u1,1) ~= nT
            error('inconsistent time series sizes');
        else
            nT = size(bootstraps.u1,1);
            nBootstraps = size(bootstraps.u1,2);
        end
    elseif  parameters(i) == ModelParameter.strain
        shouldEstimateStrain = 1;
        if nT > 0 && size(bootstraps.sigma_n,1) ~= nT
            error('inconsistent time series sizes');
        else
            nT = size(bootstraps.sigma_n,1);
            nBootstraps = size(bootstraps.sigma_n,2);
        end
    elseif  parameters(i) == ModelParameter.vorticity
        shouldEstimateVorticity = 1;
        if nT > 0 && size(bootstraps.zeta,1) ~= nT
            error('inconsistent time series sizes');
        else
            nT = size(bootstraps.zeta,1);
            nBootstraps = size(bootstraps.zeta,2);
        end
    elseif  parameters(i) == ModelParameter.divergence
        shouldEstimateDivergence = 1;
        if nT > 0 && size(bootstraps.delta,1) ~= nT
            error('inconsistent time series sizes');
        else
            nT = size(bootstraps.delta,1);
            nBootstraps = size(bootstraps.delta,2);
        end
    end
end

m = bootstraps;
jointLikelihood = zeros(length(1:stride:nT),nBootstraps);
iOutput = 1;
for iTime=1:stride:nT
    if shouldEstimateU0V0 == 1    
        model = KernelDensityEstimate.fromData([m.u0(iTime, :).', m.v0(iTime, :).']);
        pointDensity = model.densityAt([m.u0(iTime, :).', m.v0(iTime, :).']);
        jointLikelihood(iOutput,:) = jointLikelihood(iOutput,:) + log10(safeDensity(pointDensity)).';
    end
    if shouldEstimateU1V1 == 1
        model = KernelDensityEstimate.fromData([m.u1(iTime, :).', m.v1(iTime, :).']);
        pointDensity = model.densityAt([m.u1(iTime, :).', m.v1(iTime, :).']);
        jointLikelihood(iOutput,:) = jointLikelihood(iOutput,:) + log10(safeDensity(pointDensity)).';
    end
    if shouldEstimateStrain == 1
        model = KernelDensityEstimate.fromData([m.sigma_n(iTime, :).', m.sigma_s(iTime, :).']);
        pointDensity = model.densityAt([m.sigma_n(iTime, :).', m.sigma_s(iTime, :).']);
        jointLikelihood(iOutput,:) = jointLikelihood(iOutput,:) + log10(safeDensity(pointDensity)).';
    end
    if shouldEstimateVorticity == 1
        model = KernelDensityEstimate.fromData(m.zeta(iTime, :).');
        pointDensity = model.densityAt(m.zeta(iTime, :).');
        jointLikelihood(iOutput,:) = jointLikelihood(iOutput,:) + log10(safeDensity(pointDensity)).';
    end
    if shouldEstimateDivergence == 1
        model = KernelDensityEstimate.fromData(m.delta(iTime, :).');
        pointDensity = model.densityAt(m.delta(iTime, :).');
        jointLikelihood(iOutput,:) = jointLikelihood(iOutput,:) + log10(safeDensity(pointDensity)).';
    end
    iOutput = iOutput + 1;
end
jointLikelihood = sum(jointLikelihood,1);
end

function density = safeDensity(values)
density = reshape(values, [], 1);
density(~isfinite(density) | density <= 0) = realmin("double");
end
