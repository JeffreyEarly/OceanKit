profiles = cell(1,1);
profiles{1} = 'constant';
profiles{2} = 'exponential';
profiles{3} = 'pycnocline-constant';
profiles{4} = 'pycnocline-exponential';
profiles{5} = 'latmix-site1';
profiles{6} = 'latmix-site1-surface';
profiles{7} = 'latmix-site1-exponential';
profiles{8} = 'latmix-site1-constant';

methods = cell(5,1);
methods{1} = 'finiteDifference';
methods{2} = 'wkbSpectral';
methods{3} = 'densitySpectral';
methods{4} = 'spectral';
methods{5} = 'wkbAdaptiveSpectral';

for iProfile=1:length(profiles)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize the analytical solution
    n = 2*64;
    latitude = 33;
    [rhoFunction, N2Function, zIn] = InternalModes.StratificationProfileWithName(profiles{iProfile});
    z = linspace(min(zIn),max(zIn),n)';
    rho = rhoFunction(z);
    monotonicityTolerance = 1e-12*max(1, max(abs(rho)));
    isMonotonicProfile = all(diff(rho) <= monotonicityTolerance);
    if strcmp(profiles{iProfile}, 'latmix-site1-surface')
        isMonotonicProfile = false;
    end

    fprintf('\n---%s stratification profile---\n',profiles{iProfile});
    for iMethod=1:length(methods)
        methodRequiresMonotonicDensity = ...
            strcmp(methods{iMethod}, 'wkbSpectral') || ...
            strcmp(methods{iMethod}, 'densitySpectral') || ...
            strcmp(methods{iMethod}, 'wkbAdaptiveSpectral');
        if ~isMonotonicProfile && methodRequiresMonotonicDensity
            fprintf('Skipping %s because this profile is not monotonic.\n\n', methods{iMethod});
            continue;
        end

        % initialize directly from the function
        try
            im = InternalModes(rhoFunction,zIn,z,latitude,'nModes',n, 'method', methods{iMethod}, 'shouldShowDiagnostics', 1);
        catch ME
            fprintf('*******FAILED**********\n');
            if ~isempty(ME.cause)
                fprintf('%s : %s\n', ME.cause{1}.identifier, ME.cause{1}.message);
            else
                fprintf('%s : %s\n', ME.identifier, ME.message);
            end
            fprintf('***********************\n');
        end
        fprintf('\n')
        
        
        % initialize directly from a equispaced grid.
        try
            im = InternalModes(rho,z,z,latitude,'nModes',n, 'method', methods{iMethod}, 'shouldShowDiagnostics', 1);
        catch ME
            fprintf('*******FAILED**********\n');
            if ~isempty(ME.cause)
                fprintf('%s : %s\n', ME.cause{1}.identifier, ME.cause{1}.message);
            else
                fprintf('%s : %s\n', ME.identifier, ME.message);
            end
            fprintf('***********************\n');
        end
        fprintf('\n')
    end
end
