%% checkout_and_install_repos.m
% Clone/update all repos into the current working directory, then run
% mpminstall(..., Authoring=true, Force=true, InstallDependencies=false)

clear; clc;

%% OceanKit repo (checkout only)
oceanKit = struct( ...
    "name", "OceanKit", ...
    "url",  "https://github.com/JeffreyEarly/OceanKit.git");

%% Package repos
repos = [ ...
    struct("name","chebfun",                       "url","https://github.com/JeffreyEarly/chebfun.git")
    struct("name","class-annotations",             "url","https://github.com/JeffreyEarly/class-annotations.git")
    struct("name","class-docs",                    "url","https://github.com/JeffreyEarly/class-docs.git")
    struct("name","distributions",                 "url","https://github.com/JeffreyEarly/distributions.git")
    struct("name","geographic-projection",         "url","https://github.com/JeffreyEarly/geographic-projection.git")
    struct("name","internal-modes",                "url","https://github.com/JeffreyEarly/internal-modes.git")
    struct("name","netcdf",                        "url","https://github.com/JeffreyEarly/netcdf.git")
    struct("name","spline-core",                   "url","https://github.com/JeffreyEarly/spline-core.git")
    struct("name","wave-vortex-model",             "url","https://github.com/JeffreyEarly/wave-vortex-model.git")
    struct("name","wave-vortex-model-diagnostics", "url","https://github.com/Energy-Pathways-Group/wave-vortex-model-diagnostics.git")
    struct("name","AlongTrackSimulator",           "url","https://github.com/satmapkit/AlongTrackSimulator.git")
];

%% Use current working directory
baseDir = pwd;
fprintf('Using current working directory:\n  %s\n\n', baseDir);

%% Verify git exists
[status, out] = system("git --version");
if status ~= 0
    error("Git not found on system PATH.\nOutput:\n%s", out);
end
fprintf("Found %s\n\n", strtrim(out));

%% Checkout OceanKit (no installation)
oceanKitPath = fullfile(baseDir, oceanKit.name);
fprintf('============================================================\n');
fprintf('Checking out %s\n', oceanKit.name);
fprintf('============================================================\n');
cloneOrUpdateRepo(oceanKit.url, oceanKitPath);

%% Clone/update and install each repo locally
for k = 1:numel(repos)
    repoPath = fullfile(baseDir, repos(k).name);

    fprintf('\n============================================================\n');
    fprintf('Processing %s\n', repos(k).name);
    fprintf('============================================================\n');

    cloneOrUpdateRepo(repos(k).url, repoPath);

    fprintf('Installing with mpminstall (Authoring=true, Force=true, InstallDependencies=false):\n  %s\n', repoPath);
    try
        mpminstall(repoPath, Authoring=true, Prompt=false, InstallDependencies=false);
        fprintf('Installed %s successfully.\n', repos(k).name);
    catch ME
        warning("mpminstall failed for %s:\n%s", repos(k).name, ME.message);
    end
end

fprintf('\nDone.\n');

%% Helper function
function cloneOrUpdateRepo(repoUrl, repoPath)

    if isfolder(repoPath)
        fprintf('Repo already exists, updating:\n  %s\n', repoPath);

        cmd = sprintf('git -C "%s" pull --ff-only', repoPath);
        [status, out] = system(cmd);

        if status ~= 0
            warning("git pull failed for %s\nOutput:\n%s", repoPath, out);
        else
            fprintf('%s\n', strtrim(out));
        end

    else
        fprintf('Cloning:\n  %s\ninto:\n  %s\n', repoUrl, repoPath);

        cmd = sprintf('git clone "%s" "%s"', repoUrl, repoPath);
        [status, out] = system(cmd);

        if status ~= 0
            error("git clone failed for %s\nOutput:\n%s", repoUrl, out);
        else
            fprintf('%s\n', strtrim(out));
        end
    end

end