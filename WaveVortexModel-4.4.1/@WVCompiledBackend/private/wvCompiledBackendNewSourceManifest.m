function manifest = wvCompiledBackendNewSourceManifest(packageRoot)
% Compute the deterministic source identity used to build the native bridge.
arguments
    packageRoot (1,1) string
end
packageRoot = canonicalPath(packageRoot);
treeRoots = [ ...
    "CompiledKernel/src" ...
    "CompiledKernel/include" ...
    "CompiledKernel/adapters/native-fftw" ...
    "CompiledKernel/adapters/accelerate" ...
    "PortableRuntime/src" ...
    "PortableRuntime/include" ...
    "@WVCompiledBackend" ...
    "@WVCompiledTransformBackend" ...
    "@WVCompiledConstantStratificationBackend"];
rootEntries = dir(packageRoot);
rootDirectoryNames = string({rootEntries([rootEntries.isdir]).name});
publicClassRoots = rootDirectoryNames(startsWith(rootDirectoryNames,"@WVTransform") | ...
    startsWith(rootDirectoryNames,"@WVGeometry"));
treeRoots = unique([treeRoots publicClassRoots "@WVModel" "Forcing" "ObservingSystems" "Operations" "FlowComponents"],"stable");
relativePaths = strings(0,1);
for treeRoot = treeRoots
    root = fullfile(packageRoot,treeRoot);
    if ~isfolder(root), continue, end
    entries = dir(fullfile(root,"**","*"));
    entries = entries(~[entries.isdir]);
    for index = 1:numel(entries)
        relative = extractAfter(string(fullfile(entries(index).folder,entries(index).name)),strlength(packageRoot)+1);
        if ~any(endsWith(lower(relative),[".cpp" ".hpp" ".h" ".m"])), continue, end
        relativePaths(end+1,1) = replace(relative,filesep,"/"); %#ok<AGROW>
    end
end
rootMatlabFiles = dir(fullfile(packageRoot,"*.m"));
relativePaths = [relativePaths;string({rootMatlabFiles.name}).'];
providerManifest = "CompiledKernel/native-fftw-provider.env";
if isfile(fullfile(packageRoot,providerManifest))
    relativePaths(end+1,1) = providerManifest;
end
relativePaths = sort(unique(relativePaths));
if isempty(relativePaths)
    error("WaveVortexModel:CompiledBackendSourceManifestEmpty","The compiled-backend source manifest is empty under %s.",packageRoot);
end
absolutePaths = arrayfun(@(value)fullfile(packageRoot,replace(value,"/",filesep)),relativePaths);
hashes = sha256Files(absolutePaths);
files = repmat(struct("path","","sha256","","bytes",0),numel(relativePaths),1);
for index = 1:numel(relativePaths)
    information = dir(absolutePaths(index));
    files(index) = struct("path",relativePaths(index),"sha256",hashes(index),"bytes",double(information.bytes));
end
manifest = struct( ...
    "schemaVersion","1.0.0", ...
    "algorithm","sha256(path-tab-content-sha256-newline-v1)", ...
    "aggregateSHA256",aggregateHash(relativePaths,hashes), ...
    "fileCount",numel(relativePaths), ...
    "files",files);
end

function hashes = sha256Files(paths)
command = "/usr/bin/shasum -a 256 --";
for pathname = reshape(paths,1,[])
    command = command+" "+shellQuote(pathname);
end
[status,output] = system(command);
if status ~= 0
    error("WaveVortexModel:CompiledBackendHash","Unable to hash the compiled-backend source manifest.");
end
lines = splitlines(strtrim(string(output)));
if numel(lines) ~= numel(paths)
    error("WaveVortexModel:CompiledBackendHash","The source-manifest hash output was incomplete.");
end
hashes = extractBefore(lines,65);
if any(strlength(hashes) ~= 64)
    error("WaveVortexModel:CompiledBackendHash","The source-manifest hash output was invalid.");
end
end

function hash = aggregateHash(paths,hashes)
temporary = string(tempname);
cleanup = onCleanup(@()deleteIfPresent(temporary));
fileId = fopen(temporary,"w");
if fileId < 0
    error("WaveVortexModel:CompiledBackendHash","Unable to stage the source-manifest aggregate.");
end
fileCleanup = onCleanup(@()fclose(fileId));
for index = 1:numel(paths)
    fprintf(fileId,"%s\t%s\n",paths(index),hashes(index));
end
clear fileCleanup
[status,output] = system("/usr/bin/shasum -a 256 "+shellQuote(temporary));
if status ~= 0
    error("WaveVortexModel:CompiledBackendHash","Unable to hash the source-manifest aggregate.");
end
hash = extractBefore(string(strtrim(output))," ");
clear cleanup
end

function pathname = canonicalPath(pathname)
[status,output] = system("/bin/realpath "+shellQuote(pathname));
if status ~= 0
    error("WaveVortexModel:CompiledBackendPath","Unable to resolve %s.",pathname);
end
pathname = string(strtrim(output));
end

function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end

function deleteIfPresent(pathname)
if isfile(pathname), delete(pathname); end
end
