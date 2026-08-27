function support = wvCompiledBackendSupport(overrides)
release = string(version("-release"));
if ~startsWith(release,"R"), release = "R"+release; end
architecture = string(computer("arch"));
operatingSystem = string(systemIdentifier());
threadLimit = min(18,maxNumCompThreads);
if isfield(overrides,"Release"), release = string(overrides.Release); end
if ~startsWith(release,"R"), release = "R"+release; end
if isfield(overrides,"Architecture"), architecture = string(overrides.Architecture); end
if isfield(overrides,"OperatingSystem"), operatingSystem = string(overrides.OperatingSystem); end
if isfield(overrides,"MaxThreads"), threadLimit = min(18,double(overrides.MaxThreads)); end
releaseSupported = releaseAtLeast(release,"R2025b");
platformSupported = operatingSystem == "macOS" && architecture == "maca64";
reasons = strings(0,1);
if ~releaseSupported, reasons(end+1,1) = "MATLAB "+release+" is older than R2025b."; end
if operatingSystem ~= "macOS", reasons(end+1,1) = "The native provider requires macOS."; end
if architecture ~= "maca64", reasons(end+1,1) = "The native provider requires Apple-silicon maca64."; end
support = struct( ...
    "matlab",struct("release",release,"version",string(version),"minimumRelease","R2025b","isSupported",releaseSupported), ...
    "platform",struct("architecture",architecture,"operatingSystem",operatingSystem,"isSupported",platformSupported), ...
    "isSupported",releaseSupported && platformSupported, ...
    "reasons",reasons, ...
    "threadCount",threadLimit);
end

function identifier = systemIdentifier
if ismac
    identifier = "macOS";
elseif isunix
    identifier = "Linux";
elseif ispc
    identifier = "Windows";
else
    identifier = "unknown";
end
end

function tf = releaseAtLeast(actual,minimum)
actualTokens = regexp(actual,"^R(\d{4})([ab])$","tokens","once");
minimumTokens = regexp(minimum,"^R(\d{4})([ab])$","tokens","once");
if isempty(actualTokens) || isempty(minimumTokens)
    tf = false;
    return
end
actualValue = 2*str2double(actualTokens{1}) + double(actualTokens{2} == 'b');
minimumValue = 2*str2double(minimumTokens{1}) + double(minimumTokens{2} == 'b');
tf = actualValue >= minimumValue;
end
