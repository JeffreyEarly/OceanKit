function nc_truncate_unlimdim_in_group(infile, outfile, targetGroupPath, targetDimName, Ndrop)
% Copy a NetCDF-4 file (with groups) to a new file, truncating one dimension
% in one specific group by dropping Ndrop samples from the end.
%
% Example:
%   nc_truncate_unlimdim_in_group("in.nc","out.nc","/wave-vortex","t",1)

arguments
    infile (1,1) string
    outfile (1,1) string
    targetGroupPath (1,1) string
    targetDimName (1,1) string
    Ndrop (1,1) double {mustBeInteger,mustBeNonnegative} = 1
end

if isfile(outfile)
    delete(outfile);
end

rootInfo = ncinfo(infile);
ncid_out = netcdf.create(outfile, "NETCDF4");

% Inherited dimension maps:
%  - nameMap: dim name -> dimid (nearest ancestor wins)
%  - pathMap: absolute dim path (e.g. "/j" or "/wave-vortex/j") -> dimid
nameMap0 = containers.Map("KeyType","char","ValueType","double");
pathMap0 = containers.Map("KeyType","char","ValueType","double");

% PASS 1: define everything (metadata only)
node = defineGroupRecursive(ncid_out, rootInfo, "/", targetGroupPath, targetDimName, Ndrop, nameMap0, pathMap0);

% End define mode ONCE for the whole file
netcdf.endDef(ncid_out);

% PASS 2: write data
writeGroupRecursive(infile, node, targetGroupPath, targetDimName, Ndrop);

netcdf.close(ncid_out);
end

% ========================================================================
% PASS 1: DEFINE (with inherited dimension maps)
% ========================================================================
function node = defineGroupRecursive(outGid, gInfo, thisPath, targetGroupPath, targetDimName, Ndrop, inheritedNameMap, inheritedPathMap)

thisPathN   = normalizePath(thisPath);
targetPathN = normalizePath(targetGroupPath);
isTargetGroup = (thisPathN == targetPathN);

node = struct();
node.path = thisPathN;
node.gid  = outGid;

% Copy group attributes
copyAttributesToGroup(outGid, gInfo.Attributes);

% Make local copies of inherited maps so we can extend/override for this subtree
nameMap = copyMap(inheritedNameMap);
pathMap = copyMap(inheritedPathMap);

% Define dimensions in this group and extend maps
for d = 1:numel(gInfo.Dimensions)
    dim = gInfo.Dimensions(d);

    isUnlim = isfield(dim, "Unlimited") && logical(dim.Unlimited);

    if isTargetGroup && strcmp(dim.Name, targetDimName)
        oldLen = dim.Length;
        newLen = oldLen - Ndrop;
        if newLen < 0
            error("Truncation would make dimension '%s' negative (oldLen=%d, Ndrop=%d).", ...
                targetDimName, oldLen, Ndrop);
        end

        if isUnlim
            % Keep it UNLIMITED (mutable). We'll just write fewer records.
            dimid = netcdf.defDim(outGid, dim.Name, netcdf.getConstant("NC_UNLIMITED"));
        else
            % If it isn't unlimited, fall back to fixed (but your use-case said it is unlimited)
            dimid = netcdf.defDim(outGid, dim.Name, newLen);
        end
    else
        if isUnlim
            % Preserve unlimited-ness
            dimid = netcdf.defDim(outGid, dim.Name, netcdf.getConstant("NC_UNLIMITED"));
        else
            dimid = netcdf.defDim(outGid, dim.Name, dim.Length);
        end
    end

    dimid = double(dimid);

    % Update inherited maps
    nameMap(dim.Name) = dimid;
    absKey = makeDimAbsPath(thisPathN, dim.Name);
    pathMap(char(absKey)) = dimid;
end

node.nameMap = nameMap; %#ok<STRNU>
node.pathMap = pathMap; %#ok<STRNU>

% Define variables (store varids)
varIdMap = containers.Map("KeyType","char","ValueType","double");

for v = 1:numel(gInfo.Variables)
    var = gInfo.Variables(v);

    dimids = zeros(1, numel(var.Dimensions), "double");
    for k = 1:numel(var.Dimensions)
        ref = string(var.Dimensions(k).Name);

        dimids(k) = resolveDimId(ref, nameMap, pathMap);
    end

    xtype = xtypeFromNcinfoDatatype(var.Datatype);
    varid = netcdf.defVar(outGid, var.Name, xtype, dimids);
    varid = double(varid);
    varIdMap(var.Name) = varid;

    % Copy variable attributes
    for a = 1:numel(var.Attributes)
        att = var.Attributes(a);
        netcdf.putAtt(outGid, varid, att.Name, att.Value);
    end
end

node.varIdMap = varIdMap; %#ok<STRNU>
node.varInfo  = gInfo.Variables; %#ok<STRNU>

% Recurse into subgroups (still DEFINE mode)
node.sub = {};
for gg = 1:numel(gInfo.Groups)
    sub = gInfo.Groups(gg);

    subGid  = netcdf.defGrp(outGid, sub.Name);
    subPath = joinPath(thisPathN, sub.Name);

    node.sub{end+1,1} = defineGroupRecursive(subGid, sub, subPath, targetGroupPath, targetDimName, Ndrop, nameMap, pathMap);
end
end

% ========================================================================
% PASS 2: WRITE
% ========================================================================
function writeGroupRecursive(infile, node, targetGroupPath, targetDimName, Ndrop)

isTargetGroup = (normalizePath(node.path) == normalizePath(targetGroupPath));

for v = 1:numel(node.varInfo)
    var = node.varInfo(v);

    inVarPath = makeVarPath(node.path, var.Name);
    data = ncread(infile, inVarPath);

    if isTargetGroup
        dimRefs = getVarDimRefs(var);          % safe even for scalars
        baseNames = stripDimRefToBaseName(dimRefs);
        it = find(baseNames == targetDimName, 1);

        if ~isempty(it)
            oldLen = size(data, it);
            newLen = oldLen - Ndrop;
            if newLen < 0
                error("Variable '%s' would be truncated below 0 along '%s'.", inVarPath, targetDimName);
            end
            idx = repmat({':'}, 1, ndims(data));
            idx{it} = 1:newLen;
            data = data(idx{:});
        end
    end

    varid_out = node.varIdMap(var.Name);
    putVarWhole(node.gid, varid_out, data);
end

for k = 1:numel(node.sub)
    writeGroupRecursive(infile, node.sub{k}, targetGroupPath, targetDimName, Ndrop);
end
end

% ========================================================================
% HELPERS
% ========================================================================
function dimid = resolveDimId(ref, nameMap, pathMap)
% ref can be "j" or "/j" or "/some/group/j"

ref = string(ref);

if startsWith(ref, "/")
    absKey = normalizePath(ref);               % "/j" or "/some/group/j"
    if isKey(pathMap, char(absKey))
        dimid = pathMap(char(absKey));
        return;
    end

    % Fallback: sometimes ncinfo gives "/j" even though we stored "/j" already;
    % or it might give "/wave-vortex/j" but dim is actually root "/j".
    base = stripDimRefToBaseName(ref);
    absKey2 = "/" + base;
    if isKey(pathMap, char(absKey2))
        dimid = pathMap(char(absKey2));
        return;
    end

    error("Dimension reference '%s' not found in absolute dim map.", ref);
else
    % Plain name: allow inherited lookup by name
    if isKey(nameMap, char(ref))
        dimid = nameMap(char(ref));
        return;
    end

    error("Dimension '%s' not found in inherited name map.", ref);
end
end

function out = stripDimRefToBaseName(refs)
% vectorized: refs may be string array
refs = string(refs);
out = refs;
for i = 1:numel(refs)
    r = refs(i);
    if contains(r, "/")
        parts = split(r, "/");
        parts(parts=="") = [];
        out(i) = parts(end);
    else
        out(i) = r;
    end
end
end

function m2 = copyMap(m1)
m2 = containers.Map("KeyType", m1.KeyType, "ValueType", m1.ValueType);
ks = m1.keys;
for i = 1:numel(ks)
    k = ks{i};
    m2(k) = m1(k);
end
end

function copyAttributesToGroup(gid, attributes)
for a = 1:numel(attributes)
    att = attributes(a);
    netcdf.putAtt(gid, netcdf.getConstant("NC_GLOBAL"), att.Name, att.Value);
end
end

function xtype = xtypeFromNcinfoDatatype(dt)
dt = lower(string(dt));
switch dt
    case "double"
        xtype = netcdf.getConstant("NC_DOUBLE");
    case {"single","float"}
        xtype = netcdf.getConstant("NC_FLOAT");

    case "int64"
        xtype = netcdf.getConstant("NC_INT64");
    case "uint64"
        xtype = netcdf.getConstant("NC_UINT64");

    case {"int32","int"}
        xtype = netcdf.getConstant("NC_INT");
    case {"uint32","uint"}
        xtype = netcdf.getConstant("NC_UINT");

    case {"int16","short"}
        xtype = netcdf.getConstant("NC_SHORT");
    case {"uint16","ushort"}
        xtype = netcdf.getConstant("NC_USHORT");

    case {"int8","byte"}
        xtype = netcdf.getConstant("NC_BYTE");
    case {"uint8","ubyte"}
        xtype = netcdf.getConstant("NC_UBYTE");

    case "char"
        xtype = netcdf.getConstant("NC_CHAR");

    case "string"
        xtype = netcdf.getConstant("NC_STRING"); % NETCDF4 only

    otherwise
        error("Unsupported Datatype '%s' (compound/enum/opaque/vlen types need extra handling).", dt);
end
end

function p = normalizePath(p)
p = char(p);
if isempty(p), p = "/"; end
if p(1) ~= '/', p = ['/' p]; end %#ok<AGROW>
if numel(p) > 1 && p(end) == '/', p(end) = []; end
p = string(p);
end

function p = joinPath(parent, name)
parent = normalizePath(parent);
name = string(name);
if parent == "/"
    p = "/" + name;
else
    p = parent + "/" + name;
end
end

function absKey = makeDimAbsPath(groupPath, dimName)
groupPath = normalizePath(groupPath);
dimName = string(dimName);
if groupPath == "/"
    absKey = "/" + dimName;
else
    absKey = groupPath + "/" + dimName;
end
end

function varPath = makeVarPath(groupPath, varName)
groupPath = normalizePath(groupPath);
varName = string(varName);
if groupPath == "/"
    varPath = varName;
else
    varPath = groupPath + "/" + varName;
end
end

function dimRefs = getVarDimRefs(var)
% Return dimension references (names) for a variable as a string row vector.
% Robust to scalar vars where var.Dimensions may be [] (non-struct).

if ~isfield(var, "Dimensions") || isempty(var.Dimensions)
    dimRefs = strings(1,0);
    return;
end

% In some MATLAB/netcdf cases, scalar dims come in as [] not struct; guard:
if ~isstruct(var.Dimensions)
    dimRefs = strings(1,0);
    return;
end

if isempty(var.Dimensions)
    dimRefs = strings(1,0);
    return;
end

dimRefs = string({var.Dimensions.Name});
end

function putVarWhole(gid, varid, data)
% Write entire variable using explicit START/COUNT (required if any dim is unlimited).

% Query variable rank and shape from the output file
[~, ~, dimids, ~] = netcdf.inqVar(gid, varid);
nd = numel(dimids);

start = zeros(1, nd);              % netcdf uses 0-based indexing
count = zeros(1, nd);

for i = 1:nd
    [~, dimlen] = netcdf.inqDim(gid, dimids(i));
    if dimlen == netcdf.getConstant("NC_UNLIMITED")
        % In data mode, inqDim returns current length for unlimited dims in many builds,
        % but to be safe, use the size of the data we're writing.
        % (If dimlen comes back as 0, size(data,i) is what we want anyway.)
        dimlen = size(data, i);
    end
    if dimlen == 0
        dimlen = size(data, i);
    end
    count(i) = dimlen;
end

% For scalars nd=0: putVar with no start/count is OK
if nd == 0
    netcdf.putVar(gid, varid, data);
else
    % Ensure count matches actual data size (especially after truncation)
    sz = size(data);
    % MATLAB: trailing singleton dims may be omitted; pad sz to nd
    if numel(sz) < nd
        sz(end+1:nd) = 1;
    end
    count = double(sz(1:nd));
    start = double(start);

    netcdf.putVar(gid, varid, start, count, data);
end
end