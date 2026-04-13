function fit = griddedStreamfunctionFromFile(path)
arguments
    path {mustBeTextScalar}
end

path = char(path);
ncfile = NetCDFFile(path, shouldReadOnly=true);
cleanup = onCleanup(@() ncfile.close()); %#ok<NASGU>
if isKey(ncfile.attributes, 'AnnotatedClass')
    className = string(ncfile.attributes('AnnotatedClass'));
    if ncfile.hasGroupWithName(className)
        group = ncfile.groupWithName(className);
    else
        group = ncfile;
    end
else
    error('GriddedStreamfunction:MissingAnnotatedClass', 'Unable to find the AnnotatedClass attribute in %s.', path);
end

fit = griddedStreamfunctionFromGroup(group);
end
