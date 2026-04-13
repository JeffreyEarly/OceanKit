function fit = griddedStreamfunctionFromGroup(group)
arguments
    group NetCDFGroup
end

fit = GriddedStreamfunction.annotatedClassFromGroup(group);
end
