function bootstrap = griddedStreamfunctionBootstrapFromGroup(group)
arguments
    group NetCDFGroup
end

bootstrap = GriddedStreamfunctionBootstrap.annotatedClassFromGroup(group);
end
