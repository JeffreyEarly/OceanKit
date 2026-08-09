repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
buildtool("test:full","-buildFile",fullfile(repositoryRoot,"buildfile.m"));
