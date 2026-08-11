repoRoot = fileparts(fileparts(mfilename('fullpath')));
cleanupRepoPath = addPathIfNeeded(repoRoot);

helperRoot = tempname;
mkdir(helperRoot);
cleanupHelperRoot = onCleanup(@() removeDirectoryIfPresent(helperRoot));

cases = {
    'm', 'units of $$\mathrm{m}$$.'
    'm s-1', 'units of $$\mathrm{m\,s^{-1}}$$.'
    'm/s', 'units of $$\mathrm{m\,s^{-1}}$$.'
    'kg m-3', 'units of $$\mathrm{kg\,m^{-3}}$$.'
    'kg/m3', 'units of $$\mathrm{kg\,m^{-3}}$$.'
    'rad m-1', 'units of $$\mathrm{rad\,m^{-1}}$$.'
    'rad2 m-2', 'units of $$\mathrm{rad^{2}\,m^{-2}}$$.'
    'kg m-1 s-2', 'units of $$\mathrm{kg\,m^{-1}\,s^{-2}}$$.'
    'kg/m/s2', 'units of $$\mathrm{kg\,m^{-1}\,s^{-2}}$$.'
    'm^2 s^{-1}', 'units of $$\mathrm{m^{2}\,s^{-1}}$$.'
    's^-1', 'units of $$\mathrm{s^{-1}}$$.'
    'rad m^{-1}', 'units of $$\mathrm{rad\,m^{-1}}$$.'
    'm^2/rad^2', 'units of $$\mathrm{m^{2}\,rad^{-2}}$$.'
    '1', 'is dimensionless.'
    'degrees_north', 'units of degrees north ($$^\circ\mathrm{N}$$).'
    'degrees_east', 'units of degrees east ($$^\circ\mathrm{E}$$).'
    '', 'no units.'
    'mode number', 'units of `mode number`.'
    'm//s', 'units of `m//s`.'
    'percent%unit', 'units of `percent%unit`.'
    };

for iCase = 1:size(cases,1)
    rawUnits = cases{iCase,1};
    expectedText = cases{iCase,2};
    documentation = MethodDocumentation("sample" + iCase);
    documentation.shortDescription = "Sample property.";
    documentation.functionType = FunctionType.transformProperty;
    documentation.dimensions = {'sample'};
    documentation.units = rawUnits;
    documentation.isComplex = false;
    documentation.detailedDescription = "Further context.";
    documentation.pathOfOutputFile = fullfile(helperRoot,"sample" + iCase + ".md");
    documentation.writeToFile("Example",iCase,"Classes");

    generatedText = string(fileread(documentation.pathOfOutputFile));
    assert(contains(generatedText,expectedText), ...
        'Unit value "%s" did not render as expected.',rawUnits);
    expectedParagraph = "Real valued property with dimension $$sample$$ and " + expectedText + newline + newline + "## Discussion" + newline + "Further context.";
    assert(contains(generatedText,expectedParagraph), ...
        'The complete property description must retain spacing and real Markdown line breaks.');
    assert(~contains(generatedText,"andunits") && ~contains(generatedText,"\\n\\n"), ...
        'Property descriptions must not contain collapsed words or literal newline escapes.');
    assert(strcmp(string(documentation.units),string(rawUnits)), ...
        'Rendering must not modify the raw units metadata.');
end

recognizedLegacyForms = ["m/s", "kg/m3", "kg/m/s2", "m^2 s^{-1}", "s^-1", "rad m^{-1}", "m^2/rad^2"];
for iCase = 1:numel(recognizedLegacyForms)
    caseIndex = find(strcmp(string(cases(:,1)),recognizedLegacyForms(iCase)),1);
    generatedText = string(fileread(fullfile(helperRoot,"sample" + caseIndex + ".md")));
    assert(~contains(generatedText,"units of $$" + recognizedLegacyForms(iCase) + "$$"), ...
        'Recognized legacy unit forms should not survive as raw MathJax.');
end

function cleanup = addPathIfNeeded(pathToAdd)
if contains(path, [pathsep pathToAdd pathsep]) || startsWith(path, [pathToAdd pathsep]) || endsWith(path, [pathsep pathToAdd]) || strcmp(path, pathToAdd)
    cleanup = onCleanup(@() []);
else
    addpath(pathToAdd, '-begin');
    cleanup = onCleanup(@() rmpath(pathToAdd));
end
end

function removeDirectoryIfPresent(pathToRemove)
if isfolder(pathToRemove)
    rmdir(pathToRemove, 's');
end
end
