repoRoot = fileparts(fileparts(mfilename('fullpath')));
cleanupRepoPath = addPathIfNeeded(repoRoot); %#ok<NASGU>

plainProperty = CAObjectProperty('plainValue', 'Plain object property.');
assert(strlength(string(plainProperty.className)) == 0, ...
    'Legacy CAObjectProperty construction should leave className empty.');
assert(strlength(string(plainProperty.sizeText)) == 0, ...
    'Legacy CAObjectProperty construction should leave sizeText empty.');

typedProperty = CAObjectProperty('typedValue', 'Typed object property.', ...
    className='AnnotationTestClassB', sizeText='nonempty vector');
assert(strcmp(string(typedProperty.className), "AnnotationTestClassB"), ...
    'CAObjectProperty should store documentation-only class metadata.');
assert(strcmp(string(typedProperty.sizeText), "nonempty vector"), ...
    'CAObjectProperty should store documentation-only size metadata.');

annotations = AnnotationTestClass.classDefinedPropertyAnnotations();
names = string({annotations.name});
myObjsAnnotation = annotations(names == "myObjs");
assert(isscalar(myObjsAnnotation), 'Expected a single myObjs annotation.');
assert(strcmp(string(myObjsAnnotation.className), "AnnotationTestClassB"), ...
    'Class-defined object annotations should retain className metadata.');
assert(strcmp(string(myObjsAnnotation.sizeText), "vector or empty array"), ...
    'Class-defined object annotations should retain sizeText metadata.');

function cleanup = addPathIfNeeded(pathToAdd)
if contains(path, [pathsep pathToAdd pathsep]) || startsWith(path, [pathToAdd pathsep]) || endsWith(path, [pathsep pathToAdd]) || strcmp(path, pathToAdd)
    cleanup = onCleanup(@() []);
else
    addpath(pathToAdd, '-begin');
    cleanup = onCleanup(@() rmpath(pathToAdd));
end
end
