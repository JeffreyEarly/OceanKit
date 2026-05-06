testPath = fullfile(fileparts(mfilename('fullpath')), 'test.nc');
cleanupTestFile = onCleanup(@() deleteIfPresent(testPath));
deleteIfPresent(testPath);

%%
atc_orig = AnnotationTestClass();
ncfile = atc_orig.writeToFile(testPath,shouldOverwriteExisting=1);
atc_read = AnnotationTestClass.annotatedTestClassFromFile(testPath);
assert(atc_orig.isequal(atc_read),"objects not equal");

%%
atc_orig = AnnotationTestClass();
atc_orig.myObjs(2) = AnnotationTestClassB();
ncfile = atc_orig.writeToFile(testPath,shouldOverwriteExisting=1);
atc_read = AnnotationTestClass.annotatedTestClassFromFile(testPath);
assert(atc_orig.isequal(atc_read),"objects not equal");
assert(isequal(atc_read.myObjs(1).y, atc_orig.myObjs(1).y), "Nested numeric properties should round-trip together with their supporting dimensions.");

%%
atcB_orig = AnnotationTestClassB();
ncfile = atcB_orig.writeToFile(testPath,shouldOverwriteExisting=1);
atcB_read = AnnotationTestClass.annotatedTestClassFromFile(testPath);
assert(atcB_orig.isequal(atcB_read),"objects not equal");

%%
atc_empty = AnnotationTestClass();
atc_empty.myObjs = AnnotationTestClassB.empty(0,1);
ncfile = atc_empty.writeToFile(testPath,shouldOverwriteExisting=1);
atc_empty_read = AnnotationTestClass.annotatedTestClassFromFile(testPath);
assert(isempty(atc_empty_read.myObjs), "Empty object arrays should round-trip through annotated persistence.");

%%
atc_orig = AnnotationTestClass();
newProp = CANumericProperty('newVar',{'z'},'', 'nothing');
atc_orig.addPropertyAnnotation(newProp);

%%
path = strcat(tempname, '.nc');
cleanup = onCleanup(@() deleteIfPresent(path));

atc_warn = AnnotationTestClass();
cleanupToken = onCleanup(@() []);
atc_warn.f = @(z) cleanupToken;
atc_warn.writeToFile(path, shouldOverwriteExisting=1);
assert(isa(AnnotationTestClass.annotatedTestClassFromFile(path).f, 'function_handle'), 'Warning-only function handles should still round-trip through file persistence.');

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
