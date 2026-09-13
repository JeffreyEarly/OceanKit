classdef TestPortableQualificationCatalog < matlab.unittest.TestCase
    methods (Test,TestTags="full")
        function archivedEvidenceRequiresExactOwnSlice(testCase)
            root = string(fileparts(fileparts(mfilename("fullpath"))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"tools")));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            target = string(fixture.Folder);
            current = "PortableRuntime/contracts/portable-forcing-compatibility-v1.json";
            archived = "PortableRuntime/qualification/forcing-catalog-through-hydrostatic-v1.json";
            mkdir(fileparts(fullfile(target,current))); mkdir(fileparts(fullfile(target,archived)));
            copyfile(fullfile(root,current),fullfile(target,current));
            copyfile(fullfile(root,archived),fullfile(target,archived));
            original = jsondecode(fileread(fullfile(target,current)));
            for family = ["hydrostatic","stratified-qg"]
                report = jsondecode(fileread(fullfile(root,"PortableRuntime","qualification",family+"-apple-silicon-v1.json")));
                [~,valid] = portableQualificationCatalog(report,family,target); testCase.verifyTrue(valid);
                bad = report; bad.catalogSHA256="stale";
                [~,valid] = portableQualificationCatalog(bad,family,target); testCase.verifyFalse(valid);
                bad = report; bad.catalogPath="../outside.json";
                [~,valid] = portableQualificationCatalog(bad,family,target); testCase.verifyFalse(valid);
                changed = original; row=find(startsWith(string({changed.rows.configuration}),family),1);
                changed.rows(row).matlab.priority=13;
                writeJSON(fullfile(target,current),changed);
                [~,valid] = portableQualificationCatalog(report,family,target); testCase.verifyFalse(valid);
                changed = original; prefix="hydro-"; if family=="stratified-qg", prefix="sqg-"; end
                index=find(startsWith(string({changed.evidence.id}),prefix),1);
                changed.evidence(index).symbol="differentScientificEvidence";
                writeJSON(fullfile(target,current),changed);
                [~,valid] = portableQualificationCatalog(report,family,target); testCase.verifyFalse(valid);
                changed = original; row=find(startsWith(string({changed.rows.configuration}),"boussinesq"),1);
                changed.rows(row).matlab.priority=13;
                writeJSON(fullfile(target,current),changed);
                [~,valid] = portableQualificationCatalog(report,family,target); testCase.verifyTrue(valid);
                writeJSON(fullfile(target,current),original);
            end
        end
    end
end
function writeJSON(path,value)
file=fopen(path,"w"); cleanup=onCleanup(@()fclose(file)); fprintf(file,"%s\n",jsonencode(value));
end
