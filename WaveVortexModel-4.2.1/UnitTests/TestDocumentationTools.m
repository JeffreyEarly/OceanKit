classdef TestDocumentationTools < matlab.unittest.TestCase
    properties
        temporaryFolder
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(repositoryRoot,"tools")));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.temporaryFolder = string(fixture.Folder);
        end
    end

    methods (Test,TestTags="full")
        function internalRouteVariantsResolve(testCase)
            root = fullfile(testCase.temporaryFolder,"site");
            mkdir(fullfile(root,"guide"));
            mkdir(fullfile(root,"assets"));
            testCase.writeText(fullfile(root,"index.md"),[
                "[root](/target.html)"
                "[relative](target.md)"
                "[directory](guide/)"
                "[permalink](/special/)"
                "[fragment](target.md#section-name)"
                "![asset](assets/image.png)"
                "[external](https://example.com/path)"
                "<a href=""/target.html"">HTML link</a>"
                "<img src=""assets/image.png"">"
                "<a href='target.md'>Single-quoted HTML link</a>"
                ]);
            testCase.writeText(fullfile(root,"target.md"),"# Section name");
            testCase.writeText(fullfile(root,"guide","index.md"),"# Guide");
            testCase.writeText(fullfile(root,"permalink.md"),[
                "---"
                "permalink: /special/"
                "---"
                "# Permalink"
                ]);
            testCase.writeText(fullfile(root,"assets","image.png"),"fixture");

            report = validateWebsiteDocumentation(root,ShouldFail=false,ShouldCheckHierarchy=false,ShouldCheckGeneratedContent=false);
            testCase.verifyTrue(report.IsValid,strjoin(report.Diagnostics,newline));
        end

        function invalidRoutesAreReported(testCase)
            root = fullfile(testCase.temporaryFolder,"site");
            mkdir(root);
            testCase.writeText(fullfile(root,"Target.md"),"# Existing target");
            testCase.writeText(fullfile(root,"index.md"),[
                "[missing](missing.md)"
                "[case mismatch](target.md)"
                "[malformed](missing.md"
                "<a href=""unterminated>HTML link</a>"
                ]);

            report = validateWebsiteDocumentation(root,ShouldFail=false,ShouldCheckHierarchy=false,ShouldCheckGeneratedContent=false);
            testCase.verifyFalse(report.IsValid);
            testCase.verifyTrue(any(contains(report.Diagnostics,"missing.md")));
            testCase.verifyTrue(any(contains(report.Diagnostics,"target.md")));
            testCase.verifyTrue(any(contains(report.Diagnostics,"malformed Markdown")));
            testCase.verifyTrue(any(contains(report.Diagnostics,"malformed HTML")));
        end

        function missingFragmentIsReported(testCase)
            root = fullfile(testCase.temporaryFolder,"site");
            mkdir(root);
            testCase.writeText(fullfile(root,"index.md"),"[fragment](target.md#missing)");
            testCase.writeText(fullfile(root,"target.md"),"# Present");

            report = validateWebsiteDocumentation(root,ShouldFail=false,ShouldCheckHierarchy=false,ShouldCheckGeneratedContent=false);
            testCase.verifyEqual(numel(report.Diagnostics),1);
            testCase.verifySubstring(report.Diagnostics,"unresolved fragment #missing");
        end

        function generatedHierarchyIsValidated(testCase)
            root = fullfile(testCase.temporaryFolder,"site");
            classFolder = fullfile(root,"classes","example");
            mkdir(classFolder);
            testCase.writeText(fullfile(classFolder,"index.md"),[
                "---"
                "title: ""ExampleClass"""
                "parent: Transforms"
                "---"
                "## Declaration"
                ]);
            methodPath = fullfile(classFolder,"method.md");
            testCase.writeText(methodPath,[
                "---"
                "title: method"
                "parent: ExampleClass"
                "grand_parent: Transforms"
                "---"
                ]);

            report = validateWebsiteDocumentation(root,ShouldFail=false,ShouldCheckGeneratedContent=false);
            testCase.verifyTrue(report.IsValid,strjoin(report.Diagnostics,newline));

            testCase.writeText(methodPath,[
                "---"
                "title: method"
                "parent: ExampleClass"
                "grand_parent: Classes"
                "---"
                ]);
            report = validateWebsiteDocumentation(root,ShouldFail=false,ShouldCheckGeneratedContent=false);
            testCase.verifyTrue(any(contains(report.Diagnostics,"grand_parent does not match")));
        end

        function generatedContentPreservesRho(testCase)
            root = fullfile(testCase.temporaryFolder,"site");
            transformFolder = fullfile(root,"classes","transforms","wvtransform");
            mkdir(transformFolder);
            indexPath = fullfile(transformFolder,"index.md");
            testCase.writeText(indexPath,"Density is $$\rho$$.");

            report = validateWebsiteDocumentation(root,ShouldFail=false,ShouldCheckHierarchy=false);
            testCase.verifyTrue(report.IsValid,strjoin(report.Diagnostics,newline));

            testCase.writeText(indexPath,"Density is $$ho$$.");
            report = validateWebsiteDocumentation(root,ShouldFail=false,ShouldCheckHierarchy=false);
            testCase.verifyFalse(report.IsValid);
            testCase.verifyTrue(any(contains(report.Diagnostics,"historical truncated")));
        end

        function validFrontMatterAndMathMarkupPass(testCase)
            root = fullfile(testCase.temporaryFolder,"valid-source");
            mkdir(root);
            testCase.writeText(fullfile(root,"index.md"),[
                "---"
                "layout: default"
                "title: ""Guide: details"""
                "---"
                "# Guide"
                "The literal delimiter `$$` may be discussed in code."
                "Inline $$5 \\leq \\lvert x\\rvert \\leq 85$$ math."
                "$$"
                "\\begin{align}"
                "x &= y"
                "\\end{align}"
                "$$"
                ]);

            report = validateWebsiteDocumentation(root,ShouldFail=false,ShouldCheckHierarchy=false,ShouldCheckGeneratedContent=false);
            testCase.verifyTrue(report.IsValid,strjoin(report.Diagnostics,newline));
        end

        function malformedFrontMatterAndMathAreReported(testCase)
            root = fullfile(testCase.temporaryFolder,"invalid-source");
            mkdir(root);
            testCase.writeText(fullfile(root,"unterminated.md"),["---"; "title: Missing end"]);
            testCase.writeText(fullfile(root,"unquoted.md"),["---"; "title: Guide: details"; "---"]);
            testCase.writeText(fullfile(root,"malformed-quote.md"),["---"; "title: ""Missing end"; "---"]);
            testCase.writeText(fullfile(root,"duplicate.md"),["---"; "title: First"; "title: Second"; "---"]);
            testCase.writeText(fullfile(root,"unbalanced-display.md"),["# Display"; "$$"; "x=y"]);
            testCase.writeText(fullfile(root,"display-markdown.md"),["# Display"; "$$"; "[link](target.md) and `code`"; "## Heading"; "$$"]);
            testCase.writeText(fullfile(root,"inline-bar.md"),"Inline $$|x|$$ math.");

            report = validateWebsiteDocumentation(root,ShouldFail=false,ShouldCheckHierarchy=false,ShouldCheckGeneratedContent=false);
            expected = [
                "unterminated front matter"
                "must be quoted"
                "malformed quoted front matter value"
                "duplicate front matter key"
                "unbalanced display-math"
                "Markdown link or image appears inside"
                "Markdown backticks appear inside"
                "Markdown heading appears inside"
                "unescaped vertical bar"
                ];
            for diagnostic = expected'
                testCase.verifyTrue(any(contains(report.Diagnostics,diagnostic)),"Missing diagnostic: " + diagnostic);
            end
        end

        function renderedWebsiteRejectsRawMarkup(testCase)
            root = fullfile(testCase.temporaryFolder,"rendered-site");
            mkdir(root);
            testCase.writeText(fullfile(root,"index.html"),'<main><h1>Home</h1><p>Rendered $$x$$ content.</p><code>`[example](target)`</code></main>');
            cleanReport = validateRenderedWebsite(root,ShouldFail=false);
            testCase.verifyTrue(cleanReport.IsValid,strjoin(cleanReport.Diagnostics,newline));

            testCase.writeText(fullfile(root,"missing-main.html"),'<h1>Missing main</h1>');
            testCase.writeText(fullfile(root,"missing-heading.html"),'<main><p>Missing heading</p></main>');
            testCase.writeText(fullfile(root,"raw-math.html"),'<main><h1>Math</h1><p>Unbalanced $$x</p></main>');
            testCase.writeText(fullfile(root,"table-math.html"),'<main><h1>Table</h1><table><tr><td>$$|x</td><td>|$$</td></tr></table></main>');
            testCase.writeText(fullfile(root,"raw-link.html"),'<main><h1>Link</h1><p>[link](target.html)</p></main>');
            testCase.writeText(fullfile(root,"raw-fence.html"),'<main><h1>Fence</h1><p>```matlab</p></main>');
            testCase.writeText(fullfile(root,"raw-backtick.html"),'<main><h1>Code</h1><p>`code`</p></main>');

            report = validateRenderedWebsite(root,ShouldFail=false);
            expected = [
                "missing main content element"
                "missing level-one heading"
                "unbalanced MathJax delimiters"
                "MathJax expression is split across rendered table cells"
                "raw Markdown link or image"
                "raw fenced-code marker"
                "raw Markdown backtick"
                ];
            for diagnostic = expected'
                testCase.verifyTrue(any(contains(report.Diagnostics,diagnostic)),"Missing diagnostic: " + diagnostic);
            end
        end

        function treeComparisonClassifiesDifferences(testCase)
            reference = fullfile(testCase.temporaryFolder,"reference");
            candidate = fullfile(testCase.temporaryFolder,"candidate");
            mkdir(reference);
            mkdir(candidate);
            testCase.writeText(fullfile(reference,"identical.md"),"same");
            testCase.writeText(fullfile(candidate,"identical.md"),"same");
            testCase.writeText(fullfile(reference,"whitespace.md"),"alpha beta");
            testCase.writeText(fullfile(candidate,"whitespace.md"),"alpha    beta");
            testCase.writeText(fullfile(reference,"modified.md"),"alpha");
            testCase.writeText(fullfile(candidate,"modified.md"),"beta");
            testCase.writeText(fullfile(reference,"removed.md"),"removed");
            testCase.writeText(fullfile(candidate,"added.md"),"added");

            comparison = compareDocumentationTrees(reference,candidate);
            testCase.verifyFalse(comparison.IsEqual);
            testCase.verifyEqual(comparison.Added,"added.md");
            testCase.verifyEqual(comparison.Removed,"removed.md");
            testCase.verifyEqual(comparison.Modified,["modified.md"; "whitespace.md"]);
            testCase.verifyEqual(comparison.WhitespaceOnly,"whitespace.md");
            testCase.verifyEqual(comparison.Substantive,"modified.md");
        end

        function successfulReplacementRemovesStalePages(testCase)
            destination = fullfile(testCase.temporaryFolder,"docs");
            staging = fullfile(testCase.temporaryFolder,"staging");
            mkdir(destination);
            mkdir(staging);
            testCase.writeText(fullfile(destination,"stale.md"),"stale");
            testCase.writeText(fullfile(staging,"current.md"),"current");

            replaceDocumentationTree(staging,destination);

            testCase.verifyFalse(isfile(fullfile(destination,"stale.md")));
            testCase.verifyTrue(isfile(fullfile(destination,"current.md")));
        end

        function failedReplacementRestoresOriginalTree(testCase)
            destination = fullfile(testCase.temporaryFolder,"docs");
            original = fullfile(testCase.temporaryFolder,"original");
            missingStaging = fullfile(testCase.temporaryFolder,"missing-staging");
            mkdir(destination);
            mkdir(original);
            testCase.writeText(fullfile(destination,"index.md"),"original");
            testCase.writeText(fullfile(original,"index.md"),"original");

            didFail = false;
            try
                replaceDocumentationTree(missingStaging,destination);
            catch
                didFail = true;
            end

            testCase.verifyTrue(didFail);
            comparison = compareDocumentationTrees(original,destination);
            testCase.verifyTrue(comparison.IsEqual);
        end

        function validationFailureLeavesDestinationUnchanged(testCase)
            destination = fullfile(testCase.temporaryFolder,"docs");
            original = fullfile(testCase.temporaryFolder,"original");
            staging = fullfile(testCase.temporaryFolder,"staging");
            mkdir(destination);
            mkdir(original);
            mkdir(staging);
            testCase.writeText(fullfile(destination,"index.md"),"original");
            testCase.writeText(fullfile(original,"index.md"),"original");
            testCase.writeText(fullfile(staging,"index.md"),"[missing](missing.md)");

            testCase.verifyError(@()validateWebsiteDocumentation(staging,ShouldCheckHierarchy=false,ShouldCheckGeneratedContent=false),"WaveVortexModel:InvalidDocumentation");
            comparison = compareDocumentationTrees(original,destination);
            testCase.verifyTrue(comparison.IsEqual);
        end
    end

    methods (Access=private)
        function writeText(~,path,text)
            writelines(text,path);
        end
    end
end
