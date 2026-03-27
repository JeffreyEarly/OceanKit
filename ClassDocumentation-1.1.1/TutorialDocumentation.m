classdef TutorialDocumentation < handle
    properties
        sourcePath string
        sourceFile string
        sourceRelativePath string

        title string
        slug string
        description string
        nav_order (1,1) double = NaN

        sections struct = struct("Heading", {}, "Blocks", {})
        figures struct = struct("Name", {}, "Caption", {}, "RelativePath", {})

        buildFolder string
        websiteFolder string
        websiteRootURL string
        executionPaths string = string.empty(0,1)
        previousBuildFolder string = ""

        pathOfOutputFile string
        pathOfAssetFolderOnHardDrive string
        pathOfPreviousAssetFolderOnHardDrive string = ""
        assetPagePrefix string
    end

    methods
        function self = TutorialDocumentation(sourcePath, options)
            arguments
                sourcePath {mustBeTextScalar}
                options.buildFolder {mustBeTextScalar}
                options.websiteFolder {mustBeTextScalar} = "tutorials"
                options.websiteRootURL {mustBeTextScalar} = ""
                options.executionPaths string = string.empty(0,1)
                options.sourceRoot {mustBeTextScalar} = ""
                options.previousBuildFolder {mustBeTextScalar} = ""
            end

            sourcePath = TutorialDocumentation.canonicalPath(string(sourcePath));
            parsedTutorial = TutorialDocumentation.parseSourceFile(sourcePath);

            self.sourcePath = string(sourcePath);
            self.sourceFile = parsedTutorial.SourceFile;
            self.title = parsedTutorial.Title;
            self.slug = parsedTutorial.Slug;
            self.description = parsedTutorial.Description;
            self.nav_order = parsedTutorial.NavOrder;
            self.sections = parsedTutorial.Sections;

            self.buildFolder = TutorialDocumentation.canonicalPath(string(options.buildFolder));
            self.websiteFolder = string(options.websiteFolder);
            self.websiteRootURL = string(options.websiteRootURL);
            self.executionPaths = reshape(string(options.executionPaths), [], 1);
            if string(options.previousBuildFolder) ~= ""
                self.previousBuildFolder = TutorialDocumentation.canonicalPath(string(options.previousBuildFolder));
            end

            self.pathOfOutputFile = fullfile(self.buildFolder, self.websiteFolder, self.slug + ".md");
            self.pathOfAssetFolderOnHardDrive = fullfile(self.buildFolder, self.websiteFolder, self.slug);
            self.assetPagePrefix = "./" + self.slug + "/";
            if self.previousBuildFolder ~= ""
                self.pathOfPreviousAssetFolderOnHardDrive = ...
                    fullfile(self.previousBuildFolder, self.websiteFolder, self.slug);
            end

            if string(options.sourceRoot) ~= ""
                self.sourceRelativePath = TutorialDocumentation.relativePathFromRoot( ...
                    self.sourcePath, TutorialDocumentation.canonicalPath(string(options.sourceRoot)));
            else
                self.sourceRelativePath = self.sourceFile;
            end
        end

        function writeToFile(self)
            runtime = TutorialBuildRuntime(self.pathOfAssetFolderOnHardDrive, ...
                assetPagePrefix=self.assetPagePrefix, ...
                comparisonAssetFolder=self.pathOfPreviousAssetFolderOnHardDrive);

            oldPath = path;
            pathCleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
            if ~isempty(self.executionPaths)
                executionPaths = cellfun(@char, cellstr(self.executionPaths), 'UniformOutput', false);
                addpath(executionPaths{:});
            end

            oldFigureVisibility = get(groot, "defaultFigureVisible");
            figureCleanup = onCleanup(@() set(groot, "defaultFigureVisible", oldFigureVisibility)); %#ok<NASGU>
            set(groot, "defaultFigureVisible", "off");

            tutorialFigureCapture = @(name, varargin) runtime.captureFigure(name, varargin{:}); %#ok<NASGU>
            close all force hidden
            try
                run(char(self.sourcePath));
            catch ME
                close all force hidden
                buildError = MException("TutorialDocumentation:TutorialExecutionFailed", ...
                    "Failed while building tutorial '%s' from '%s'.", ...
                    self.title, self.sourcePath);
                throw(addCause(buildError, ME));
            end

            self.figures = runtime.getFigureRecords();
            self.writeMarkdownPage();
            close all force hidden
        end
    end

    methods (Static)
        function tutorialDocumentation = documentationFromSourceFiles(sourceFiles, options)
            arguments
                sourceFiles
                options.buildFolder {mustBeTextScalar}
                options.websiteFolder {mustBeTextScalar} = "tutorials"
                options.websiteRootURL {mustBeTextScalar} = ""
                options.executionPaths string = string.empty(0,1)
                options.sourceRoot {mustBeTextScalar} = ""
                options.previousBuildFolder {mustBeTextScalar} = ""
            end

            sourceFiles = reshape(string(sourceFiles), [], 1);
            tutorialDocumentation = TutorialDocumentation.empty(numel(sourceFiles), 0);
            for iSource = 1:numel(sourceFiles)
                tutorialDocumentation(iSource) = TutorialDocumentation(sourceFiles(iSource), ...
                    buildFolder=options.buildFolder, ...
                    websiteFolder=options.websiteFolder, ...
                    websiteRootURL=options.websiteRootURL, ...
                    executionPaths=options.executionPaths, ...
                    sourceRoot=options.sourceRoot, ...
                    previousBuildFolder=options.previousBuildFolder);
            end

            tutorialDocumentation = TutorialDocumentation.assignNavOrder(tutorialDocumentation);
        end

        function writeMarkdownIndex(tutorialDocumentation, options)
            arguments
                tutorialDocumentation
                options.buildFolder {mustBeTextScalar}
                options.websiteFolder {mustBeTextScalar} = "tutorials"
                options.nav_order (1,1) double = 4
                options.title {mustBeTextScalar} = "Tutorials"
                options.description {mustBeTextScalar} = ...
                    "These examples are written as plain MATLAB scripts and are rendered into website pages during the documentation build."
            end

            tutorialFolder = fullfile(string(options.buildFolder), string(options.websiteFolder));
            if ~isfolder(tutorialFolder)
                mkdir(tutorialFolder);
            end

            indexPath = fullfile(tutorialFolder, "index.md");
            fileID = fopen(indexPath, "w");
            assert(fileID ~= -1, "Could not open tutorials index for writing.");

            fprintf(fileID, "---\n");
            fprintf(fileID, "layout: default\n");
            fprintf(fileID, "title: %s\n", char(options.title));
            fprintf(fileID, "nav_order: %d\n", options.nav_order);
            fprintf(fileID, "has_children: true\n");
            fprintf(fileID, "mathjax: true\n");
            fprintf(fileID, "permalink: /tutorials\n");
            fprintf(fileID, "---\n\n");
            fprintf(fileID, "# %s\n\n", char(options.title));
            fprintf(fileID, "%s\n\n", char(options.description));

            if isempty(tutorialDocumentation)
                fprintf(fileID, "No tutorials are currently available.\n");
                fclose(fileID);
                return;
            end

            fprintf(fileID, "## Available Tutorials\n\n");
            for iTutorial = 1:numel(tutorialDocumentation)
                fprintf(fileID, "- [%s](./%s)\n", ...
                    char(tutorialDocumentation(iTutorial).title), ...
                    char(tutorialDocumentation(iTutorial).slug));
                fprintf(fileID, "  %s\n", char(tutorialDocumentation(iTutorial).description));
            end

            fclose(fileID);
        end
    end

    methods (Access = private)
        function writeMarkdownPage(self)
            if numel(self.figures) ~= TutorialDocumentation.countFigureMarkers(self.sections)
                error("TutorialDocumentation:FigureMismatch", ...
                    "Tutorial '%s' registered %d figures but contains %d tutorialFigureCapture markers.", ...
                    self.sourceFile, numel(self.figures), TutorialDocumentation.countFigureMarkers(self.sections));
            end

            fileID = fopen(self.pathOfOutputFile, "w");
            assert(fileID ~= -1, "Could not open tutorial page for writing.");

            fprintf(fileID, "---\n");
            fprintf(fileID, "layout: default\n");
            fprintf(fileID, "title: %s\n", char(self.title));
            fprintf(fileID, "parent: Tutorials\n");
            fprintf(fileID, "nav_order: %d\n", self.nav_order);
            fprintf(fileID, "mathjax: true\n");
            fprintf(fileID, "permalink: /tutorials/%s\n", char(self.slug));
            fprintf(fileID, "---\n\n");

            fprintf(fileID, "# %s\n\n", char(self.title));
            fprintf(fileID, "%s\n\n", char(self.description));
            fprintf(fileID, "Source: `%s`\n\n", char(self.sourceRelativePath));

            for iSection = 1:numel(self.sections)
                fprintf(fileID, "## %s\n\n", char(self.sections(iSection).Heading));
                blocks = self.sections(iSection).Blocks;
                for iBlock = 1:numel(blocks)
                    block = blocks(iBlock);
                    switch block.Type
                        case "text"
                            fprintf(fileID, "%s\n\n", char(block.Text));
                        case "code"
                            fprintf(fileID, "```matlab\n%s\n```\n\n", char(block.Text));
                        case "figure"
                            figureRecord = self.figures(block.FigureIndex);
                            altText = figureRecord.Caption;
                            if altText == ""
                                altText = figureRecord.Name;
                            end
                            fprintf(fileID, "![%s](%s)\n\n", char(altText), char(figureRecord.RelativePath));
                            if figureRecord.Caption ~= ""
                                fprintf(fileID, "*%s*\n\n", char(figureRecord.Caption));
                            end
                        otherwise
                            error("TutorialDocumentation:UnknownBlockType", ...
                                "Unknown tutorial block type '%s'.", block.Type);
                    end
                end
            end

            fclose(fileID);
        end
    end

    methods (Static, Access = private)
        function tutorialDocumentation = assignNavOrder(tutorialDocumentation)
            if isempty(tutorialDocumentation)
                return;
            end

            sortTable = table((1:numel(tutorialDocumentation))', [tutorialDocumentation.nav_order]', ...
                string({tutorialDocumentation.title})', ...
                VariableNames={'OriginalIndex', 'NavOrder', 'Title'});
            missingOrder = isnan(sortTable.NavOrder);
            sortTable.NavOrder(missingOrder) = inf;
            sortTable = sortrows(sortTable, {'NavOrder', 'Title', 'OriginalIndex'});
            tutorialDocumentation = tutorialDocumentation(sortTable.OriginalIndex);

            for iTutorial = 1:numel(tutorialDocumentation)
                tutorialDocumentation(iTutorial).nav_order = iTutorial;
            end
        end

        function parsedTutorial = parseSourceFile(sourcePath)
            parsedTutorial = struct( ...
                "SourceFile", "", ...
                "Title", "", ...
                "Slug", "", ...
                "Description", "", ...
                "NavOrder", NaN, ...
                "Sections", struct("Heading", {}, "Blocks", {}));

            sourceText = string(fileread(sourcePath));
            sourceText = replace(sourceText, sprintf("\r\n"), newline);
            sourceText = replace(sourceText, sprintf("\r"), newline);
            lines = splitlines(sourceText);

            metadataLines = strings(0,1);
            sections = struct("Heading", {}, "Blocks", {});
            currentHeading = "";
            currentBlocks = struct("Type", {}, "Text", {}, "FigureIndex", {});
            textBuffer = strings(0,1);
            codeBuffer = strings(0,1);
            nextFigureIndex = 0;
            inMetadataSection = false;

            for iLine = 1:numel(lines)
                line = lines(iLine);

                headingMatch = regexp(char(line), '^\s*%%\s*(.*)$', 'tokens', 'once');
                if ~isempty(headingMatch)
                    [currentBlocks, textBuffer, codeBuffer] = TutorialDocumentation.flushBuffers(currentBlocks, textBuffer, codeBuffer);
                    if strlength(currentHeading) > 0 && ~strcmpi(currentHeading, "Tutorial Metadata")
                        sections(end+1) = struct("Heading", currentHeading, "Blocks", currentBlocks); %#ok<AGROW>
                    end

                    currentHeading = string(strtrim(headingMatch{1}));
                    currentBlocks = struct("Type", {}, "Text", {}, "FigureIndex", {});
                    textBuffer = strings(0,1);
                    codeBuffer = strings(0,1);
                    inMetadataSection = strcmpi(currentHeading, "Tutorial Metadata");
                    continue;
                end

                if ~isempty(regexp(char(line), '^\s*function\b', 'once'))
                    break;
                end

                if inMetadataSection
                    metadataLines(end+1) = line; %#ok<AGROW>
                    continue;
                end

                if strlength(currentHeading) == 0
                    continue;
                end

                commentMatch = regexp(char(line), '^\s*%(?!%)\s?(.*)$', 'tokens', 'once');
                if ~isempty(commentMatch)
                    [currentBlocks, codeBuffer] = TutorialDocumentation.flushCodeBuffer(currentBlocks, codeBuffer);
                    textBuffer(end+1) = string(commentMatch{1}); %#ok<AGROW>
                    continue;
                end

                [currentBlocks, textBuffer] = TutorialDocumentation.flushTextBuffer(currentBlocks, textBuffer);
                if ~isempty(regexp(char(line), 'tutorialFigureCapture\(', 'once'))
                    [currentBlocks, codeBuffer] = TutorialDocumentation.flushCodeBuffer(currentBlocks, codeBuffer);
                    nextFigureIndex = nextFigureIndex + 1;
                    currentBlocks(end+1) = struct( ...
                        "Type", "figure", ...
                        "Text", "", ...
                        "FigureIndex", nextFigureIndex); %#ok<AGROW>
                    continue;
                end

                codeBuffer(end+1) = line; %#ok<AGROW>
            end

            [currentBlocks, textBuffer, codeBuffer] = TutorialDocumentation.flushBuffers(currentBlocks, textBuffer, codeBuffer);
            if strlength(currentHeading) > 0 && ~strcmpi(currentHeading, "Tutorial Metadata")
                sections(end+1) = struct("Heading", currentHeading, "Blocks", currentBlocks); %#ok<AGROW>
            end

            metadata = TutorialDocumentation.parseTutorialMetadata(metadataLines, TutorialDocumentation.fileNameFromPath(sourcePath));
            parsedTutorial.SourceFile = TutorialDocumentation.fileNameFromPath(sourcePath);
            parsedTutorial.Title = metadata.Title;
            parsedTutorial.Slug = metadata.Slug;
            parsedTutorial.Description = metadata.Description;
            parsedTutorial.NavOrder = metadata.NavOrder;
            parsedTutorial.Sections = sections;
        end

        function metadata = parseTutorialMetadata(metadataLines, sourceFile)
            metadata = struct("Title", "", "Slug", "", "Description", "", "NavOrder", NaN);

            for iLine = 1:numel(metadataLines)
                metadataMatch = regexp(char(metadataLines(iLine)), ...
                    '^\s*%\s*(?<Key>[A-Za-z]+)\s*:\s*(?<Value>.*)$', ...
                    'names', 'once');
                if isempty(metadataMatch)
                    continue;
                end

                key = lower(string(metadataMatch.Key));
                value = string(strtrim(metadataMatch.Value));
                switch key
                    case "title"
                        metadata.Title = value;
                    case "slug"
                        metadata.Slug = TutorialDocumentation.slugify(value);
                    case "description"
                        metadata.Description = value;
                    case "navorder"
                        metadata.NavOrder = str2double(value);
                    otherwise
                        error("TutorialDocumentation:UnknownMetadataKey", ...
                            "Unknown tutorial metadata key '%s' in '%s'.", key, sourceFile);
                end
            end

            if metadata.Title == ""
                error("TutorialDocumentation:MissingTitle", ...
                    "Tutorial '%s' is missing a Title entry in the Tutorial Metadata section.", ...
                    sourceFile);
            end
            if metadata.Description == ""
                error("TutorialDocumentation:MissingDescription", ...
                    "Tutorial '%s' is missing a Description entry in the Tutorial Metadata section.", ...
                    sourceFile);
            end
            if metadata.Slug == ""
                metadata.Slug = TutorialDocumentation.slugify(erase(sourceFile, ".m"));
            end
        end

        function count = countFigureMarkers(sections)
            count = 0;
            for iSection = 1:numel(sections)
                blockTypes = string({sections(iSection).Blocks.Type});
                count = count + sum(blockTypes == "figure");
            end
        end

        function [blocks, textBuffer, codeBuffer] = flushBuffers(blocks, textBuffer, codeBuffer)
            [blocks, textBuffer] = TutorialDocumentation.flushTextBuffer(blocks, textBuffer);
            [blocks, codeBuffer] = TutorialDocumentation.flushCodeBuffer(blocks, codeBuffer);
        end

        function [blocks, textBuffer] = flushTextBuffer(blocks, textBuffer)
            textBuffer = TutorialDocumentation.trimBlankEdges(textBuffer);
            if isempty(textBuffer)
                textBuffer = strings(0,1);
                return;
            end

            blocks(end+1) = struct("Type", "text", "Text", join(textBuffer, newline), "FigureIndex", NaN); %#ok<AGROW>
            textBuffer = strings(0,1);
        end

        function [blocks, codeBuffer] = flushCodeBuffer(blocks, codeBuffer)
            codeBuffer = TutorialDocumentation.trimBlankEdges(codeBuffer);
            if isempty(codeBuffer)
                codeBuffer = strings(0,1);
                return;
            end

            blocks(end+1) = struct("Type", "code", "Text", join(codeBuffer, newline), "FigureIndex", NaN); %#ok<AGROW>
            codeBuffer = strings(0,1);
        end

        function lines = trimBlankEdges(lines)
            if isempty(lines)
                return;
            end

            isBlank = @(value) strlength(strtrim(string(value))) == 0;
            while ~isempty(lines) && isBlank(lines(1))
                lines(1) = [];
            end
            while ~isempty(lines) && isBlank(lines(end))
                lines(end) = [];
            end
        end

        function slug = slugify(textValue)
            slug = lower(strtrim(string(textValue)));
            slug = regexprep(slug, "[^a-z0-9]+", "-");
            slug = regexprep(slug, "^-+|-+$", "");
        end

        function fileName = fileNameFromPath(filePath)
            [~, name, ext] = fileparts(char(filePath));
            fileName = string(name) + string(ext);
        end

        function relativePath = relativePathFromRoot(filePath, rootPath)
            normalizedFilePath = replace(string(filePath), "\", "/");
            normalizedRootPath = replace(string(rootPath), "\", "/");
            if ~endsWith(normalizedRootPath, "/")
                normalizedRootPath = normalizedRootPath + "/";
            end

            if startsWith(normalizedFilePath, normalizedRootPath)
                relativePath = extractAfter(normalizedFilePath, strlength(normalizedRootPath));
            else
                relativePath = string(filePath);
            end
        end

        function pathValue = canonicalPath(pathValue)
            pathValue = string(java.io.File(char(pathValue)).getCanonicalPath());
        end
    end
end
