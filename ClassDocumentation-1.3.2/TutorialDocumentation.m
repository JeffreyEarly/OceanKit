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
        movies struct = struct("Name", {}, "Caption", {}, "RelativePath", {}, "PosterRelativePath", {})
        outputs struct = struct("Caption", {}, "Language", {}, "Text", {})

        buildFolder string
        websiteFolder string
        websiteRootURL string
        executionPaths string = string.empty(0,1)
        previousBuildFolder string = ""
        rebuildTutorials (1,1) logical = false
        sourceHash string

        pathOfOutputFile string
        pathOfBuildStampFile string
        pathOfAssetFolderOnHardDrive string
        pathOfPreservedAssetFolderOnHardDrive string
        pathOfPreviousOutputFile string = ""
        pathOfPreviousBuildStampFile string = ""
        pathOfPreviousAssetFolderOnHardDrive string = ""
        assetPagePrefix string
        preservedAssetDirectoryRelativeToBuildFolder string
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
                options.rebuildTutorials (1,1) logical = false
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
            self.sourceHash = TutorialDocumentation.hashNormalizedSourceText( ...
                TutorialDocumentation.normalizedSourceTextFromFile(sourcePath));

            self.buildFolder = TutorialDocumentation.canonicalPath(string(options.buildFolder));
            self.websiteFolder = string(options.websiteFolder);
            self.websiteRootURL = string(options.websiteRootURL);
            self.executionPaths = reshape(string(options.executionPaths), [], 1);
            self.rebuildTutorials = options.rebuildTutorials;
            if string(options.previousBuildFolder) ~= ""
                self.previousBuildFolder = TutorialDocumentation.canonicalPath(string(options.previousBuildFolder));
            end

            self.pathOfOutputFile = fullfile(self.buildFolder, self.websiteFolder, self.slug + ".md");
            self.pathOfBuildStampFile = fullfile(self.buildFolder, self.websiteFolder, self.slug + ".build-stamp.json");
            self.pathOfAssetFolderOnHardDrive = fullfile(self.buildFolder, self.websiteFolder, self.slug);
            self.pathOfPreservedAssetFolderOnHardDrive = fullfile(self.pathOfAssetFolderOnHardDrive, "preserved");
            self.assetPagePrefix = "./" + self.slug + "/";
            self.preservedAssetDirectoryRelativeToBuildFolder = fullfile(self.websiteFolder, self.slug, "preserved");
            if self.previousBuildFolder ~= ""
                self.pathOfPreviousOutputFile = fullfile(self.previousBuildFolder, self.websiteFolder, self.slug + ".md");
                self.pathOfPreviousBuildStampFile = ...
                    fullfile(self.previousBuildFolder, self.websiteFolder, self.slug + ".build-stamp.json");
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
            buildStampText = self.tutorialBuildStampText();
            if self.tryReuseWholeTutorialFromPreviousBuild(buildStampText)
                return;
            end

            runtime = TutorialBuildRuntime(self.pathOfAssetFolderOnHardDrive, ...
                assetPagePrefix=self.assetPagePrefix, ...
                comparisonAssetFolder=self.pathOfPreviousAssetFolderOnHardDrive, ...
                sourcePath=self.sourcePath);

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
            tutorialMovieCapture = @(name, varargin) runtime.captureMovie(name, varargin{:}); %#ok<NASGU>
            tutorialOutputCapture = @(source, varargin) runtime.captureOutput(source, varargin{:}); %#ok<NASGU>
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
            self.movies = runtime.getMovieRecords();
            self.outputs = runtime.getOutputRecords();
            self.writeMarkdownPage();
            TutorialDocumentation.writeTextFile(self.pathOfBuildStampFile, buildStampText);
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
                options.rebuildTutorials (1,1) logical = false
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
                    previousBuildFolder=options.previousBuildFolder, ...
                    rebuildTutorials=options.rebuildTutorials);
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
        function buildStampText = tutorialBuildStampText(self)
            stamp = struct( ...
                "Version", 1, ...
                "Slug", char(self.slug), ...
                "SourceHash", char(self.sourceHash));
            buildStampText = string(jsonencode(stamp));
        end

        function didReuse = tryReuseWholeTutorialFromPreviousBuild(self, buildStampText)
            didReuse = false;
            if self.rebuildTutorials || self.previousBuildFolder == ""
                return;
            end
            if ~isfile(self.pathOfPreviousOutputFile) || ~isfile(self.pathOfPreviousBuildStampFile)
                return;
            end

            previousBuildStampText = string(fileread(self.pathOfPreviousBuildStampFile));
            if previousBuildStampText ~= buildStampText
                return;
            end

            TutorialDocumentation.ensureParentFolderExists(self.pathOfOutputFile);
            copyfile(self.pathOfPreviousOutputFile, self.pathOfOutputFile, "f");

            if isfolder(self.pathOfPreviousAssetFolderOnHardDrive)
                if isfolder(self.pathOfAssetFolderOnHardDrive)
                    rmdir(self.pathOfAssetFolderOnHardDrive, "s");
                end
                TutorialDocumentation.ensureParentFolderExists(self.pathOfAssetFolderOnHardDrive);
                copyfile(self.pathOfPreviousAssetFolderOnHardDrive, self.pathOfAssetFolderOnHardDrive);
            end

            TutorialDocumentation.writeTextFile(self.pathOfBuildStampFile, buildStampText);
            didReuse = true;
        end

        function writeMarkdownPage(self)
            if numel(self.figures) ~= TutorialDocumentation.countBlockType(self.sections, "figure")
                error("TutorialDocumentation:FigureMismatch", ...
                    "Tutorial '%s' registered %d figures but contains %d tutorialFigureCapture markers.", ...
                    self.sourceFile, numel(self.figures), TutorialDocumentation.countBlockType(self.sections, "figure"));
            end
            if numel(self.movies) ~= TutorialDocumentation.countBlockType(self.sections, "movie")
                error("TutorialDocumentation:MovieMismatch", ...
                    "Tutorial '%s' registered %d movies but contains %d tutorialMovieCapture markers.", ...
                    self.sourceFile, numel(self.movies), TutorialDocumentation.countBlockType(self.sections, "movie"));
            end
            if numel(self.outputs) ~= TutorialDocumentation.countBlockType(self.sections, "output")
                error("TutorialDocumentation:OutputMismatch", ...
                    "Tutorial '%s' registered %d outputs but contains %d tutorialOutputCapture markers.", ...
                    self.sourceFile, numel(self.outputs), TutorialDocumentation.countBlockType(self.sections, "output"));
            end

            outputFolder = fileparts(self.pathOfOutputFile);
            if ~isfolder(outputFolder)
                mkdir(outputFolder);
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
                        case "movie"
                            movieRecord = self.movies(block.MovieIndex);
                            fprintf(fileID, "<video\n");
                            fprintf(fileID, "  controls\n");
                            fprintf(fileID, "  preload=""metadata""\n");
                            if movieRecord.PosterRelativePath ~= ""
                                fprintf(fileID, "  poster=""%s""\n", char(movieRecord.PosterRelativePath));
                            end
                            fprintf(fileID, "  style=""max-width: 100%%; height: auto;"">\n");
                            fprintf(fileID, "  <source src=""%s"" type=""%s"">\n", ...
                                char(movieRecord.RelativePath), ...
                                char(TutorialDocumentation.movieMimeType(movieRecord.RelativePath)));
                            fprintf(fileID, "  Your browser does not support the HTML5 video tag.\n");
                            fprintf(fileID, "</video>\n\n");
                            if movieRecord.Caption ~= ""
                                fprintf(fileID, "*%s*\n\n", char(movieRecord.Caption));
                            end
                        case "output"
                            outputRecord = self.outputs(block.OutputIndex);
                            fprintf(fileID, "```%s\n%s\n```\n\n", char(outputRecord.Language), char(outputRecord.Text));
                            if outputRecord.Caption ~= ""
                                fprintf(fileID, "*%s*\n\n", char(outputRecord.Caption));
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

            sourceText = TutorialDocumentation.normalizedSourceTextFromFile(sourcePath);
            lines = splitlines(sourceText);

            metadataLines = strings(0,1);
            sections = struct("Heading", {}, "Blocks", {});
            currentHeading = "";
            currentBlocks = TutorialDocumentation.emptyBlocks();
            textBuffer = strings(0,1);
            codeBuffer = strings(0,1);
            nextFigureIndex = 0;
            nextMovieIndex = 0;
            nextOutputIndex = 0;
            inMetadataSection = false;

            iLine = 1;
            while iLine <= numel(lines)
                line = lines(iLine);

                headingMatch = regexp(char(line), '^\s*%%\s*(.*)$', 'tokens', 'once');
                if ~isempty(headingMatch)
                    [currentBlocks, textBuffer, codeBuffer] = TutorialDocumentation.flushBuffers(currentBlocks, textBuffer, codeBuffer);
                    if strlength(currentHeading) > 0 && ~strcmpi(currentHeading, "Tutorial Metadata")
                        sections(end+1) = struct("Heading", currentHeading, "Blocks", currentBlocks); %#ok<AGROW>
                    end

                    currentHeading = string(strtrim(headingMatch{1}));
                    currentBlocks = TutorialDocumentation.emptyBlocks();
                    textBuffer = strings(0,1);
                    codeBuffer = strings(0,1);
                    inMetadataSection = strcmpi(currentHeading, "Tutorial Metadata");
                    iLine = iLine + 1;
                    continue;
                end

                if ~isempty(regexp(char(line), '^\s*function(?:\s|$)', 'once'))
                    break;
                end

                if inMetadataSection
                    metadataLines(end+1) = line; %#ok<AGROW>
                    iLine = iLine + 1;
                    continue;
                end

                if strlength(currentHeading) == 0
                    iLine = iLine + 1;
                    continue;
                end

                commentMatch = regexp(char(line), '^\s*%(?!%)\s?(.*)$', 'tokens', 'once');
                if ~isempty(commentMatch)
                    [currentBlocks, codeBuffer] = TutorialDocumentation.flushCodeBuffer(currentBlocks, codeBuffer);
                    textBuffer(end+1) = string(commentMatch{1}); %#ok<AGROW>
                    iLine = iLine + 1;
                    continue;
                end

                captureType = TutorialDocumentation.captureStatementType(line);
                if captureType ~= ""
                    [currentBlocks, textBuffer] = TutorialDocumentation.flushTextBuffer(currentBlocks, textBuffer);
                    codeBuffer = TutorialDocumentation.stripTrailingCaptureGuard(codeBuffer);
                    [currentBlocks, codeBuffer] = TutorialDocumentation.flushCodeBuffer(currentBlocks, codeBuffer);

                    switch captureType
                        case "figure"
                            nextFigureIndex = nextFigureIndex + 1;
                            currentBlocks(end+1) = TutorialDocumentation.assetBlock("figure", nextFigureIndex, NaN, NaN); %#ok<AGROW>
                        case "movie"
                            nextMovieIndex = nextMovieIndex + 1;
                            currentBlocks(end+1) = TutorialDocumentation.assetBlock("movie", NaN, nextMovieIndex, NaN); %#ok<AGROW>
                        case "output"
                            nextOutputIndex = nextOutputIndex + 1;
                            currentBlocks(end+1) = TutorialDocumentation.assetBlock("output", NaN, NaN, nextOutputIndex); %#ok<AGROW>
                    end

                    iLine = TutorialDocumentation.advancePastCaptureStatement(lines, iLine);
                    if iLine <= numel(lines) && TutorialDocumentation.lineClosesCaptureGuard(lines(iLine))
                        iLine = iLine + 1;
                    end
                    continue;
                end

                [currentBlocks, textBuffer] = TutorialDocumentation.flushTextBuffer(currentBlocks, textBuffer);
                codeBuffer(end+1) = line; %#ok<AGROW>
                iLine = iLine + 1;
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

        function count = countBlockType(sections, blockType)
            count = 0;
            for iSection = 1:numel(sections)
                blockTypes = string({sections(iSection).Blocks.Type});
                count = count + sum(blockTypes == blockType);
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

            blocks(end+1) = TutorialDocumentation.textOrCodeBlock("text", join(textBuffer, newline)); %#ok<AGROW>
            textBuffer = strings(0,1);
        end

        function [blocks, codeBuffer] = flushCodeBuffer(blocks, codeBuffer)
            codeBuffer = TutorialDocumentation.trimBlankEdges(codeBuffer);
            if isempty(codeBuffer)
                codeBuffer = strings(0,1);
                return;
            end

            blocks(end+1) = TutorialDocumentation.textOrCodeBlock("code", join(codeBuffer, newline)); %#ok<AGROW>
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

        function sourceText = normalizedSourceTextFromFile(sourcePath)
            sourceText = string(fileread(sourcePath));
            sourceText = replace(sourceText, sprintf("\r\n"), newline);
            sourceText = replace(sourceText, sprintf("\r"), newline);
        end

        function hashValue = hashNormalizedSourceText(sourceText)
            hashValue = TutorialDocumentation.hashBytes(unicode2native(char(sourceText), "UTF-8"));
        end

        function hashValue = hashBytes(fileBytes)
            messageDigest = java.security.MessageDigest.getInstance("SHA-256");
            messageDigest.update(int8(fileBytes));
            digestBytes = typecast(messageDigest.digest(), "uint8");
            hexMatrix = dec2hex(digestBytes, 2).';
            hashValue = lower(string(reshape(hexMatrix, 1, [])));
        end

        function ensureParentFolderExists(filePath)
            parentFolder = fileparts(filePath);
            if strlength(parentFolder) == 0 || isfolder(parentFolder)
                return;
            end

            mkdir(parentFolder);
        end

        function writeTextFile(filePath, fileText)
            TutorialDocumentation.ensureParentFolderExists(filePath);
            fileID = fopen(filePath, "w");
            assert(fileID ~= -1, "Could not open file for writing: %s", filePath);
            fwrite(fileID, char(fileText));
            fclose(fileID);
        end

        function blocks = emptyBlocks()
            blocks = struct("Type", {}, "Text", {}, "FigureIndex", {}, "MovieIndex", {}, "OutputIndex", {});
        end

        function block = textOrCodeBlock(blockType, blockText)
            block = struct("Type", blockType, "Text", string(blockText), "FigureIndex", NaN, "MovieIndex", NaN, "OutputIndex", NaN);
        end

        function block = assetBlock(blockType, figureIndex, movieIndex, outputIndex)
            block = struct("Type", blockType, "Text", "", "FigureIndex", figureIndex, "MovieIndex", movieIndex, "OutputIndex", outputIndex);
        end

        function captureType = captureStatementType(line)
            captureType = "";
            if ~isempty(regexp(char(line), 'tutorialFigureCapture\(', 'once'))
                captureType = "figure";
                return;
            end
            if ~isempty(regexp(char(line), 'tutorialMovieCapture\(', 'once'))
                captureType = "movie";
                return;
            end
            if ~isempty(regexp(char(line), 'tutorialOutputCapture\(', 'once'))
                captureType = "output";
            end
        end

        function codeBuffer = stripTrailingCaptureGuard(codeBuffer)
            if isempty(codeBuffer)
                return;
            end

            iLine = numel(codeBuffer);
            while iLine >= 1 && strlength(strtrim(codeBuffer(iLine))) == 0
                iLine = iLine - 1;
            end
            if iLine < 1
                codeBuffer = strings(0,1);
                return;
            end

            if TutorialDocumentation.isCaptureGuardLine(codeBuffer(iLine))
                codeBuffer(iLine:end) = [];
            end
        end

        function tf = isCaptureGuardLine(line)
            tf = ~isempty(regexp(char(line), ...
                '^\s*if\b.*(?:tutorialFigureCapture|tutorialMovieCapture|tutorialOutputCapture).*$', ...
                'once'));
        end

        function tf = lineClosesCaptureGuard(line)
            tf = ~isempty(regexp(char(line), '^\s*end\s*;?\s*$', 'once'));
        end

        function nextLineIndex = advancePastCaptureStatement(lines, startLineIndex)
            nextLineIndex = startLineIndex + 1;
            parenthesisBalance = 0;
            didSeeOpeningParenthesis = false;
            iLine = startLineIndex;
            while iLine <= numel(lines)
                sanitizedLine = TutorialDocumentation.sanitizedLineForParenthesisCounting(lines(iLine));
                parenthesisBalance = parenthesisBalance + count(sanitizedLine, "(") - count(sanitizedLine, ")");
                didSeeOpeningParenthesis = didSeeOpeningParenthesis || contains(sanitizedLine, "(");
                if didSeeOpeningParenthesis && parenthesisBalance <= 0
                    nextLineIndex = iLine + 1;
                    return;
                end
                iLine = iLine + 1;
            end
        end

        function sanitizedLine = sanitizedLineForParenthesisCounting(line)
            lineCharacters = char(line);
            sanitizedCharacters = lineCharacters;
            inSingleQuotes = false;
            inDoubleQuotes = false;

            iCharacter = 1;
            while iCharacter <= numel(lineCharacters)
                currentCharacter = lineCharacters(iCharacter);
                if inSingleQuotes
                    sanitizedCharacters(iCharacter) = ' ';
                    if currentCharacter == ''''
                        if iCharacter < numel(lineCharacters) && lineCharacters(iCharacter + 1) == ''''
                            sanitizedCharacters(iCharacter + 1) = ' ';
                            iCharacter = iCharacter + 1;
                        else
                            inSingleQuotes = false;
                        end
                    end
                elseif inDoubleQuotes
                    sanitizedCharacters(iCharacter) = ' ';
                    if currentCharacter == '"'
                        if iCharacter < numel(lineCharacters) && lineCharacters(iCharacter + 1) == '"'
                            sanitizedCharacters(iCharacter + 1) = ' ';
                            iCharacter = iCharacter + 1;
                        else
                            inDoubleQuotes = false;
                        end
                    end
                else
                    if currentCharacter == '%'
                        sanitizedCharacters(iCharacter:end) = ' ';
                        break;
                    elseif currentCharacter == ''''
                        inSingleQuotes = true;
                        sanitizedCharacters(iCharacter) = ' ';
                    elseif currentCharacter == '"'
                        inDoubleQuotes = true;
                        sanitizedCharacters(iCharacter) = ' ';
                    end
                end
                iCharacter = iCharacter + 1;
            end

            sanitizedLine = string(sanitizedCharacters);
        end

        function mimeType = movieMimeType(relativePath)
            [~, ~, extension] = fileparts(char(relativePath));
            switch lower(string(extension))
                case ".mp4"
                    mimeType = "video/mp4";
                case ".webm"
                    mimeType = "video/webm";
                case ".mov"
                    mimeType = "video/quicktime";
                otherwise
                    mimeType = "video/mp4";
            end
        end
    end
end
