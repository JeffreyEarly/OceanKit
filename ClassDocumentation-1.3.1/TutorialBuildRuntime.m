classdef TutorialBuildRuntime < handle
    properties
        assetFolder string
        assetPagePrefix string
        comparisonAssetFolder string = ""
        sourcePath string = ""
        figureRecords struct = struct("Name", {}, "Caption", {}, "RelativePath", {})
        movieRecords struct = struct("Name", {}, "Caption", {}, "RelativePath", {}, "PosterRelativePath", {})
        outputRecords struct = struct("Caption", {}, "Language", {}, "Text", {})
    end

    methods
        function self = TutorialBuildRuntime(assetFolder, options)
            arguments
                assetFolder (1,1) string
                options.assetPagePrefix (1,1) string = "./"
                options.comparisonAssetFolder (1,1) string = ""
                options.sourcePath (1,1) string = ""
            end

            self.assetFolder = assetFolder;
            self.assetPagePrefix = options.assetPagePrefix;
            self.comparisonAssetFolder = options.comparisonAssetFolder;
            self.sourcePath = options.sourcePath;
        end

        function captureFigure(self, name, options)
            arguments
                self (1,1) TutorialBuildRuntime
                name {mustBeTextScalar}
                options.Caption string = ""
                options.Figure = gcf
                options.Resolution (1,1) double {mustBePositive} = 200
            end

            figureName = TutorialBuildRuntime.slugify(string(name));
            if any(strcmp(string({self.figureRecords.Name}), figureName))
                error("TutorialBuildRuntime:DuplicateFigureName", ...
                    "Duplicate tutorial figure name '%s'.", figureName);
            end

            TutorialBuildRuntime.ensureFolderExists(self.assetFolder);

            fileName = figureName + ".png";
            targetPath = fullfile(self.assetFolder, fileName);
            temporaryPath = string(tempname()) + ".png";
            exportgraphics(options.Figure, temporaryPath, Resolution=options.Resolution);

            comparisonPath = "";
            if self.comparisonAssetFolder ~= ""
                comparisonPath = fullfile(self.comparisonAssetFolder, fileName);
            end

            if comparisonPath ~= "" && isfile(comparisonPath) ...
                    && TutorialBuildRuntime.imagesMatch(temporaryPath, comparisonPath)
                copyfile(comparisonPath, targetPath, "f");
                delete(temporaryPath);
            else
                movefile(temporaryPath, targetPath, "f");
            end

            self.figureRecords(end+1) = struct( ...
                "Name", figureName, ...
                "Caption", string(options.Caption), ...
                "RelativePath", self.assetPagePrefix + fileName);
        end

        function captureMovie(self, name, options)
            arguments
                self (1,1) TutorialBuildRuntime
                name {mustBeTextScalar}
                options.Caption string = ""
                options.Poster {mustBeTextScalar} = ""
                options.Dependencies string = string.empty(0,1)
                options.Settings = struct()
                options.Build (1,1) function_handle
                options.Force (1,1) logical = false
            end

            movieName = TutorialBuildRuntime.normalizeGeneratedAssetFileName(string(name), "movie");
            if any(strcmp(string({self.movieRecords.Name}), movieName))
                error("TutorialBuildRuntime:DuplicateMovieName", ...
                    "Duplicate tutorial movie name '%s'.", movieName);
            end

            posterName = "";
            posterPath = "";
            posterRelativePath = "";
            if string(options.Poster) ~= ""
                posterName = TutorialBuildRuntime.normalizeGeneratedAssetFileName(string(options.Poster), "poster");
                posterPath = fullfile(self.assetFolder, posterName);
                posterRelativePath = self.assetPagePrefix + posterName;
            end

            TutorialBuildRuntime.ensureFolderExists(self.assetFolder);

            targetPath = fullfile(self.assetFolder, movieName);
            stampPath = targetPath + ".stamp.json";
            stampText = self.movieStampText(movieName, posterName, options.Dependencies, options.Settings);

            didReuse = false;
            if ~options.Force
                didReuse = self.tryReuseMovieFromPreviousBuild(targetPath, posterPath, stampPath, movieName, posterName, stampText);
            end

            if ~didReuse
                TutorialBuildRuntime.ensureParentFolderExists(targetPath);
                if posterPath ~= ""
                    TutorialBuildRuntime.ensureParentFolderExists(posterPath);
                end

                options.Build(targetPath, posterPath);

                if ~isfile(targetPath)
                    error("TutorialBuildRuntime:MissingMovieOutput", ...
                        "tutorialMovieCapture did not create the expected movie file '%s'.", targetPath);
                end
                if posterPath ~= "" && ~isfile(posterPath)
                    error("TutorialBuildRuntime:MissingPosterOutput", ...
                        "tutorialMovieCapture did not create the expected poster file '%s'.", posterPath);
                end

                TutorialBuildRuntime.writeTextFile(stampPath, stampText);
            end

            self.movieRecords(end+1) = struct( ...
                "Name", movieName, ...
                "Caption", string(options.Caption), ...
                "RelativePath", self.assetPagePrefix + movieName, ...
                "PosterRelativePath", posterRelativePath);
        end

        function captureOutput(self, source, options)
            arguments
                self (1,1) TutorialBuildRuntime
                source
                options.Caption string = ""
                options.Language {mustBeTextScalar} = "text"
            end

            if isa(source, "function_handle")
                outputText = string(evalc("source()"));
            elseif isstring(source) || ischar(source)
                mustBeTextScalar(source)
                outputText = string(source);
            else
                error("TutorialBuildRuntime:InvalidOutputSource", ...
                    "tutorialOutputCapture source must be a function handle or a text scalar.");
            end

            self.outputRecords(end+1) = struct( ...
                "Caption", string(options.Caption), ...
                "Language", string(options.Language), ...
                "Text", TutorialBuildRuntime.normalizeCapturedOutput(outputText));
        end

        function figureRecords = getFigureRecords(self)
            figureRecords = self.figureRecords;
        end

        function movieRecords = getMovieRecords(self)
            movieRecords = self.movieRecords;
        end

        function outputRecords = getOutputRecords(self)
            outputRecords = self.outputRecords;
        end
    end

    methods (Access = private)
        function didReuse = tryReuseMovieFromPreviousBuild(self, targetPath, posterPath, stampPath, movieName, posterName, stampText)
            didReuse = false;
            if self.comparisonAssetFolder == ""
                return;
            end

            previousMoviePath = fullfile(self.comparisonAssetFolder, movieName);
            previousStampPath = previousMoviePath + ".stamp.json";
            if ~isfile(previousMoviePath) || ~isfile(previousStampPath)
                return;
            end

            previousStampText = string(fileread(previousStampPath));
            if previousStampText ~= stampText
                return;
            end

            previousPosterPath = "";
            if posterName ~= ""
                previousPosterPath = fullfile(self.comparisonAssetFolder, posterName);
                if ~isfile(previousPosterPath)
                    return;
                end
            end

            TutorialBuildRuntime.ensureParentFolderExists(targetPath);
            copyfile(previousMoviePath, targetPath, "f");
            if posterPath ~= ""
                TutorialBuildRuntime.ensureParentFolderExists(posterPath);
                copyfile(previousPosterPath, posterPath, "f");
            end
            TutorialBuildRuntime.writeTextFile(stampPath, stampText);
            didReuse = true;
        end

        function stampText = movieStampText(self, movieName, posterName, dependencies, settings)
            dependencyRecords = self.dependencyStampRecords(dependencies);
            stamp = struct( ...
                "Version", 1, ...
                "Movie", movieName, ...
                "Poster", posterName, ...
                "Dependencies", dependencyRecords, ...
                "Settings", TutorialBuildRuntime.canonicalizeForJSON(settings));
            stampText = string(jsonencode(stamp));
        end

        function dependencyRecords = dependencyStampRecords(self, dependencies)
            dependencyPaths = reshape(string(dependencies), [], 1);
            if self.sourcePath ~= ""
                dependencyPaths = [self.sourcePath; dependencyPaths];
            end
            dependencyPaths = dependencyPaths(strlength(strtrim(dependencyPaths)) > 0);
            dependencyPaths = unique(dependencyPaths, "stable");

            dependencyRecords = repmat(struct("Label", "", "Hash", ""), numel(dependencyPaths), 1);
            for iPath = 1:numel(dependencyPaths)
                dependencyPath = TutorialBuildRuntime.resolveDependencyPath(dependencyPaths(iPath));
                if ~isfile(dependencyPath)
                    error("TutorialBuildRuntime:MissingMovieDependency", ...
                        "Movie dependency '%s' does not exist.", dependencyPaths(iPath));
                end

                dependencyRecords(iPath) = struct( ...
                    "Label", TutorialBuildRuntime.dependencyLabelForStamp(dependencyPaths(iPath)), ...
                    "Hash", TutorialBuildRuntime.hashFile(dependencyPath));
            end
        end
    end

    methods (Static, Access = private)
        function tf = imagesMatch(firstPath, secondPath)
            [firstImage, ~, firstAlpha] = imread(firstPath);
            [secondImage, ~, secondAlpha] = imread(secondPath);

            tf = isequal(size(firstImage), size(secondImage)) ...
                && strcmp(class(firstImage), class(secondImage)) ...
                && isequal(firstImage, secondImage) ...
                && isequal(firstAlpha, secondAlpha);
        end

        function slug = slugify(textValue)
            slug = lower(strtrim(string(textValue)));
            slug = regexprep(slug, "[^a-z0-9]+", "-");
            slug = regexprep(slug, "^-+|-+$", "");
            if slug == ""
                error("TutorialBuildRuntime:InvalidFigureName", ...
                    "Tutorial figure names must contain letters or numbers.");
            end
        end

        function fileName = normalizeGeneratedAssetFileName(textValue, assetKind)
            fileName = strtrim(string(textValue));
            if fileName == ""
                error("TutorialBuildRuntime:InvalidGeneratedAssetName", ...
                    "Tutorial %s names must not be empty.", assetKind);
            end

            normalizedPath = replace(fileName, "\", "/");
            if contains(normalizedPath, "/")
                error("TutorialBuildRuntime:NestedGeneratedAssetPath", ...
                    "Tutorial %s '%s' must be a file name in the tutorial asset root, not a nested path.", ...
                    assetKind, fileName);
            end
            if contains(normalizedPath, "..")
                error("TutorialBuildRuntime:InvalidGeneratedAssetPath", ...
                    "Tutorial %s '%s' must not contain parent-directory segments.", assetKind, fileName);
            end

            [~, ~, extension] = fileparts(char(fileName));
            if strlength(string(extension)) == 0
                error("TutorialBuildRuntime:MissingGeneratedAssetExtension", ...
                    "Tutorial %s '%s' must include a file extension.", assetKind, fileName);
            end
        end

        function ensureFolderExists(folderPath)
            if ~isfolder(folderPath)
                mkdir(folderPath);
            end
        end

        function ensureParentFolderExists(filePath)
            parentFolder = fileparts(filePath);
            if strlength(parentFolder) == 0 || isfolder(parentFolder)
                return;
            end

            mkdir(parentFolder);
        end

        function writeTextFile(filePath, fileText)
            fileID = fopen(filePath, "w");
            assert(fileID ~= -1, "Could not open file for writing: %s", filePath);
            fwrite(fileID, char(fileText));
            fclose(fileID);
        end

        function dependencyPath = resolveDependencyPath(dependencyPath)
            dependencyPath = string(dependencyPath);
            if isfile(dependencyPath)
                dependencyPath = TutorialBuildRuntime.canonicalPath(dependencyPath);
            end
        end

        function label = dependencyLabelForStamp(pathValue)
            normalizedPath = replace(strtrim(string(pathValue)), "\", "/");
            pathSegments = split(normalizedPath, "/");
            pathSegments = pathSegments(strlength(pathSegments) > 0);
            if numel(pathSegments) > 4
                pathSegments = pathSegments(end-3:end);
            end
            label = join(pathSegments, "/");
        end

        function hashValue = hashFile(filePath)
            fileID = fopen(filePath, "r");
            assert(fileID ~= -1, "Could not open dependency file for reading: %s", filePath);
            fileBytes = fread(fileID, inf, "*uint8");
            fclose(fileID);
            hashValue = TutorialBuildRuntime.hashBytes(fileBytes);
        end

        function hashValue = hashBytes(fileBytes)
            messageDigest = java.security.MessageDigest.getInstance("SHA-256");
            messageDigest.update(int8(fileBytes));
            digestBytes = typecast(messageDigest.digest(), "uint8");
            hexMatrix = dec2hex(digestBytes, 2).';
            hashValue = lower(string(reshape(hexMatrix, 1, [])));
        end

        function canonicalValue = canonicalizeForJSON(value)
            if isstruct(value)
                fieldNames = sort(fieldnames(value));
                canonicalValue = value;
                for iElement = 1:numel(value)
                    elementStruct = struct();
                    for iField = 1:numel(fieldNames)
                        fieldName = fieldNames{iField};
                        elementStruct.(fieldName) = TutorialBuildRuntime.canonicalizeForJSON(value(iElement).(fieldName));
                    end
                    canonicalValue(iElement) = elementStruct;
                end
                return;
            end

            if iscell(value)
                canonicalValue = cell(size(value));
                for iValue = 1:numel(value)
                    canonicalValue{iValue} = TutorialBuildRuntime.canonicalizeForJSON(value{iValue});
                end
                return;
            end

            canonicalValue = value;
        end

        function pathValue = canonicalPath(pathValue)
            pathValue = string(java.io.File(char(pathValue)).getCanonicalPath());
        end

        function outputText = normalizeCapturedOutput(outputText)
            outputText = replace(string(outputText), sprintf("\r\n"), newline);
            outputText = replace(outputText, sprintf("\r"), newline);
            outputText = string(TutorialBuildRuntime.applyBackspaces(char(outputText)));
            outputText = regexprep(char(outputText), "</?strong>", "");
            outputLines = splitlines(string(outputText));
            outputLines = TutorialBuildRuntime.trimBlankEdges(outputLines);
            outputText = join(outputLines, newline);
        end

        function textValue = applyBackspaces(textValue)
            characters = char(textValue);
            cleanedCharacters = strings(0, 1);
            for iCharacter = 1:numel(characters)
                currentCharacter = characters(iCharacter);
                if currentCharacter == char(8)
                    if ~isempty(cleanedCharacters)
                        cleanedCharacters(end) = [];
                    end
                else
                    cleanedCharacters(end+1, 1) = string(currentCharacter); %#ok<AGROW>
                end
            end
            textValue = char(join(cleanedCharacters, ""));
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
    end
end
