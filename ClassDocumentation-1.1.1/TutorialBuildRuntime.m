classdef TutorialBuildRuntime < handle
    properties
        assetFolder string
        assetPagePrefix string
        comparisonAssetFolder string = ""
        figureRecords struct = struct("Name", {}, "Caption", {}, "RelativePath", {})
    end

    methods
        function self = TutorialBuildRuntime(assetFolder, options)
            arguments
                assetFolder (1,1) string
                options.assetPagePrefix (1,1) string = "./"
                options.comparisonAssetFolder (1,1) string = ""
            end

            self.assetFolder = assetFolder;
            self.assetPagePrefix = options.assetPagePrefix;
            self.comparisonAssetFolder = options.comparisonAssetFolder;
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

            if ~isfolder(self.assetFolder)
                mkdir(self.assetFolder);
            end

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

        function figureRecords = getFigureRecords(self)
            figureRecords = self.figureRecords;
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
    end
end
