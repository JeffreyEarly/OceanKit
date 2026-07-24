classdef WVAlongTrackObservingSystem < WVObservingSystem
    % An model observing system for along track data
    %
    % The observing system is automatically added when a
    % WVModelOutputGroupAlongTrack is added 
    %
    % - Declaration: classdef WVAlongTrackObservingSystem < WVObservingSystem
    properties (WeakHandle)
        alongTrackGroup WVModelOutputGroupAlongTrack
    end

    methods
        function self = WVAlongTrackObservingSystem(model,alongTrackGroup)
            %create a new along track observing system
            %
            % This class is initialized when a
            % `WVModelOutputGroupAlongTrack` is initialized, and thus does
            % not need to be initialized directly.
            %
            % - Topic: Initialization
            % - Declaration: self = WVAlongTrackObservingSystem(model,alongTrackGroup)
            % - Parameter model: the WVModel instance
            % - Parameter alongTrackGroup: name of the observing system
            % - Returns self: a new instance of WVAlongTrackObservingSystem
            arguments
                model WVModel
                alongTrackGroup
            end
            self@WVObservingSystem(model,alongTrackGroup.missionName);
            self.alongTrackGroup = alongTrackGroup;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function initializeStorage(self,group)
            spatialDimensionNames = self.model.wvt.spatialDimensionNames;
            for iVar=1:2
                attributes = containers.Map(KeyType='char',ValueType='any');
                attributes('units') = self.model.wvt.dimensionAnnotationWithName(spatialDimensionNames{iVar}).units;
                attributes('long_name') = strcat(self.model.wvt.dimensionAnnotationWithName(spatialDimensionNames{iVar}).description,' position of observation');
                group.addVariable(strcat("track",'_',spatialDimensionNames{iVar}),{"t"},type="double",attributes=attributes,isComplex=false);
            end

            for iField = 1:length(self.alongTrackGroup.fieldNames)
                fieldName = self.alongTrackGroup.fieldNames(iField);
                varAnnotation = self.model.wvt.propertyAnnotationWithName(fieldName);
                attributes = containers.Map(KeyType='char',ValueType='any');
                annotationAttributeNames = varAnnotation.attributes.keys;
                for iAttribute = 1:length(annotationAttributeNames)
                    attributes(annotationAttributeNames{iAttribute}) = varAnnotation.attributes(annotationAttributeNames{iAttribute});
                end
                attributes('units') = varAnnotation.units;
                attributes('long_name') = string(varAnnotation.description) + ", sampled along the " + string(self.alongTrackGroup.missionName) + " mission ground track.";
                group.addVariable(fieldName,{'t'},type="double",attributes=attributes,isComplex=false);
            end
        end

        function writeTimeStepToFile(self,group,outputIndices)
            if ~isinf(self.alongTrackGroup.repeatCycle)
                iPassover = find( abs(self.alongTrackGroup.firstPassoverTime - mod(self.wvt.t,self.alongTrackGroup.repeatCycle)) < 1, 1, 'first');
            else
                iPassover = find( abs(self.alongTrackGroup.firstPassoverTime - self.wvt.t) < 1, 1, 'first');
            end
            tracks = self.alongTrackGroup.tracks;
            group.variableWithName("track_x").setValueAlongDimensionAtIndex(tracks{iPassover}.x,'t',outputIndices);
            group.variableWithName("track_y").setValueAlongDimensionAtIndex(tracks{iPassover}.y,'t',outputIndices);
            fieldNames = cellstr(self.alongTrackGroup.fieldNames);
            fieldValues = cell(size(fieldNames));
            [fieldValues{:}] = self.model.wvt.variableAtPositionWithName(tracks{iPassover}.x,tracks{iPassover}.y,[],fieldNames{:});
            for iField = 1:length(fieldNames)
                values = reshape(fieldValues{iField},[],1);
                group.variableWithName(fieldNames{iField}).setValueAlongDimensionAtIndex(values,'t',outputIndices);
            end
        end

        function os = observingSystemWithResolutionOfTransform(self,wvtX2)
            %create a new WVObservingSystem with a new resolution
            %
            % Subclasses to should override this method an implement the
            % correct logic.
            %
            % - Topic: Initialization
            % - Declaration: os = observingSystemWithResolutionOfTransform(self,wvtX2)
            % - Parameter wvtX2: the WVTransform with increased resolution
            % - Returns force: a new instance of WVObservingSystem
            os = WVAlongTrackObservingSystem(wvtX2,self.name);
            error('this needs to be implemented');
        end
    end
end
