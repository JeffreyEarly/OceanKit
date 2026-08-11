classdef Topic < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name char
        methodAnnotations
        subtopics
    end

    methods
        function self = Topic(name)
            self.name = name;
            self.subtopics = Topic.empty(0,0);
            self.methodAnnotations = MethodDocumentation.empty(0,0);
        end

        function addMethod(self,methodAnnotation)
            arguments
                self Topic
                methodAnnotation {mustBeNonempty}
            end
            self.methodAnnotations(end+1) = methodAnnotation;
            [~,indices] = sort([self.methodAnnotations.nav_order]);
            self.methodAnnotations = self.methodAnnotations(indices);
        end

        function addSubtopic(self,subtopic)
            arguments
                self Topic
                subtopic Topic
            end
            self.subtopics(end+1) = subtopic;
        end

        function subtopic = subtopicWithName(self,subtopicName)
            arguments (Input)
                self Topic
                subtopicName char
            end
            subtopic = [];
            for iSubtopic=1:length(self.subtopics)
                if strcmpi(self.subtopics(iSubtopic).name,subtopicName)
                    subtopic = self.subtopics(iSubtopic);
                end
            end
        end

        function description(self,options)
            arguments
                self Topic
                options.indent char = '';
            end
            fprintf('%s- %s\n',options.indent,self.name);
            for iMethod=1:length(self.methodAnnotations)
                fprintf('%s  - %s\n',options.indent,self.methodAnnotations(iMethod).name);
            end
            for iSubtopic=1:length(self.subtopics)
                self.subtopics(iSubtopic).description(indent=[options.indent,'  ']);
            end
        end
    end

    methods (Static)
        function detailedDescription = trimTopicsFromString(detailedDescription)
            % Remove any topic metadata from a string
            arguments
                detailedDescription
            end
            % extract topics and the detailed description (minus those topics)
            subsubtopicExpression = '- topic:([ \t]*)(?<topicName>[^—\r\n]+)—([ \t]*)(?<subtopicName>[^\r\n]+)—([ \t]*)(?<subsubtopicName>[^\r\n]+)(?:$|\n)';
            detailedDescription = regexprep(detailedDescription,subsubtopicExpression,'','ignorecase');
            subtopicExpression = '- topic:([ \t]*)(?<topicName>[^—\r\n]+)—([ \t]*)(?<subtopicName>[^\r\n]+)(?:$|\n)';
            detailedDescription = regexprep(detailedDescription,subtopicExpression,'','ignorecase');
            topicExpression = '- topic:([ \t]*)(?<topicName>[^\r\n]+)(?:$|\n)';
            detailedDescription = regexprep(detailedDescription,topicExpression,'','ignorecase');
        end

        function rootTopic = topicsFromString(detailedDescription)
            % Extracts topics (and nested subtopics) from a detailedDescription and
            % creates a Topic tree useful for creating an ordered topic index.
            %
            % Each line of the form
            %   - Topic: A
            %   - Topic: A — B
            %   - Topic: A — B — C
            % is interpreted as a path A -> B -> C in the topic tree.
            %
            % The metadata lines themselves are stripped from detailedDescription
            % (locally in this function; if you need the cleaned text returned, we
            % can add that as a second output).

            arguments
                detailedDescription
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Step 1) Extract all "- Topic: ..." lines from the description
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Match whole lines that start with "- Topic:" (case-insensitive),
            % capturing everything after "Topic:" up to the end of the line.
            topicLineExpression = '^[ \t]*- topic:([^\r\n]+)\r?(?:\n|$)';

            % tokens: cell array where each cell has { <text after Topic:> }
            topicLines = regexp(detailedDescription, topicLineExpression, ...
                'tokens', 'ignorecase', 'lineanchors');

            % Remove those lines from the detailed description (locally)
            detailedDescription = regexprep(detailedDescription, topicLineExpression, ...
                '', 'ignorecase', 'lineanchors');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Step 2) Build the Topic tree from the parsed lines
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            rootTopic = Topic('Root');

            % Helper: get or create a subtopic under a given parent
            function child = getOrCreateSubtopic(parent, name)
                name = strtrim(name);
                if isempty(name)
                    child = parent;
                    return
                end
                child = parent.subtopicWithName(name);
                if isempty(child)
                    child = Topic(name);
                    parent.addSubtopic(child);
                end
            end

            % Each topic line may represent A, A—B, A—B—C, or deeper
            for iLine = 1:numel(topicLines)
                % topicLines{iLine} is a 1x1 cell, containing the captured string
                fullPathStr = topicLines{iLine}{1};

                % Split on em dash to get path components
                parts = regexp(fullPathStr, '—', 'split');

                % Walk down the tree, creating nodes as needed
                currentNode = rootTopic;
                for iPart = 1:numel(parts)
                    partName = strtrim(parts{iPart});
                    if isempty(partName)
                        continue
                    end
                    currentNode = getOrCreateSubtopic(currentNode, partName);
                end
            end
        end
    end
end