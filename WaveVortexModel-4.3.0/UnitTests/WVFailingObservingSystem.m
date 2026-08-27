classdef WVFailingObservingSystem < WVObservingSystem
    properties (SetAccess=private)
        failurePhase string
    end

    methods
        function self = WVFailingObservingSystem(model,options)
            arguments
                model WVModel
                options.name (1,1) string = "failing observer"
                options.failurePhase (1,1) string {mustBeMember(options.failurePhase,["initialize","write"])} = "initialize"
            end
            self@WVObservingSystem(model,options.name);
            self.failurePhase = options.failurePhase;
        end

        function initializeStorage(self,group)
            if self.failurePhase == "initialize"
                error('The test observing system failed during storage initialization.');
            end
            group.addVariable(self.name,{'t'},type="double",isComplex=false);
        end

        function writeTimeStepToFile(self,~,~)
            if self.failurePhase == "write"
                error('The test observing system failed while writing a time step.');
            end
        end
    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {'name','failurePhase'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CAPropertyAnnotation('name','name of the deliberately failing observing system');
            propertyAnnotations(end+1) = CAPropertyAnnotation('failurePhase','output phase in which the test observer throws');
        end
    end
end
