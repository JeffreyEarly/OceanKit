classdef AnnotationTestClassB < handle & CAAnnotatedClass
    properties
        x
        y
    end

    methods
        function self = AnnotationTestClassB(options)
            arguments
                options.x
                options.y
            end
            requiredProperties = feval(strcat(class(self),'.classRequiredPropertyNames'));
            canInitializeDirectly = all(isfield(options,requiredProperties));
            if canInitializeDirectly
                for iVar = 1:length(requiredProperties)
                    name = requiredProperties{iVar};
                    self.(name) = options.(name);
                end
            else
                self.x = linspace(-100,0,101);
                self.y = self.x.^2;
            end
        end
    end

    methods (Static)

        function vars = classRequiredPropertyNames()
            vars = {'x','y'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CADimensionProperty('x', 'm', 'x coordinate');
            propertyAnnotations(end+1) = CANumericProperty('y', {'x'}, '', 'A variable quadratic in x.');
        end
        
        function atc = annotatedTestClassFromFile(path)
            atc = CAAnnotatedClass.annotatedClassFromFile(path);
        end

    end
end
