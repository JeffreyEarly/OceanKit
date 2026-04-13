classdef GriddedStreamfunctionBootstrapFitMetadata < CAAnnotatedClass
    % Canonical restart metadata for one bootstrap replicate fit.
    %
    % This helper persists the resolved fast and mesoscale knot vectors
    % needed to reconstruct one bootstrap replicate exactly.
    %
    % - Topic: Developer helpers
    % - Declaration: classdef GriddedStreamfunctionBootstrapFitMetadata < CAAnnotatedClass
    % - Developer: true

    properties (SetAccess = private)
        fastKnotPoints (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
        psiQKnotPoints (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
        psiRKnotPoints (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
        psiTKnotPoints (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
    end

    properties (Dependent)
        psiKnotPoints
    end

    methods
        function self = GriddedStreamfunctionBootstrapFitMetadata(options)
            arguments
                options.fastKnotPoints (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.psiQKnotPoints (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.psiRKnotPoints (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
                options.psiTKnotPoints (:,1) double {mustBeReal,mustBeFinite} = zeros(0, 1)
            end

            self@CAAnnotatedClass();

            self.fastKnotPoints = options.fastKnotPoints;
            self.psiQKnotPoints = options.psiQKnotPoints;
            self.psiRKnotPoints = options.psiRKnotPoints;
            self.psiTKnotPoints = options.psiTKnotPoints;
        end

        function psiKnotPoints = get.psiKnotPoints(self)
            psiKnotPoints = {self.psiQKnotPoints, self.psiRKnotPoints, self.psiTKnotPoints};
        end
    end

    methods (Static)
        function self = fromFit(fit)
            arguments
                fit (1,1) GriddedStreamfunction
            end

            self = GriddedStreamfunctionBootstrapFitMetadata( ...
                fastKnotPoints=fit.fastKnotPoints, ...
                psiQKnotPoints=fit.psiKnotPoints{1}, ...
                psiRKnotPoints=fit.psiKnotPoints{2}, ...
                psiTKnotPoints=fit.psiKnotPoints{3});
        end
    end

    methods (Static, Hidden)
        function self = annotatedClassFromGroup(group)
            vars = CAAnnotatedClass.propertyValuesFromGroup(group, GriddedStreamfunctionBootstrapFitMetadata.classRequiredPropertyNames());
            self = GriddedStreamfunctionBootstrapFitMetadata();
            self.fastKnotPoints = vars.fastKnotPoints;
            self.psiQKnotPoints = vars.psiQKnotPoints;
            self.psiRKnotPoints = vars.psiRKnotPoints;
            self.psiTKnotPoints = vars.psiTKnotPoints;
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0, 0);
            propertyAnnotations(end+1) = CADimensionProperty('fastKnotPoints', '', 'Fast temporal knot vector.');
            propertyAnnotations(end+1) = CADimensionProperty('psiQKnotPoints', '', 'Mesoscale q-knot vector.');
            propertyAnnotations(end+1) = CADimensionProperty('psiRKnotPoints', '', 'Mesoscale r-knot vector.');
            propertyAnnotations(end+1) = CADimensionProperty('psiTKnotPoints', '', 'Mesoscale t-knot vector.');
        end

        function names = classRequiredPropertyNames()
            names = {'fastKnotPoints', 'psiQKnotPoints', 'psiRKnotPoints', 'psiTKnotPoints'};
        end
    end
end
