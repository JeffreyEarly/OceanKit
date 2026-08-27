classdef WVTestPortablePointDiagnostic < WVObservingSystem
    properties (SetAccess=immutable)
        x (:,1) double
        y (:,1) double
        z (:,1) double
        fieldName (1,1) string
        scale (1,1) double
        offset (1,1) double
        interpolation (1,1) string
    end

    methods
        function self = WVTestPortablePointDiagnostic(model,options)
            arguments
                model WVModel
                options.name {mustBeText} = "point diagnostic"
                options.x (:,1) double
                options.y (:,1) double
                options.z (:,1) double
                options.fieldName (1,1) string
                options.scale (1,1) double = 1
                options.offset (1,1) double = 0
                options.interpolation (1,1) string {mustBeMember(options.interpolation,["linear" "spline"])} = "linear"
            end
            self@WVObservingSystem(model,options.name);
            if isempty(options.x) || ~isequal(size(options.x),size(options.y),size(options.z))
                error("WVTestPortablePointDiagnostic:InvalidCoordinates","x, y, and z must be equal nonempty column vectors.");
            end
            self.x = options.x;
            self.y = options.y;
            self.z = options.z;
            self.fieldName = options.fieldName;
            self.scale = options.scale;
            self.offset = options.offset;
            self.interpolation = options.interpolation;
        end

        function contract = portableImplementationContract(self)
            payload = struct("name",string(self.name),"fieldNames",self.fieldName,"x",self.x,"y",self.y,"z",self.z,"outputScale",self.scale,"outputOffset",self.offset,"trackedFieldInterpolation",self.interpolation);
            contract = self.supportedPortableImplementationContract("WVTestPortablePointDiagnostic",payload);
        end
    end
end
