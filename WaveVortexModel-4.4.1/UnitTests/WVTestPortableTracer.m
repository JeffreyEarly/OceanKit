classdef WVTestPortableTracer < WVTracer
    methods
        function self = WVTestPortableTracer(model,options)
            arguments
                model WVModel
                options.name {mustBeText} = "test portable tracer"
                options.isXYOnly (1,1) logical = false
                options.phi double
                options.absTolerance = 1e-5
                options.shouldAntialias (1,1) logical = true
            end
            self@WVTracer(model,name=options.name,isXYOnly=options.isXYOnly,phi=options.phi,absTolerance=options.absTolerance,shouldAntialias=options.shouldAntialias);
        end

        function contract = portableImplementationContract(self)
            payload = struct("name",string(self.name),"isXYOnly",logical(self.isXYOnly),"stateShape",uint64(size(self.phi)),"absTolerance",double(self.absTolerance),"shouldAntialias",logical(self.shouldAntialias));
            contract = self.supportedPortableImplementationContract("WVTestPortableTracer",payload);
        end
    end
end
