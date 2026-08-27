classdef WVTestUnregisteredTracer < WVTracer
    methods
        function self = WVTestUnregisteredTracer(model,options)
            arguments
                model WVModel
                options.name {mustBeText} = "unregistered tracer"
                options.phi double
            end
            self@WVTracer(model,name=options.name,phi=options.phi);
        end
    end
end
