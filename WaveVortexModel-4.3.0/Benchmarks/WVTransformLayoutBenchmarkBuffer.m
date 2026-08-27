classdef WVTransformLayoutBenchmarkBuffer < handle
    % Hold one persistent full-complex benchmark buffer.
    properties
        value
        physicalSize
        storageShape
    end

    methods
        function self = WVTransformLayoutBenchmarkBuffer(sz,storageShape)
            arguments
                sz (1,3) double {mustBeInteger,mustBePositive}
                storageShape (1,1) string {mustBeMember(storageShape,["full-3d","rows-2d"])}
            end
            self.physicalSize = sz;
            self.storageShape = storageShape;
            if storageShape == "rows-2d"
                self.value = complex(zeros(prod(sz(1:2)),sz(3)));
            else
                self.value = complex(zeros(sz));
            end
        end

        function reset(self)
            self.value(:) = 0;
        end
    end
end
