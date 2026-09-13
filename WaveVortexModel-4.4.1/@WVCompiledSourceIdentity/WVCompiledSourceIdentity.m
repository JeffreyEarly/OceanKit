classdef (Hidden, Sealed) WVCompiledSourceIdentity < handle
    % Opaque runtime identity for one transform's compiled scientific source.
    %
    % Identity follows MATLAB handle equality and survives clear functions.
    % Transforms retain this token in transient storage. It counts native
    % attachments without retaining the transform or entering persistence.
    %
    % - Developer: true
    % - Topic: Compiled transform internals
    % - Declaration: classdef (Hidden, Sealed) WVCompiledSourceIdentity < handle

    properties (Access=private)
        nativeLeaseCount (1,1) uint64 = uint64(0)
    end

    methods (Hidden)
        function acquireLease(self)
            if self.nativeLeaseCount == intmax("uint64")
                error("WaveVortexModel:CompiledTransformIdentityLeaseOverflow","The compiled source identity has too many active native leases.")
            end
            self.nativeLeaseCount = self.nativeLeaseCount + uint64(1);
        end

        function releaseLease(self)
            if self.nativeLeaseCount == 0
                error("WaveVortexModel:CompiledTransformIdentityLeaseUnderflow","The compiled source identity has no active native lease to release.")
            end
            self.nativeLeaseCount = self.nativeLeaseCount - uint64(1);
        end

        function flag = hasLease(self)
            flag = self.nativeLeaseCount ~= 0;
        end
    end
end
