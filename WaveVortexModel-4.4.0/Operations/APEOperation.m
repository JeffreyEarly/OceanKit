classdef APEOperation < WVOperation
    % Compute APE with the same profile and material height as eta_true.
    methods
        function self = APEOperation(wvt)
            arguments
                wvt WVTransform
            end
            outputVariables = WVVariableAnnotation('ape',{'x','y','z'},'m2 s-2', 'available potential energy density');
            self@WVOperation('ape',outputVariables,@disp);
        end

        function varargout = compute(self,wvt,varargin)
            profile = EtaTrueOperation.profileForDiagnostics(wvt);
            materialHeight = wvt.Z-wvt.eta_true;
            varargout = {profile.availablePotentialEnergy(wvt.Z,materialHeight,wvt.g,wvt.rho0)};
        end
    end
end
