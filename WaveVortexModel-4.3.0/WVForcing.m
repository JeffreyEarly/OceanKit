classdef WVForcing < handle & matlab.mixin.Heterogeneous & CAAnnotatedClass
    % Add forcing or dissipation to a wave-vortex transform.
    %
    % `WVForcing` is the base class for forcing and closure objects
    % attached to a `WVTransform`. Use one of the supplied subclasses or
    % subclass this interface to implement a custom forcing. Each instance
    % belongs to the transform supplied at construction and is registered by
    % its unique `name` through `WVTransform.addForcing`.
    %
    % Forcing is applied in three stages. Physical-space forcing contributes
    % tendencies to $$(u,v,\eta)$$ for hydrostatic flow or
    % $$(u,v,w,\eta)$$ for nonhydrostatic flow. Those tendencies are projected
    % into wave-vortex space before spectral forcing contributes directly to
    % $$(F_+,F_-,F_0)$$. Spectral-amplitude forcing may then update `Ap`, `Am`,
    % and `A0` directly. Potential-vorticity variants contribute only to the
    % zero-frequency tendency or amplitude.
    %
    % A custom subclass declares one or more stages with `forcingType` and
    % overrides the corresponding evaluation methods. Spatial forcing is
    % evaluated before projection, spectral forcing is evaluated after
    % projection, and spectral-amplitude forcing is evaluated last. Within
    % each stage, smaller `priority` values are evaluated first.
    %
    % `WVTransform.addForcing` registers forcing objects, orders compatible
    % contributions by priority, and exposes their energy transfers through
    % the transform diagnostics.
    %
    % - Topic: Create the forcing
    % - Topic: Inspect forcing configuration
    % - Topic: Configure the forcing
    % - Topic: Inspect forcing or damping scales
    % - Topic: Evaluate prescribed forcing
    % - Topic: Generate forcing inputs
    % - Topic: Implement forcing evaluation
    % - Topic: Convert forcing resolution
    % - Topic: Forcing persistence
    % - Topic: Forcing internals
    % - Declaration: classdef WVForcing < handle

    % Can't set wvt here because WeakHandle cannot use an abstract class
    % properties (WeakHandle, GetAccess=public, SetAccess=protected)
    % 
    % end

    properties (GetAccess=public, SetAccess=protected)
        % Transform to which this forcing belongs.
        %
        % A forcing can be registered only with the same transform instance
        % supplied to its constructor.
        %
        % - Topic: Inspect forcing configuration
        wvt

        % Name used to register this forcing with its transform.
        %
        % Names identify forcing objects in transform lookup, replacement,
        % removal, diagnostics, and persistence operations.
        %
        % - Topic: Inspect forcing configuration
        name

        % Evaluation stages implemented by this forcing.
        %
        % Each `WVForcingType` value declares the method or methods that a
        % subclass implements:
        %
        % | Type | Evaluation method |
        % | --- | --- |
        % | `HydrostaticSpatial` | `addHydrostaticSpatialForcing` |
        % | `NonhydrostaticSpatial` | `addNonhydrostaticSpatialForcing` |
        % | `PVSpatial` | `addPotentialVorticitySpatialForcing` |
        % | `Spectral` | `addSpectralForcing` |
        % | `PVSpectral` | `addPotentialVorticitySpectralForcing` |
        % | `SpectralAmplitude` | `setSpectralForcing` and `setSpectralAmplitude` |
        % | `PVSpectralAmplitude` | `setPotentialVorticitySpectralForcing` and `setPotentialVorticitySpectralAmplitude` |
        %
        % Spectral-amplitude forcing first modifies the coefficient tendency
        % with the corresponding `*SpectralForcing` method. After an
        % integration step, the corresponding `*SpectralAmplitude` method
        % restores the constrained coefficient values exactly.
        %
        % - Topic: Inspect forcing configuration
        forcingType WVForcingType = WVForcingType.empty(0,0)

        % Whether this forcing is a small-scale closure.
        %
        % Closure objects remove variance at unresolved scales and cause
        % `WVTransform.hasClosure` to return `true`.
        %
        % - Topic: Inspect forcing configuration
        isClosure logical =  false

        % Order within a forcing stage, from 0 first to 255 last.
        %
        % The default is 255. Priority is compared only among forcing objects
        % in the same evaluation stage: all spatial forcing is evaluated
        % before spectral forcing, regardless of priority. Nonlinear advection
        % and explicit antialiasing use priority 127 so they precede ordinary
        % default-priority forcing in their respective stages.
        %
        % - Topic: Inspect forcing configuration
        priority uint8 = 255
    end

    methods

        function self = WVForcing(wvt,name,forcingType)
            % Initialize the base state for a forcing subclass.
            %
            % Subclass constructors call this constructor with their owning
            % transform, unique registry name, and implemented evaluation
            % stages. Users normally construct a concrete supplied subclass.
            %
            % - Topic: Create the forcing
            % - Declaration: self = WVForcing(wvt,name,forcingType)
            % - Parameter wvt: transform that owns and evaluates the forcing
            % - Parameter name: unique forcing registry name
            % - Parameter forcingType: one or more `WVForcingType` evaluation stages implemented by the subclass
            % - Returns self: initialized `WVForcing` base instance
            arguments
                wvt WVTransform
                name {mustBeText}
                forcingType WVForcingType {mustBeNonempty}
            end
            self@CAAnnotatedClass(inheritedClassName=class(wvt));
            self.wvt = wvt;
            self.name = name;
            self.forcingType = forcingType;
        end

        function contract = portableImplementationContract(self)
            % Describe availability of the paired portable C++ implementation.
            %
            % The base class is intentionally unavailable. A supported
            % concrete forcing overrides this method with a versioned,
            % data-only contract. Portable execution still requires an exact
            % source-level C++ registration.
            %
            % - Topic: Forcing internals
            % - Declaration: contract = portableImplementationContract(self)
            % - Parameter self: forcing to inspect
            % - Returns contract: scalar portable implementation contract
            % - Developer: true
            contract = WVInternal.portableImplementationContract(string(class(self)),string(class(self)),"unavailable","No matching portable C++ forcing implementation is declared.",struct());
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
            % Add hydrostatic physical-space tendencies.
            %
            % Subclasses declaring `HydrostaticSpatial` override this hook.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: [Fu,Fv,Feta] = addHydrostaticSpatialForcing(wvt,Fu,Fv,Feta)
            % - Parameter wvt: transform evaluating the forcing
            % - Parameter Fu: accumulated zonal-velocity tendency
            % - Parameter Fv: accumulated meridional-velocity tendency
            % - Parameter Feta: accumulated isopycnal-displacement tendency
            % - Returns Fu: updated zonal-velocity tendency
            % - Returns Fv: updated meridional-velocity tendency
            % - Returns Feta: updated isopycnal-displacement tendency
            % - Developer: true
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
            % Add nonhydrostatic physical-space tendencies.
            %
            % Subclasses declaring `NonhydrostaticSpatial` override this hook.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: [Fu,Fv,Fw,Feta] = addNonhydrostaticSpatialForcing(wvt,Fu,Fv,Fw,Feta)
            % - Parameter wvt: transform evaluating the forcing
            % - Parameter Fu: accumulated zonal-velocity tendency
            % - Parameter Fv: accumulated meridional-velocity tendency
            % - Parameter Fw: accumulated vertical-velocity tendency
            % - Parameter Feta: accumulated isopycnal-displacement tendency
            % - Returns Fu: updated zonal-velocity tendency
            % - Returns Fv: updated meridional-velocity tendency
            % - Returns Fw: updated vertical-velocity tendency
            % - Returns Feta: updated isopycnal-displacement tendency
            % - Developer: true
        end

        function [Fp, Fm, F0] = addSpectralForcing(self, wvt, Fp, Fm, F0)
            % Add wave-vortex coefficient tendencies in spectral space.
            %
            % Subclasses declaring `Spectral` override this hook.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: [Fp,Fm,F0] = addSpectralForcing(wvt,Fp,Fm,F0)
            % - Parameter wvt: transform evaluating the forcing
            % - Parameter Fp: accumulated `Ap` tendency
            % - Parameter Fm: accumulated `Am` tendency
            % - Parameter F0: accumulated `A0` tendency
            % - Returns Fp: updated `Ap` tendency
            % - Returns Fm: updated `Am` tendency
            % - Returns F0: updated `A0` tendency
            % - Developer: true
        end

        function [Ap, Am, A0] = setSpectralAmplitude(self, wvt, Ap, Am, A0)
            % Restore selected wave-vortex coefficients after a model step.
            %
            % Subclasses declaring `SpectralAmplitude` override this hook to
            % enforce their constrained coefficient values exactly.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: [Ap,Am,A0] = setSpectralAmplitude(wvt,Ap,Am,A0)
            % - Parameter wvt: transform evaluating the forcing
            % - Parameter Ap: positive-frequency coefficients
            % - Parameter Am: negative-frequency coefficients
            % - Parameter A0: zero-frequency coefficients
            % - Returns Ap: updated positive-frequency coefficients
            % - Returns Am: updated negative-frequency coefficients
            % - Returns A0: updated zero-frequency coefficients
            % - Developer: true
        end

        function [Fp, Fm, F0] = setSpectralForcing(self, wvt, Fp, Fm, F0)
            % Modify tendencies for a spectral-amplitude constraint.
            %
            % Subclasses declaring `SpectralAmplitude` override this hook to
            % cancel or replace the tendency of constrained coefficients.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: [Fp,Fm,F0] = setSpectralForcing(wvt,Fp,Fm,F0)
            % - Parameter wvt: transform evaluating the forcing
            % - Parameter Fp: accumulated `Ap` tendency
            % - Parameter Fm: accumulated `Am` tendency
            % - Parameter F0: accumulated `A0` tendency
            % - Returns Fp: updated `Ap` tendency
            % - Returns Fm: updated `Am` tendency
            % - Returns F0: updated `A0` tendency
            % - Developer: true
        end



        function F0 = addPotentialVorticitySpatialForcing(self, wvt, F0)
            % Add a physical-space QGPV tendency.
            %
            % Subclasses declaring `PVSpatial` override this hook.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: F0 = addPotentialVorticitySpatialForcing(wvt,F0)
            % - Parameter wvt: QG transform evaluating the forcing
            % - Parameter F0: accumulated physical-space QGPV tendency
            % - Returns F0: updated physical-space QGPV tendency
            % - Developer: true
        end

        function F0 = addPotentialVorticitySpectralForcing(self, wvt, F0)
            % Add a spectral QGPV tendency.
            %
            % Subclasses declaring `PVSpectral` override this hook.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: F0 = addPotentialVorticitySpectralForcing(wvt,F0)
            % - Parameter wvt: QG transform evaluating the forcing
            % - Parameter F0: accumulated spectral QGPV tendency
            % - Returns F0: updated spectral QGPV tendency
            % - Developer: true
        end

        function A0 = setPotentialVorticitySpectralAmplitude(self, wvt, A0)
            % Restore selected QG coefficients after a model step.
            %
            % Subclasses declaring `PVSpectralAmplitude` override this hook.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: A0 = setPotentialVorticitySpectralAmplitude(wvt,A0)
            % - Parameter wvt: QG transform evaluating the forcing
            % - Parameter A0: zero-frequency coefficients
            % - Returns A0: updated zero-frequency coefficients
            % - Developer: true
        end

        function F0 = setPotentialVorticitySpectralForcing(self, wvt, F0)
            % Modify QGPV tendencies for a spectral-amplitude constraint.
            %
            % Subclasses declaring `PVSpectralAmplitude` override this hook.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: F0 = setPotentialVorticitySpectralForcing(wvt,F0)
            % - Parameter wvt: QG transform evaluating the forcing
            % - Parameter F0: accumulated zero-frequency tendency
            % - Returns F0: updated zero-frequency tendency
            % - Developer: true
        end

        function didGetRemovedFromTransform(self, wvt)
            % Release resources when a forcing is removed from its transform.
            %
            % Subclasses with listeners or other transform-owned resources
            % override this lifecycle callback.
            %
            % - Topic: Forcing internals
            % - Declaration: didGetRemovedFromTransform(wvt)
            % - Parameter wvt: transform from which the forcing was removed
            % - Developer: true
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            % Rebuild a forcing for a compatible transform resolution.
            %
            % Subclasses that support transform resolution conversion override
            % this method and preserve their user configuration.
            %
            % - Topic: Convert forcing resolution
            % - Declaration: force = forcingWithResolutionOfTransform(wvtX2)
            % - Parameter wvtX2: compatible transform at the target resolution
            % - Returns force: equivalent forcing owned by `wvtX2`
            % - Developer: true
            force = WVForcing(wvtX2,self.name);
        end

    end

    methods (Access=protected)
        function contract = supportedPortableImplementationContract(self,typeIdentifier,payload)
            contract = WVInternal.portableImplementationContract(string(class(self)),string(typeIdentifier),"supported","",payload);
        end
    end

    methods (Sealed)
        function tf = eq(obj1, obj2)
            tf = eq@handle(obj1, obj2);
        end

        function tf = ne(obj1, obj2)
            tf = ne@handle(obj1, obj2);
        end
    end

    methods (Static)

        function forceTypes = spatialFluxTypes()
            % Return the physical-space forcing types.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: forceTypes = spatialFluxTypes()
            % - Returns forceTypes: `HydrostaticSpatial`, `NonhydrostaticSpatial`, and `PVSpatial`
            % - Developer: true
            forceTypes = WVForcingType(["HydrostaticSpatial","NonhydrostaticSpatial","PVSpatial"]);
        end

        function forceTypes = spectralFluxTypes()
            % Return the spectral-tendency forcing types.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: forceTypes = spectralFluxTypes()
            % - Returns forceTypes: `Spectral` and `PVSpectral`
            % - Developer: true
            forceTypes = WVForcingType(["Spectral","PVSpectral"]);
        end

        function forceTypes = spectralAmplitudeTypes()
            % Return the spectral-amplitude forcing types.
            %
            % - Topic: Implement forcing evaluation
            % - Declaration: forceTypes = spectralAmplitudeTypes()
            % - Returns forceTypes: `SpectralAmplitude` and `PVSpectralAmplitude`
            % - Developer: true
            forceTypes = WVForcingType(["SpectralAmplitude","PVSpectralAmplitude"]);
        end

        function force = forcingFromGroup(group,wvt)
            % Restore a concrete forcing from a NetCDF group.
            %
            % The annotated class name and required properties select and
            % reconstruct the concrete forcing subclass.
            %
            % - Topic: Forcing persistence
            % - Declaration: force = forcingFromGroup(group,wvt)
            % - Parameter group: NetCDF group containing annotated forcing state
            % - Parameter wvt: transform that will own the restored forcing
            % - Returns force: restored concrete `WVForcing` instance
            % - Developer: true
            arguments
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            className = group.attributes('AnnotatedClass');
            vars = CAAnnotatedClass.requiredPropertiesFromGroup(group);
            if isempty(vars)
                force = feval(className,wvt);
            else
                options = namedargs2cell(vars);
                force = feval(className,wvt,options{:});
            end
        end
    end
end
