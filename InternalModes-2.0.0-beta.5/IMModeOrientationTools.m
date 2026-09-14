classdef (Hidden, Sealed) IMModeOrientationTools
    % Define the canonical phase convention for internal-mode families.

    properties (Constant)
        % Stable identifier for the canonical mode-orientation convention.
        convention = "shallowInteriorGPositive-v1"
    end

    methods (Static)
        function signs = shallowInteriorGPositive(options)
            % Return signs that make each `G` mode positive below the surface.
            %
            % A resolved nonzero surface value sets the sign directly. If
            % the surface value is zero, the one-sided Taylor coefficient
            % $$-D G_z(z_s)$$ sets the sign because the ocean interior is
            % at $$z_s-\delta$$. A known `F`-form zero-eigenvalue mode may
            % instead use `F` when its diagnostic `G` structure vanishes
            % identically.
            %
            % - Topic: Developer topics
            % - Declaration: signs = shallowInteriorGPositive(options)
            % - Parameter options.GValues: raw `G` values on a representative grid
            % - Parameter options.GzSurface: raw surface derivatives of `G`
            % - Parameter options.FValues: aligned raw `F` values on the same grid
            % - Parameter options.depth: positive domain depth
            % - Parameter options.surfaceIndex: surface row of the sampled arrays
            % - Parameter options.allowFFallback: modes allowed to use `F` when `G` vanishes
            % - Returns signs: row vector containing one `-1` or `+1` per mode
            % - Developer: true
            arguments
                options.GValues (:,:) double {mustBeReal, mustBeFinite}
                options.GzSurface (1,:) double {mustBeReal, mustBeFinite}
                options.FValues (:,:) double {mustBeReal, mustBeFinite}
                options.depth (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
                options.surfaceIndex (1,1) double {mustBeInteger, mustBePositive}
                options.allowFFallback (1,:) logical
            end

            nModes = size(options.GValues,2);
            if size(options.FValues,1) ~= size(options.GValues,1) || size(options.FValues,2) ~= nModes
                error("IMModeOrientationTools:InconsistentModeArrays", "FValues and GValues must have the same size.");
            end
            if length(options.GzSurface) ~= nModes || length(options.allowFFallback) ~= nModes
                error("IMModeOrientationTools:InconsistentModeCount", "GzSurface and allowFFallback must contain one value per mode column.");
            end
            if options.surfaceIndex > size(options.GValues,1)
                error("IMModeOrientationTools:InvalidSurfaceIndex", "surfaceIndex must select a row of GValues and FValues.");
            end

            relativeTolerance = 1e-10;
            GScale = max(abs(options.GValues),[],1);
            FScale = max(abs(options.FValues),[],1);
            GSurface = options.GValues(options.surfaceIndex,:);
            inwardDerivative = -options.depth*options.GzSurface;

            GIsAbsent = options.allowFFallback & GScale <= relativeTolerance.*FScale;
            surfaceIsResolved = ~GIsAbsent & abs(GSurface) > relativeTolerance.*GScale;
            derivativeIsResolved = ~GIsAbsent & ~surfaceIsResolved & abs(inwardDerivative) > relativeTolerance.*GScale;

            reference = nan(1,nModes);
            reference(surfaceIsResolved) = GSurface(surfaceIsResolved);
            reference(derivativeIsResolved) = inwardDerivative(derivativeIsResolved);

            fallbackModes = GIsAbsent;
            if any(fallbackModes)
                [~,maximumIndices] = max(abs(options.FValues(:,fallbackModes)),[],1);
                fallbackColumns = find(fallbackModes);
                linearIndices = sub2ind(size(options.FValues),maximumIndices,fallbackColumns);
                reference(fallbackModes) = options.FValues(linearIndices);
            end

            unresolved = ~isfinite(reference) | reference == 0;
            if any(unresolved)
                error("IMModeOrientationTools:IndeterminateOrientation", ...
                    "Could not determine the shallow-interior G orientation for mode column(s) %s. A nontrivial second-order G mode must have a resolved surface value or surface derivative.", ...
                    join(string(find(unresolved)),", "));
            end

            signs = ones(1,nModes);
            signs(reference < 0) = -1;
        end
    end
end
