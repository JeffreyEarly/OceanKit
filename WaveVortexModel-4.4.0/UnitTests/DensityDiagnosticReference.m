classdef DensityDiagnosticReference
    % Independent mathematical fixtures for density diagnostics, not a sorter.
    % Empirical distributions describe finite parcels with explicit volumes.
    % The twist describes a continuous map; its sampled grid is only a
    % quadrature approximation to the conserved continuous distribution.

    methods (Static)
        function distribution = empiricalDistribution(values,weights)
            arguments (Input)
                values double {mustBeReal,mustBeFinite}
                weights double {mustBeReal,mustBeFinite,mustBeNonnegative}
            end
            if numel(values) ~= numel(weights) || isempty(values) || sum(weights(:)) <= 0
                error('DensityDiagnosticReference:InvalidWeights','Supply one nonnegative volume per value and a positive total volume.');
            end
            retained = weights(:) > 0;
            values = values(:);
            weights = weights(:);
            [parcelValues,~,indices] = unique(values(retained));
            volumes = accumarray(indices,weights(retained));
            distribution = struct('values',parcelValues,'masses',volumes/sum(volumes), ...
                'cumulativeMass',cumsum(volumes)/sum(volumes),'totalVolume',sum(volumes));
        end

        function weights = verticalCellWeights(z,bounds)
            arguments (Input)
                z (:,1) double {mustBeReal,mustBeFinite}
                bounds (1,2) double {mustBeReal,mustBeFinite}
            end
            % Midpoint cells, clipped to the physical lower and upper walls.
            % These are cell volumes per unit horizontal area, not a claim
            % about the quadrature rule used by a production transform.
            if isempty(z) || any(diff(z)<=0) || bounds(1)>=bounds(2) || z(1)<bounds(1) || z(end)>bounds(2)
                error('DensityDiagnosticReference:InvalidVerticalGrid','Supply increasing heights within increasing physical bounds.');
            end
            edges = [bounds(1);(z(1:end-1)+z(2:end))/2;bounds(2)];
            weights = diff(edges);
        end

        function [mappedX,mappedZ] = polarTwist(x,z,center,radius,angle)
            arguments (Input)
                x double {mustBeReal,mustBeFinite}
                z double {mustBeReal,mustBeFinite}
                center (1,2) double {mustBeReal,mustBeFinite}
                radius (1,1) double {mustBePositive,mustBeFinite}
                angle (1,1) double {mustBeReal,mustBeFinite}
            end
            if ~isequal(size(x),size(z))
                error('DensityDiagnosticReference:CoordinateSizeMismatch','Supply x and z arrays of the same size.');
            end
            % In polar coordinates r'=r, theta'=theta+f(r), so the Jacobian
            % preserves r dr dtheta exactly. The disk must lie inside the
            % caller's domain. Keeping y fixed extends this to a volume map.
            dx = x-center(1);
            dz = z-center(2);
            theta = angle*max(0,1-(dx.^2+dz.^2)/radius^2).^3;
            mappedX = center(1)+cos(theta).*dx-sin(theta).*dz;
            mappedZ = center(2)+sin(theta).*dx+cos(theta).*dz;
            outside = dx.^2+dz.^2 >= radius^2;
            mappedX(outside) = x(outside);
            mappedZ(outside) = z(outside);
        end

        function [materialX,materialHeight] = inversePolarTwist(x,z,center,radius,angle)
            % Radius is invariant, so the inverse requires no root solve.
            [materialX,materialHeight] = DensityDiagnosticReference.polarTwist(x,z,center,radius,-angle);
        end

        function [eta,ape] = linearReference(z,materialHeight,N2)
            arguments (Input)
                z double {mustBeReal,mustBeFinite}
                materialHeight double {mustBeReal,mustBeFinite}
                N2 (1,1) double {mustBePositive,mustBeFinite}
            end
            if ~isequal(size(z),size(materialHeight))
                error('DensityDiagnosticReference:CoordinateSizeMismatch','Supply current and material heights of the same size.');
            end
            eta = z-materialHeight;
            ape = N2*eta.^2/2; % Available potential energy per unit mass.
        end

        function profile = nonlinearRestProfile(options)
            arguments (Input)
                options.N2 (1,1) double {mustBePositive,mustBeFinite} = .04
                options.curvature (1,1) double {mustBeReal,mustBeFinite} = .3
                options.depth (1,1) double {mustBePositive,mustBeFinite} = 1
                options.rho0 (1,1) double {mustBePositive,mustBeFinite} = 1025
                options.gravity (1,1) double {mustBePositive,mustBeFinite} = 9.81
            end
            % Stable quadratic profile on [-depth,depth]. Its density and
            % inverse are analytic, independent of production interpolation.
            if abs(options.curvature) >= 1
                error('DensityDiagnosticReference:UnstableProfile','Require abs(curvature)<1 to keep N2 positive throughout the domain.');
            end
            profile = options;
            profile.bounds = [-options.depth options.depth];
            profile.buoyancy = @(z) options.N2*(z+options.curvature*z.^2/(2*options.depth));
            profile.density = @(z) options.rho0*(1-profile.buoyancy(z)/options.gravity);
            profile.localN2 = @(z) options.N2*(1+options.curvature*z/options.depth);
            profile.availablePotentialEnergy = @(z,s) options.N2*(z-s).^2/2 .* ...
                (1+options.curvature*(s+(z-s)/3)/options.depth);
            % Rationalized quadratic root also handles curvature=0.
            profile.heightForDensity = @(rho) 2*(options.gravity*(1-rho/options.rho0)/options.N2) ./ ...
                (1+sqrt(1+2*options.curvature*(options.gravity*(1-rho/options.rho0)/options.N2)/options.depth));
        end
    end
end
