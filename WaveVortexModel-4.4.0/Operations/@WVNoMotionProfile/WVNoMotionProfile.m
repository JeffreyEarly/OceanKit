classdef WVNoMotionProfile
    % Evaluate a monotone density profile and its displacement energetics.
    %
    % A shape-preserving cubic interpolant supplies the density, its inverse,
    % and the exact polynomial integral used by APE. Density plateaus do not
    % have a unique inverse and are rejected explicitly.
    %
    % - Topic: Internal
    % - Developer: true
    properties (SetAccess=private)
        z
        rho
    end
    properties (Access=private)
        densityOffset
        densityScale
        coefficients
        normalizedDensity
    end
    methods
        function self = WVNoMotionProfile(z,rho)
            arguments
                z (:,1) double {mustBeFinite,mustBeReal}
                rho (:,1) double {mustBeFinite,mustBeReal}
            end
            if numel(z) < 2 || numel(z) ~= numel(rho) || any(diff(z) <= 0)
                error('WVNoMotionProfile:InvalidGrid','Use matching density and strictly increasing height vectors with at least two entries.');
            end
            if any(diff(rho) >= 0)
                error('WVNoMotionProfile:NonInvertibleDensity','Density must strictly decrease with height; plateaus have no unique inverse material height.');
            end
            self.z = z;
            self.rho = rho;
            self.densityOffset = rho(end);
            self.densityScale = rho(1)-rho(end);
            self.normalizedDensity = (rho-self.densityOffset)/self.densityScale;
            pp = pchip(z,self.normalizedDensity);
            self.coefficients = [zeros(pp.pieces,4-pp.order),pp.coefs];
        end

        function rho = density(self,z)
            arguments
                self
                z double {mustBeFinite,mustBeReal}
            end
            self.validateHeight(z);
            rho = self.densityOffset + self.densityScale * ppval(mkpp(self.z,self.coefficients),z);
        end

        function z = inverse(self,rho)
            arguments
                self
                rho double {mustBeFinite,mustBeReal}
            end
            shape = size(rho);
            tolerance = 8*eps(max(abs(self.rho)));
            if any(rho(:) < self.rho(end)-tolerance | rho(:) > self.rho(1)+tolerance)
                error('WVNoMotionProfile:DensityOutsideProfile','Total density lies outside the no-motion profile range; recompute a profile that spans the current density distribution.');
            end
            heightTolerance = 8*eps(max(abs(self.z))+diff(self.z([1 end])));
            z = zeros(shape);
            % Bound coefficient and iteration arrays while keeping contiguous
            % arithmetic within each block. Converged heights stay frozen
            % under the same safeguarded Newton steps and stopping criteria.
            blockSize = 1048576;
            for first = 1:blockSize:numel(rho)
                positions = first:min(first+blockSize-1,numel(rho));
                values = reshape(rho(positions),[],1);
                target = (min(max(values,self.rho(end)),self.rho(1))-self.densityOffset)/self.densityScale;
                interval = discretize(-target,-self.normalizedDensity);
                base = self.z(interval);
                c = self.coefficients(interval,:);
                lower = zeros(size(target));
                upper = self.z(interval+1)-base;
                height = upper .* (self.normalizedDensity(interval)-target) ./ (self.normalizedDensity(interval)-self.normalizedDensity(interval+1));
                active = true(size(target));
                for iteration = 1:64
                    residual = ((c(:,1).*height+c(:,2)).*height+c(:,3)).*height+c(:,4)-target;
                    lower(residual > 0) = height(residual > 0);
                    upper(residual < 0) = height(residual < 0);
                    derivative = (3*c(:,1).*height+2*c(:,2)).*height+c(:,3);
                    next = height-residual./derivative;
                    useMidpoint = ~isfinite(next) | next <= lower | next >= upper;
                    next(useMidpoint) = (lower(useMidpoint)+upper(useMidpoint))/2;
                    converged = residual == 0 | upper-lower <= heightTolerance;
                    active = active & ~converged;
                    if ~any(active)
                        break
                    end
                    height(active) = next(active);
                end
                if any(active)
                    error('WVNoMotionProfile:InverseDidNotConverge','Density inversion failed to reach the height tolerance.');
                end
                height = base+height;
                % Endpoint arithmetic can round a converged local inverse
                % just outside the original knots. Preserve every interior
                % value and the strict public query bounds used by APE.
                lowerRoundoff = height < self.z(1) & height >= self.z(1)-heightTolerance;
                upperRoundoff = height > self.z(end) & height <= self.z(end)+heightTolerance;
                height(lowerRoundoff) = self.z(1);
                height(upperRoundoff) = self.z(end);
                z(positions) = height;
            end
        end

        function ape = availablePotentialEnergy(self,z,materialHeight,g,rho0)
            arguments
                self
                z double {mustBeFinite,mustBeReal}
                materialHeight double {mustBeFinite,mustBeReal}
                g (1,1) double {mustBeFinite,mustBeReal,mustBePositive}
                rho0 (1,1) double {mustBeFinite,mustBeReal,mustBePositive}
            end
            if ~isequal(size(z),size(materialHeight))
                error('WVNoMotionProfile:ShapeMismatch','Height and material height must have the same shape.');
            end
            self.validateHeight(z);
            self.validateHeight(materialHeight);
            % Exchange the two integrals to evaluate (r-z)*rho'(r)
            % from material height s to z. This uses density derivatives
            % directly, avoiding cancellation even for sub-ulp density
            % changes associated with small representable displacements.
            ape = zeros(size(z));
            % Bound temporary arrays and locate each starting interval once.
            % Only parcels that cross a knot remain active for another pass;
            % work scales with actual crossings, not every vertical interval.
            blockSize = 65536;
            for first = 1:blockSize:numel(z)
                positions = first:min(first+blockSize-1,numel(z));
                height = reshape(z(positions),[],1);
                target = reshape(materialHeight(positions),[],1);
                low = min(height,target);
                high = max(height,target);
                interval = discretize(low,self.z);
                integral = zeros(size(low));
                active = find(high > low);
                while ~isempty(active)
                    index = interval(active);
                    a = low(active);
                    b = min(high(active),self.z(index+1));
                    t = a-self.z(index);
                    d = b-a;
                    c = self.coefficients(index,:);
                    A = (3*c(:,1).*t+2*c(:,2)).*t+c(:,3);
                    B = 6*c(:,1).*t+2*c(:,2);
                    C = 3*c(:,1);
                    D = height(active)-a;
                    piece = -D.*A.*d + (A-D.*B).*d.^2/2 + (B-D.*C).*d.^3/3 + C.*d.^4/4;
                    integral(active) = integral(active)+piece;
                    remaining = b < high(active);
                    active = active(remaining);
                    low(active) = b(remaining);
                    interval(active) = interval(active)+1;
                end
                ape(positions) = (g*self.densityScale/rho0)*sign(height-target).*integral;
            end
        end
    end
    methods (Access=private)
        function validateHeight(self,z)
            if any(z(:) < self.z(1) | z(:) > self.z(end))
                error('WVNoMotionProfile:HeightOutsideProfile','Height must lie inside the no-motion profile domain.');
            end
        end
    end
end
