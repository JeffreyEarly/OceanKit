function varargout = interpolatedFieldAtPosition(self,x,y,z,method,varargin)
    % interpolate gridded fields at arbitrary positions
    %
    % Two-dimensional fields are interpolated periodically in (x,y), with
    % z ignored. Three-dimensional fields retain the existing (x,y,z)
    % interpolation behavior.
    %
    % - Topic: State Variables
    if nargin-5 ~= nargout
        error('You must have the same number of input variables as output variables');
    end

    % (x,y) are periodic for the gridded solution
    x_tilde = mod(x,self.Lx);
    y_tilde = mod(y,self.Ly);

    dx = self.x(2)-self.x(1);
    x_index = floor(x_tilde/dx);
    dy = self.y(2)-self.y(1);
    y_index = floor(y_tilde/dy);

    if strcmp(method,'linear') && nargout > 0 && all(~cellfun(@ismatrix,varargin)) && all(isfinite(x(:))) && all(isfinite(y(:)))
        x_weight = x_tilde/dx-x_index;
        x_lower = x_index+1;
        x_upper = mod(x_lower,self.Nx)+1;
        y_weight = y_tilde/dy-y_index;
        y_lower = y_index+1;
        y_upper = mod(y_lower,self.Ny)+1;
        horizontalIndex = [x_lower(:)+self.Nx*(y_lower(:)-1) x_upper(:)+self.Nx*(y_lower(:)-1) x_lower(:)+self.Nx*(y_upper(:)-1) x_upper(:)+self.Nx*(y_upper(:)-1)];
        horizontalWeight = [(1-x_weight(:)).*(1-y_weight(:)) x_weight(:).*(1-y_weight(:)) (1-x_weight(:)).*y_weight(:) x_weight(:).*y_weight(:)];

        z_lower = discretize(z,self.z);
        outsideVerticalDomain = isnan(z_lower);
        z_lower(outsideVerticalDomain) = 1;
        z_weight = (z-reshape(self.z(z_lower),size(z)))./reshape(self.z(z_lower+1)-self.z(z_lower),size(z));
        z_weight(outsideVerticalDomain) = 0;
        lowerIndex = horizontalIndex+self.Nx*self.Ny*(z_lower(:)-1);
        upperIndex = lowerIndex+self.Nx*self.Ny;

        varargout = cell(1,nargout);
        for i = 1:nargout
            U = varargin{i};
            lowerPlane = sum(horizontalWeight.*U(lowerIndex),2);
            upperPlane = sum(horizontalWeight.*U(upperIndex),2);
            u = reshape(lowerPlane+(upperPlane-lowerPlane).*z_weight(:),size(x));
            u(outsideVerticalDomain) = 0;
            u(isnan(z)) = NaN;
            varargout{i} = u;
            if any(isnan(u(:)))
                fprintf('bad apple.\n')
            end
        end
        return
    end

    % Identify the particles along the interpolation boundary
    if strcmp(method,'spline')
        S = 3+1; % cubic spline, plus buffer
    elseif strcmp(method,'linear')
        S = 1+1;
    end
    bpx = x_index < S-1 | x_index > self.Nx-S;
    bpy = y_index < S-1 | y_index > self.Ny-S;

    varargout = cell(1,nargout);
    for i = 1:nargout
        U = varargin{i}; % gridded field
        u = zeros(size(x)); % interpolated value
        if ismatrix(U)
            X = self.X(:,:,end);
            Y = self.Y(:,:,end);
            u(~bpx & ~bpy) = interpn(X,Y,U,x_tilde(~bpx & ~bpy),y_tilde(~bpx & ~bpy),method,0);
            if any(bpx & bpy)
                x_tildeS = mod(x(bpx & bpy)+S*dx,self.Lx);
                y_tildeS = mod(y(bpx & bpy)+S*dy,self.Ly);
                u(bpx & bpy) = interpn(X,Y,circshift(U,[S S]),x_tildeS,y_tildeS,method,0);
            end

            if any(bpx & ~bpy)
                x_tildeS = mod(x(bpx & ~bpy)+S*dx,self.Lx);
                u(bpx & ~bpy) = interpn(X,Y,circshift(U,[S 0]),x_tildeS,y_tilde(bpx & ~bpy),method,0);
            end

            if any(~bpx & bpy)
                y_tildeS = mod(y(~bpx & bpy)+S*dy,self.Ly);
                u(~bpx & bpy) = interpn(X,Y,circshift(U,[0 S]),x_tilde(~bpx & bpy),y_tildeS,method,0);
            end
        else
            u(~bpx & ~bpy) = interpn(self.X,self.Y,self.Z,U,x_tilde(~bpx & ~bpy),y_tilde(~bpx & ~bpy),z(~bpx & ~bpy),method,0);
            if any(bpx & bpy)
                x_tildeS = mod(x(bpx & bpy)+S*dx,self.Lx);
                y_tildeS = mod(y(bpx & bpy)+S*dy,self.Ly);
                u(bpx & bpy) = interpn(self.X,self.Y,self.Z,circshift(U,[S S 0]),x_tildeS,y_tildeS,z(bpx & bpy),method,0);
            end

            if any(bpx & ~bpy)
                x_tildeS = mod(x(bpx & ~bpy)+S*dx,self.Lx);
                u(bpx & ~bpy) = interpn(self.X,self.Y,self.Z,circshift(U,[S 0 0]),x_tildeS,y_tilde(bpx & ~bpy),z(bpx & ~bpy),method,0);
            end

            if any(~bpx & bpy)
                y_tildeS = mod(y(~bpx & bpy)+S*dy,self.Ly);
                u(~bpx & bpy) = interpn(self.X,self.Y,self.Z,circshift(U,[0 S 0]),x_tilde(~bpx & bpy),y_tildeS,z(~bpx & bpy),method,0);
            end
        end
        varargout{i} = u;
        if any(isnan(u(:)))
            fprintf('bad apple.\n')
        end
    end
end
