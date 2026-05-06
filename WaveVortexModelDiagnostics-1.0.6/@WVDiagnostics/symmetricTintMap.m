function cmap = symmetricTintMap(c,options)
%SYMMETRICTINTMAP Colormap going from color -> tinted white -> color
%
% This function is used by plotPoissonFlowOverContours
%
%   c   : 1x3 RGB color in [0 1]
%   N   : total number of levels (even recommended)
%   tintStrength : fraction of c mixed into white at the center
%
% Example:
%   c = [0.2 0.6 0.8];
%   colormap(symmetricTintMap(c,256,0.05)); colorbar
arguments
    c
    options.N = 256
    options.tintStrength = 0.05;
end

N = options.N;
tintStrength = options.tintStrength;

% midpoint color: mostly white with a hint of c
c_mid = (1 - tintStrength)*[1 1 1] + tintStrength*c;

% split N into two halves
n1 = floor(N/2);
n2 = N - n1;

% first half: c -> c_mid
half1 = [linspace(c(1), c_mid(1), n1)' ...
    linspace(c(2), c_mid(2), n1)' ...
    linspace(c(3), c_mid(3), n1)'];

% second half: c_mid -> c
half2 = [linspace(c_mid(1), c(1), n2)' ...
    linspace(c_mid(2), c(2), n2)' ...
    linspace(c_mid(3), c(3), n2)'];

cmap = [half1; half2];
end