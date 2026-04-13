function D = FiniteDifferenceMatrix(numDerivs, x, leftBCDerivs, rightBCDerivs, orderOfAccuracy)
% Assemble a finite-difference matrix for $$\frac{d^m y}{dx^m}$$.
%
% `FiniteDifferenceMatrix` builds the discrete operator `D` on the grid
% `x` such that
%
% $$
% D y \approx \frac{d^m y}{dx^m},
% $$
%
% where `m = numDerivs` in the interior rows. The first and last rows use
% derivative orders `leftBCDerivs` and `rightBCDerivs` to encode boundary
% derivative conditions on the same grid. The returned matrix follows the
% existing arbitrary-grid Fornberg-weight construction used in this
% repository.
%
% - Topic: Integrators
% - Declaration: D = FiniteDifferenceMatrix(numDerivs,x,leftBCDerivs,rightBCDerivs,orderOfAccuracy)
% - Parameter numDerivs: nonnegative derivative order $$m$$ approximated in the interior rows
% - Parameter x: strictly increasing grid vector with length `n`
% - Parameter leftBCDerivs: nonnegative derivative order enforced in the first row
% - Parameter rightBCDerivs: nonnegative derivative order enforced in the last row
% - Parameter orderOfAccuracy: positive stencil-accuracy target
% - Returns D: `n x n` finite-difference matrix on the grid `x`
arguments
    numDerivs (1,1) double {mustBeInteger, mustBeNonnegative}
    x {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector}
    leftBCDerivs (1,1) double {mustBeInteger, mustBeNonnegative}
    rightBCDerivs (1,1) double {mustBeInteger, mustBeNonnegative}
    orderOfAccuracy (1,1) double {mustBeInteger, mustBePositive}
end

x = x(:);
n = length(x);
if n < 2
    error("FiniteDifferenceMatrix:GridTooShort", "x must contain at least two grid points.");
end

if any(diff(x) <= 0)
    error("FiniteDifferenceMatrix:GridNotStrictlyIncreasing", "x must be strictly increasing.");
end

if orderOfAccuracy + leftBCDerivs > n
    error("FiniteDifferenceMatrix:LeftStencilTooWide", "The left boundary stencil requires %d points, but x only contains %d points.", orderOfAccuracy + leftBCDerivs, n);
end

if orderOfAccuracy + rightBCDerivs > n
    error("FiniteDifferenceMatrix:RightStencilTooWide", "The right boundary stencil requires %d points, but x only contains %d points.", orderOfAccuracy + rightBCDerivs, n);
end

D = zeros(n,n);

% left boundary condition
range = 1:(orderOfAccuracy+leftBCDerivs); % not +1 because we're computing inclusive
c = weights( x(1), x(range), leftBCDerivs );
D(1,range) = c(leftBCDerivs+1,:);

% central derivatives, including possible weird end points
centralBandwidth = ceil(numDerivs/2)+ceil(orderOfAccuracy/2)-1;
for i=2:(n-1)
    rangeLength = 2*centralBandwidth; % not +1 because we're computing inclusive
    startIndex = max(i-centralBandwidth, 1);
    endIndex = startIndex+rangeLength;
    if (endIndex > n)
        endIndex = n;
        startIndex = endIndex-rangeLength;
    end
    range = startIndex:endIndex;
    c = weights( x(i), x(range), numDerivs );
    D(i,range) = c(numDerivs+1,:);
end

% right boundary condition
range = (n-(orderOfAccuracy+rightBCDerivs-1)):n; % not +1 because we're computing inclusive
c = weights( x(n), x(range), rightBCDerivs );
D(n,range) = c(rightBCDerivs+1,:);
end

function c = weights(z,x,m)
% Calculates FD weights. The parameters are:
%  z   location where approximations are to be accurate,
%  x   vector with x-coordinates for grid points,
%  m   highest derivative that we want to find weights for
%  c   array size m+1,lentgh(x) containing (as output) in 
%      successive rows the weights for derivatives 0,1,...,m.
%
% Taken from Bengt Fornberg
%
    n=length(x); c=zeros(m+1,n); c1=1; c4=x(1)-z; c(1,1)=1;
    for i=2:n
       mn=min(i,m+1); c2=1; c5=c4; c4=x(i)-z;
       for j=1:i-1
          c3=x(i)-x(j);  c2=c2*c3;
          if j==i-1 
             c(2:mn,i)=c1*((1:mn-1)'.*c(1:mn-1,i-1)-c5*c(2:mn,i-1))/c2;
             c(1,i)=-c1*c5*c(1,i-1)/c2;
          end
          c(2:mn,j)=(c4*c(2:mn,j)-(1:mn-1)'.*c(1:mn-1,j))/c3;
          c(1,j)=c4*c(1,j)/c3;
       end
       c1=c2;
    end

end
