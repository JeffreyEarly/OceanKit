function B = matrixForDataPoints(dataPoints, options)
% Evaluate terminated B-spline basis functions and optional derivatives.
%
% Use this to assemble a design matrix for interpolation, regression, or
% direct inspection of the basis functions.
%
% For `D=0`, the returned matrix satisfies
%
% $$
% B_{ij} = B_{j,S}(t_i;\tau), \qquad x(t_i) \approx \sum_{j=1}^{M} B_{ij}\,\xi_j.
% $$
%
% When `D > 0`, slice `B(:,:,d+1)` stores the basis values for derivative
% order `d`, so `B(:,:,d+1) * xi` evaluates the `d`th derivative at
% `dataPoints`.
%
% ```matlab
% B = BSpline.matrixForDataPoints(t, knotPoints=knotPoints, S=3);
% xi = B \ x;
% spline = BSpline(S=3, knotPoints=knotPoints, xi=xi);
% ```
%
% - Topic: Build spline bases
% - Declaration: B = matrixForDataPoints(dataPoints, options)
% - Parameter dataPoints: points at which to evaluate the splines
% - Parameter options.knotPoints: spline knot points
% - Parameter options.S: spline degree
% - Parameter options.D: (optional) number of spline derivatives to return, max(D)=S
% - Returns B: array of size `numel(dataPoints) x M x (D+1)` where `M = numel(knotPoints) - S - 1`
arguments
    dataPoints (:,1) double {mustBeNumeric,mustBeReal}
    options.knotPoints (:,1) double {mustBeNumeric,mustBeReal}
    options.S (1,1) double {mustBeInteger,mustBeNonnegative}
    options.D (1,1) double {mustBeInteger,mustBeNonnegative} = 0
end
knotPoints = options.knotPoints;
S = options.S;
K = S + 1;

if any(diff(knotPoints) < 0)
    error('BSpline:InvalidKnotPoints', 'knotPoints must be non-decreasing.');
end

D = options.D;

% number of knots
M = length(knotPoints);

nl = find(knotPoints <= knotPoints(1),1,'last');
nr = M - find(knotPoints == knotPoints(end),1,'first')+1;
if (nl < K || nr < K)
    error('Your splines are not terminated. You need to have K repeat knot points at the beginning and end.');
end

% This is true assuming the original knotPoints were strictly monotonically
% increasing (no repeat knots) and we added repeat knots at the beginning
% and end of the sequences.
N_splines = M - K;

% number of collocation points
N = length(dataPoints);

% 1st index is the N collocation points
% 2nd index is the the M splines
B = zeros(N,N_splines,D+1); % This will contain all splines and their derivatives
delta_r = zeros(N,K);
delta_l = zeros(N,K);
knot_indices = discretize(dataPoints,knotPoints(1:(M-K+1)));

% XB will contain all splines from (K-D) through order (K-1).
% These are needed to compute the derivatives of the spline, if
% requested by the user.
%
% The indexing is such that the spline of order m, is located
% at index m-(K-D)+1
if D > 0
    XB = zeros(N,N_splines,D);
    if D + 1 == K % if we go tho the max derivative, we need to manually create the 0th order spline.
        to_indices = ((1:N)' + size(XB,1)*( (knot_indices-1) + size(XB,2)*(1-1)));
        XB(to_indices) = 1;
    end
end

b = zeros(N,K);
b(:,1) = 1;
for j=1:(K-1) % loop through splines of increasing order: j+1
    delta_r(:,j) = knotPoints(knot_indices+j) - dataPoints;
    delta_l(:,j) = dataPoints - knotPoints(knot_indices+1-j);

    saved = zeros(N,1);
    for r=1:j % loop through the nonzero splines
        term = b(:,r)./(delta_r(:,r)+delta_l(:,j+1-r));
        b(:,r) = saved + delta_r(:,r).*term;
        saved = delta_l(:,j+1-r).*term;
    end
    b(:,j+1) = saved;

    % Save this info for later use in computing the derivatives
    % have to loop through one index.
    if j+1 >= K-D && j+1 <= K-1 % if K-j == 1, we're at the end, j+1 is the spline order, which goes into slot j+1-(K-1-D). Thus, if j+1=K, and D=1, this goes in slot 2
        for r = 1:(j+1)
            % (i,j,k) = (1,knot_indices-j+(r-1),j+1) --- converted to linear indices
            to_indices = ((1:N)' + size(XB,1)*( (knot_indices-j+(r-1)-1) + size(XB,2)*(j+1-(K-1-D)-1)));
            from_indices = (1:N)' + size(b,1) * (r-1);
            XB(to_indices) = b(from_indices);
        end
    end
end

for r = 1:K
    to_indices = ((1:N)' + size(B,1) * ( (knot_indices-(K-1) + (r-1)-1) + size(B,2)*(1-1) ));
    from_indices = (1:N)' + size(b,1) * (r-1);
    B(to_indices) = b(from_indices);
end

diff_coeff = @(a,r,m) (K-m)*(a(2)-a(1))/(knotPoints(r+K-m) - knotPoints(r));

if D > 0
    for r=1:N_splines
        % alpha mimics equation X.16 in deBoor's PGS, but localized to avoid
        % the zero elements.
        alpha = zeros(K+1,K+1); % row is the coefficient, column is the derivative (1=0 derivatives)
        alpha(2,1) = 1;
        for m=1:D % loop over derivatives
            for i=1:(m+1) % loop over coefficients
                a = alpha(:,m);
                alpha(i+1,m+1) = diff_coeff(a(i:end),r+i-1,m);
                if isinf(alpha(i+1,m+1)) || isnan(alpha(i+1,m+1))
                    alpha(i+1,m+1) = 0;
                end
                if r+i-1>N_splines
                    B0 = zeros(N,1);
                else
                    B0 = XB(:,r+i-1,D+1-m); % want the K-m order spline, in position D+1-m
                end
                B(:,r,m+1) = B(:,r,m+1) + alpha(i+1,m+1)*B0;
            end
        end
    end
end

end
