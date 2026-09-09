function fit = fitFlatToPowerLaw(x, y, opts)
% fitFlatToPowerLaw  Fit y(x) that transitions from flat to power law.
%
% Model:
%   y(x) = y0 + A * x^p * sigmoid( (log(x)-log(x0))/Delta )
%   sigmoid(t) = 1/(1+exp(-t))
%
% Fit is done by minimizing residuals in log-space:
%   r = log(yhat) - log(y)
%
% Inputs
%   x, y : vectors (x>0, y>0 recommended)
%   opts : (optional) struct with fields:
%          .UseLog10   (default false) use log10 instead of natural log
%          .Robust     (default true)  robust weighting (bisquare) via IRLS
%          .Display    (default 'off') lsqcurvefit display option
%
% Output
%   fit : struct with fields:
%         .theta       [y0, A, p, x0, Delta]
%         .model       function handle yhat = model(x)
%         .residuals   residual vector (log-space)
%         .exitflag, .output, .lambda, .jacobian, .cov
%         .ci95        95% approx CI for parameters (Wald)

arguments
    x (:,1) double
    y (:,1) double
    opts.UseLog10 (1,1) logical = false
    opts.Robust   (1,1) logical = true
    opts.Display  (1,1) string  = "off"
end

% --- sanitize
ok = isfinite(x) & isfinite(y) & (x > 0) & (y > 0);
x = x(ok); y = y(ok);

if opts.UseLog10
    lg = @log10;
else
    lg = @log;
end

% --- model pieces
sigmoid = @(t) 1 ./ (1 + exp(-t));

model = @(th, xx) ...
    th(1) + th(2) .* (xx.^th(3)) .* sigmoid( (lg(xx) - lg(th(4))) ./ th(5) );

% --- objective in log space (fit multiplicative errors better)
% We minimize log(yhat) - log(y). Ensure yhat stays positive by bounding.
resfun = @(th, xx) lg(model(th, xx)) ;

% For lsqcurvefit we fit v = log(y) directly:
v = lg(y);

% --- initial guesses (reasonable defaults)
% baseline y0 ~ lower quantile
y0_0 = quantile(y, 0.1);

% estimate exponent from upper tail slope in log-log (ignore baseline roughly)
% pick top 30% x values
[~, idx] = sort(x);
k = max(5, round(0.3*numel(x)));
tail = idx(end-k+1:end);
p0 = max(0, polyfit(lg(x(tail)), lg(max(y(tail)-y0_0, eps)), 1));
p0 = p0(1);

% transition location near where y rises above baseline by some fraction
yr = y - y0_0;
thr = max(quantile(yr(yr>0), 0.3), eps);
i0 = find(yr > thr, 1, "first");
if isempty(i0), i0 = round(numel(x)/2); end
x0_0 = x(i0);

% amplitude A from one point in tail
A0 = median( max(y(tail) - y0_0, eps) ./ (x(tail).^max(p0, 0)) );

% transition width in log-x units (0.3 ~ factor ~2)
Delta0 = 0.5;

theta0 = [y0_0, A0, max(p0, 0), x0_0, Delta0];

% --- bounds
% y0 >= 0, A >= 0, p >= 0, x0 within x-range, Delta positive
lb = [0,    0,    0,   min(x),  1e-3];
ub = [Inf,  Inf,  Inf, max(x),  10];

% --- solver options
options = optimoptions("lsqcurvefit", ...
    "Display", opts.Display, ...
    "MaxFunctionEvaluations", 5e4, ...
    "MaxIterations", 2e3);

% --- robust IRLS (optional)
w = ones(size(x));
theta = theta0;
for it = 1 : (opts.Robust*10 + ~opts.Robust)
    % weighted target: sqrt(w)*v, and weighted model: sqrt(w)*resfun
    v_w = sqrt(w) .* v;
    f_w = @(th, xx) sqrt(w) .* resfun(th, xx);

    [theta, resnorm, resid, exitflag, output, lambda, J] = ...
        lsqcurvefit(f_w, theta, x, v_w, lb, ub, options);

    r = resid; % already weighted residuals in log-space
    if ~opts.Robust
        break
    end

    % update bisquare weights on *unweighted* residuals
    r0 = lg(model(theta, x)) - v;

    s = 1.4826 * median(abs(r0 - median(r0))) + eps; % robust scale (MAD)
    c = 4.685; % bisquare tuning
    u = r0 ./ (c*s);
    w = (abs(u) < 1) .* (1 - u.^2).^2;

    % stop if weights stabilize
    if it > 2 && norm(r0,2) / max(norm(v,2),eps) < 1e-6
        break
    end
end

% --- covariance / CI (approx; for log-residuals)
dof = max(numel(v) - numel(theta), 1);
sigma2 = resnorm / dof;
covTheta = sigma2 * inv(J.'*J);

se = sqrt(diag(covTheta));
ci95 = [theta(:) - 1.96*se, theta(:) + 1.96*se];

fit.theta = theta;
fit.model = @(xx) model(theta, xx);
fit.residuals = lg(fit.model(x)) - lg(y);
fit.exitflag = exitflag;
fit.output = output;
fit.lambda = lambda;
fit.jacobian = J;
fit.cov = covTheta;
fit.ci95 = ci95;

end