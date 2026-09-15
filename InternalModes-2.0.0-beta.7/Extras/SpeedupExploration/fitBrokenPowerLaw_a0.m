function fit = fitBrokenPowerLaw_a0(x, y, opts)
% Fit y = A * (1 + (x/x0)^k)^(b/k)   (a=0 smooth broken power law)
%
% Params: theta = [A, x0, b, k]
%
% Fits in log-space by default: minimize log(yhat) - log(y).

arguments
    x (:,1) double
    y (:,1) double
    opts.UseLog10 (1,1) logical = false
    opts.Display  (1,1) string  = "off"
end

% --- keep positive, finite
ok = isfinite(x) & isfinite(y) & (x > 0) & (y > 0);
x = x(ok); y = y(ok);

if opts.UseLog10
    lg = @log10;
else
    lg = @log;
end

% --- numerically stable evaluation of log(yhat)
% log yhat = log A + (b/k)*log(1 + (x/x0)^k)
logModel = @(th, xx) ...
    lg(th(1)) + (th(3)./th(4)) .* lg( 1 + (xx./th(2)).^th(4) );

v = lg(y);

% --- initial guesses
A0 = quantile(y, 0.1);                       % plateau ~ low quantile
A0 = max(A0, eps);

% x0 guess: where y rises above plateau by ~50%
yr = y / A0;
i0 = find(yr > 1.5, 1, "first");
if isempty(i0)
    x0_0 = median(x);
else
    x0_0 = x(i0);
end

% b guess from high-x slope in log-log
[~, idx] = sort(x);
ktail = max(5, round(0.3*numel(x)));
tail = idx(end-ktail+1:end);
p = polyfit(lg(x(tail)), lg(y(tail)), 1);
b0 = max(p(1), 0);

k0 = 2; % moderate smoothness

theta0 = [A0, x0_0, b0, k0];

% --- bounds
lb = [min(y)*1e-3, min(x),  0,  0.1];
ub = [max(y)*1e3,  max(x), 50, 50];

optsLSQ = optimoptions("lsqcurvefit", ...
    "Display", opts.Display, ...
    "MaxIterations", 2000, ...
    "MaxFunctionEvaluations", 5e4);

% --- fit: lsqcurvefit expects model(theta,x) ~ v
[theta, resnorm, resid, exitflag, output, lambda, J] = ...
    lsqcurvefit(@(th,xx) logModel(th,xx), theta0, x, v, lb, ub, optsLSQ);

% --- pack outputs
fit.theta = theta;                 % [A, x0, b, k]
fit.model = @(xx) theta(1) .* (1 + (xx./theta(2)).^theta(4)).^(theta(3)/theta(4));
fit.logResiduals = resid;
fit.resnorm = resnorm;
fit.exitflag = exitflag;
fit.output = output;
fit.lambda = lambda;
fit.jacobian = J;

% quick approximate covariance in log-space
dof = max(numel(v) - numel(theta), 1);
sigma2 = resnorm / dof;
fit.cov = sigma2 * inv(J.'*J);
fit.se  = sqrt(diag(fit.cov));
fit.ci95 = [theta(:) - 1.96*fit.se, theta(:) + 1.96*fit.se];

end