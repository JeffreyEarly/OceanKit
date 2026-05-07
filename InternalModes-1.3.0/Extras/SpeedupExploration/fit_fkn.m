function fit = fit_fkn(k, N, y, opts)
%FIT_FKN Fit f(k,N)=A N^m (1 + k/(k0 N))^n with lsqcurvefit.
%
% Inputs (vectors same length):
%   k, N, y
% opts:
%   .UseLogResiduals (default true)  fit in log space
%   .theta0          initial guess [A k0 m n] (optional)
%   .lb, .ub         bounds (optional)
%   .Display         lsqcurvefit display (default 'off')
%
% Output:
%   fit.theta, fit.model, fit.residuals, etc.

arguments
    k (:,1) double
    N (:,1) double
    y (:,1) double
    opts.UseLogResiduals (1,1) logical = true
    opts.theta0 (1,4) double = [NaN NaN NaN NaN]
    opts.lb (1,4) double = [0,   eps, -Inf, -Inf]
    opts.ub (1,4) double = [Inf, Inf,  Inf,  Inf]
    opts.Display (1,1) string = "off"
end

% --- keep valid points
ok = isfinite(k) & isfinite(N) & isfinite(y) & (N > 0);
if opts.UseLogResiduals
    ok = ok & (y > 0);
end
k = k(ok); N = N(ok); y = y(ok);

X = [k, N];

% --- choose fitting target (linear y or log y)
if opts.UseLogResiduals
    target = log(y);
    fun = @(th, X) log(model_fkn(th, X(:,1), X(:,2)));
else
    target = y;
    fun = @(th, X) model_fkn(th, X(:,1), X(:,2));
end

% --- initial guess (if not provided)
theta0 = opts.theta0;
if any(isnan(theta0))
    % crude but often workable:
    % plateau-ish A from smallest-k per-N neighborhood (fallback: median(y))
    A0 = median(y, "omitnan");
    A0 = max(A0, eps);

    % k0 scale: typical k/N ratio (since transition uses k/(k0 N))
    r = k ./ N;
    r = r(isfinite(r) & r>0);
    if isempty(r), k0_0 = 1; else, k0_0 = median(r); end

    % m: slope vs N at small k (rough guess 0)
    m0 = 0;

    % n: slope vs k at large k/N (rough guess 1)
    n0 = 1;

    theta0 = [A0, k0_0, m0, n0];
end

% --- bounds: enforce k0>0 and A>=0 by default
lb = opts.lb;
ub = opts.ub;
lb(2) = max(lb(2), eps);

% --- optimize
opt = optimoptions("lsqcurvefit", ...
    "Display", opts.Display, ...
    "MaxIterations", 2000, ...
    "MaxFunctionEvaluations", 5e4);

[theta, resnorm, resid, exitflag, output, lambda, J] = ...
    lsqcurvefit(fun, theta0, X, target, lb, ub, opt);

% --- pack
fit.theta = theta;                         % [A k0 m n]
fit.model = @(kk,NN) model_fkn(theta, kk, NN);
fit.resnorm = resnorm;
fit.residuals = resid;
fit.exitflag = exitflag;
fit.output = output;
fit.lambda = lambda;
fit.jacobian = J;

end