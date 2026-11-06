function r = residual_fkn(theta, k, N, y, mode)
%RESIDUAL_FKN residual vector for fitting.
% mode: "linear" -> r = yhat - y
%       "log"    -> r = log(yhat) - log(y)   (good for multiplicative noise)

arguments
    theta (1,4) double
    k (:,1) double
    N (:,1) double
    y (:,1) double
    mode (1,1) string = "log"
end

yhat = model_fkn(theta, k, N);

switch lower(mode)
    case "linear"
        r = yhat - y;
    case "log"
        % Guard: require positive y and yhat
        if any(y <= 0) || any(yhat <= 0)
            r = inf(size(y));
        else
            r = log(yhat) - log(y);
        end
    otherwise
        error("Unknown mode: %s", mode);
end
end