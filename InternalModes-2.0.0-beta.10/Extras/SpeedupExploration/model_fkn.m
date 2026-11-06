function yhat = model_fkn(theta, k, N)
%MODEL_FKN  f(k,N) = A N^m (1 + k/(k0 N))^n
% theta = [A, k0, m, n]

A  = theta(1);
k0 = theta(2);
m  = theta(3);
n  = theta(4);

% Assumes N>0 and (1 + k/(k0 N))>0 (true if k>=0, k0>0, N>0)
yhat = A .* (N.^m) .* (1 + (k ./ (k0.*N))).^n;
end