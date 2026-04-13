function x = leftSolve(A, b)
% Solve a linear system, falling back to a pseudoinverse if needed.
%
% This is a defensive wrapper around left division. Sparse systems are
% solved directly with `A\b`. Dense systems first check the reciprocal
% condition number; if the matrix is effectively singular, the function
% falls back to `pinv(A)*b`.
%
% - Topic: Solve fit systems
% - Developer: true
% - Declaration: x = leftSolve(A,b)
% - Parameter A: system matrix
% - Parameter b: right-hand side
% - Returns x: solution vector or matrix
if isempty(A)
    x = zeros(size(b));
    return;
end

if issparse(A)
    x = A\b;
    return;
end

reciprocalCondition = rcond(A);
if reciprocalCondition < eps(class(A))
    x = pinv(A) * b;
else
    x = A\b;
end
