function varargout = deim(u)
%DEIM   Compute the discrete empirical interpolation method on the POD modes.
%   IND = DEIM(U) returns the indices of the data array that the DEIM algorithm
%   determines should be kept.  U is the matrix of POD modes on the nonlinear
%   function, as given by POD_NONLINEAR.
%
%   [IND, P] = DEIM(U) also returns the projection matrix P equal to the zero
%   matrix with the same size as U, except with 1 at the indices given by IND.
%
%   See also POD_NONLINEAR.

n_modes = size(u, 2);

ind = zeros(n_modes, 1); % The DEIM indices.

[~, ind(1)] = max(abs(u(:, 1)));

for k = 2:n_modes
    a = u(ind(1:k-1), 1:k-1);
    b = u(ind(1:k-1), k);
    % The coefficient vector mapping the modes to the nonlinear evaluation.
    c = a \ b;
    
    r = u(:, k) - u(:, 1:k-1) * c; % The residual.
    [~, ind(k)] = max(abs(r));
end

varargout(1) = {ind};

% Compute the projection matrix if needed.
if nargout == 2
    p = zeros(size(u)); % The projection matrix.
    
    for k = 1:n_modes
        p(ind(k), k) = 1;
    end
    
    varargout(2) = {p};
end
end