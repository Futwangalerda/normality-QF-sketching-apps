% sketching methods
classdef sketchingMethods
    methods(Static)
        % orthogonal projection
        function [sigmai, squareprojveci] = orth(i, m, X, c)
            [n, ~] = size(X);
            S = orth(randn(n, m))'; % S: m * n
            X_t = S * X;
            [~, D_t, V_t] = svd(X_t, 'econ');
            lambda = diag(D_t);
            sigmai = lambda(i).^2;
            RSVi = V_t(:, i);
            squareprojveci = sum(c * RSVi);
        end
        % iid Gaussian projection
        function [sigmai, squareprojveci] = gauss(i, m, X, c)
            [n, ~] = size(X);
            S = randn(m, n) / sqrt(m); % S: m * n
            X_t = S * X;
            [~, D_t, V_t] = svd(X_t, 'econ');
            lambda = diag(D_t);
            sigmai = lambda(i).^2;
            RSVi = V_t(:, i);
            squareprojveci = sum(c * RSVi);
        end
        % Sparse iid projection
        function [sigmai, squareprojveci] = iidsparse(i, X, c, S)
            X_t = S * X;
            [~, D_t, V_t] = svd(X_t, 'econ');
            lambda = diag(D_t);
            sigmai = lambda(i).^2;
            RSVi = V_t(:, i);
            squareprojveci = sum(c * RSVi);
        end
        % uniform subsampling
        function [sigmai, squareprojveci] = unifsamp(i, m, X, c)
            [n, ~] = size(X);
            sampled_ind = binornd(1, m / n, n, 1);
            X_t = X(sampled_ind == 1, :) / sqrt(m / n);
            [~, D_t, V_t] = svd(X_t, 'econ');
            lambda = diag(D_t);
            sigmai = lambda(i).^2;
            RSVi = V_t(:, i);
            squareprojveci = sum(c * RSVi);
        end
        % SRHT
        function [sigmai, squareprojveci] = srht(i, m, X, c)
            [n, p] = size(X);
            n = pow2(floor(log2(n)));
            X = (2 * binornd(1, 1 / 2, n, 1) - 1) .* X;
            X_t = zeros(n , p);
            for j = 1:p
                X_t(:, j) =  fwht(X(1:n, j)) * sqrt(n) / sqrt(m / n);
            end
            X_t = X_t(binornd(1, m / n, n, 1) == 1, :);
            [~, D_t, V_t] = svd(X_t, 'econ');
            lambda = diag(D_t);
            sigmai = lambda(i).^2;
            RSVi = V_t(:, i);
            squareprojveci = sum(c * RSVi);
        end
        
        % Sparse Sign Embedding
        % use "mex sparsesign.c" to input sparsesign
        function [sigmai, squareprojveci] = sse(i, m, X, zeta, c)
            [n, ~] = size(X);
            S = sparsesign(m, n, zeta);
            X_t = S * X;
            [~, D_t, V_t] = svd(X_t, 'econ');
            lambda = diag(D_t);
            sigmai = lambda(i).^2;
            RSVi = V_t(:, i);
            squareprojveci = sum(c * RSVi);
        end
        
    end
end