function [AA, Sigma, res, logL] = estimateVAR_LagY(Y, p, c)
    % Estimate VAR(p) model using the [Y X] = LagY(y,p) method
    T = size(Y, 1);
    n = size(Y, 2);
    
    % Create lagged matrices for Y and X
    [Ydep, X] = LagY(Y, p);
    
    % Add constant if specified
    if c == 1
        X = [ones(size(X, 1), 1), X];
    end
    
    % Estimate coefficients using [AA = inv(X'X)X'Y]
    AA = (X' * X) \ (X' * Ydep);
    
    % Calculate residuals
    res = Ydep - X * AA;
    
    % Compute the covariance matrix of residuals
    Sigma = (res' * res) / (T - p - c);
    
    % Compute log-likelihood
    logL = -0.5 * (T - p) * (n * log(2 * pi) + log(det(Sigma)) + 1);
end

function [Ydep, X] = LagY(Y, p)
    % Create lagged data matrices
    T = size(Y, 1);
    n = size(Y, 2);
    
    % Create lagged observations for dependent variable Ydep and regressors X
    Ydep = Y(p+1:end, :);
    X = [];
    for lag = 1:p
        X = [X, Y(p+1-lag:end-lag, :)];
    end
end
