function Y_boot = reconstruct_data(Y, res_boot, p, c)
    % Initialize bootstrapped data
    T = size(Y, 1);
    n = size(Y, 2);
    Y_boot = zeros(size(Y));
    Y_boot(1:p, :) = Y(1:p, :); % Use the original data for initial lags

    % Reconstruct data using resampled residuals
    for t = p+1:T
        % Estimate VAR model coefficients dynamically for the bootstrap sample
        [A_hat, C_hat] = estimate_VAR(Y_boot(1:t-1, :), p, c);

        % Stack lagged values as a column vector
        lagged_Y = reshape(flipud(Y_boot(t-p:t-1, :))', [], 1); % [n*p, 1]

        % Compute the t-th observation of Y_boot
        Y_boot(t, :) = (C_hat + A_hat * lagged_Y + res_boot(t-p, :)')'; % Ensure [1, n]
    end
end
