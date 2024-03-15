clc
clear

% Define your parameters
theta = [0.3, -0.5, 0.5, 0.2, 0.4];
u_mean = 1;
u_variance = 0.12;

% Initialize a matrix to store all the Y values
Y_all = zeros(100, 100);


    % Generate half of u values from N(u_mean, u_variance)
    half_size = 100/2;
    u_half1 = u_mean + sqrt(u_variance) * randn(half_size, 1);

    % Generate half of u values from N(-u_mean, u_variance)
    u_half2 = -u_mean + sqrt(u_variance) * randn(half_size, 1);
    u = [u_half1; u_half2];

    % Generate X matrix
    X = [ones(100,1), u, u.^2/2, u.^3/3, u.^4/4];

for i = 1:100
    % Generate noise N(mu, sigma)
    mu = zeros(100, 1);
    sigma =  ones(100, 1);
    N = mu + sigma .* randn(100, 1);

    % Calculate Y
    Y = X * theta' + N;
    
    % Store Y in the matrix
    Y_all(:, i) = Y;
end

% Perform LS estimation
theta_hat_LS = inv(X' * X) * X' * Y_all;

 % Calculate Y_hat
 Y_hat = X * theta_hat_LS;

% Sort X
[~, idx] = sort(X(:, 2)); % Sorting according to the second column (u)
X_sorted = X(idx, :);
 
% Sort Y according to the indexes of sorted X
Y_sorted = Y_hat(idx, :);
 
% Calculate the average of each column of Y_sorted
Y_avg = mean(Y_sorted, 2); % Take the mean along the rows (each column)
 
% Calculate the variance along columns (variance of each row)
variance_along_columns = var(Y_hat, 0, 2);

% Plot the average of Y_sorted with respect to u
figure;
plot(u(idx, :), Y_avg, '-');
xlabel('u');
ylabel('Average of Y\_sorted');
title('Plot of Average of Y\_sorted vs. u');
hold on
plot(u(idx, :), Y_avg + variance_along_columns , '-');
plot(u(idx, :), Y_avg - variance_along_columns , '-');
%scatter(u(idx, :),Y_sorted)
errorbar(u(idx(1:10:end), :), Y_sorted(1:10), variance_along_columns(1:10), 'o-');
