% Generate data from multivariate Gaussian pdfs
% Author : Wei Chu @ UCLA
% Date : 06/26/2008

% n   - # of samples to generate
% dim - # of dimensionsf
% mix - # of mixtures
% wt  - weight       (1 x mix)
% mu  - mean         (dim x mix)
% cov - covariance   (dim x dim x mix)

function [y] = mgdistri(n, dim, mix, wt, mu, cov)


% check input params
[row, col] = size(wt);
if row ~= 1 || col ~= mix
    error('wt  - weight should be: a mix x 1 matrix')
end
[row, col] = size(mu);
if row ~= dim || col ~= mix
    error('mu  - mean should be: a mix x dim matrix')
end
[id, row, col] = size(cov);
if id ~= dim || row ~= dim || col ~= mix
    error('cov - covariance should be: a mix x dim x dim matrix')
end

% rising ladder
ladder = zeros(mix, 1);
for j = 1:mix
    if j == 1
        ladder(j) = wt(j);
    else
        ladder(j) = ladder(j-1) + wt(j);
    end
end

% cholesky decomposition for all the covariance matrices: A = LL'
lcho = zeros(dim, dim, mix);
for j = 1:mix
    [lcho(:,:,j), p] = chol(cov(:,:,j), 'lower');
    if p ~= 0
        error('cov - covariance should be PSD (actually symmetric)')
    end
end

y = zeros(dim, n);
% for each sample
for i = 1:n
    % generate from which mixture ? 
    r = rand();
    for j = 1:mix
        if r < ladder(j)
            break;
        end
    end
    % generate from the jth mixture
    % N (x; 0, I) -> N (y; mu, Sigma) ; Sigma = LL'
    % y = mu + Lx 
    x = randn(dim,1);
    y(:,i) = lcho(:,:,j)*x + mu(:,j);
end



