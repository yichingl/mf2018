% draw contours for multivariate Gaussian pdfs
% Author : Wei Chu @ UCLA
% Date : 06/26/2008

% dim - # of dimensionsf
% mix - # of mixtures
% wt  - weight       (1 x mix)
% mu  - mean         (dim x mix)
% cov - covariance   (dim x dim x mix)

function [flag] = plotmgd(dim, mix, wt, mu, cov)
x1 = linspace(-2, .1, 2);
x2 = linspace(-2, .1, 2);
[X1,X2] = meshgrid(x1, x2);
num = size(X1(:));
U = [X1(:) X2(:)]; % N * dim
V = zeros(num, 1);
figure;
for j = 1:mix
    Z = mvnpdf(U, mu(:, j)', cov(:,:,j));
    V = V + wt(j)*Z;
end
coutour(X1, X2, V); hold on;
        
