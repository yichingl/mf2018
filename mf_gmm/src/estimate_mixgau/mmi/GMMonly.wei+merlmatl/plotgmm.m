function p = plotgmm(g)
% Read proto and write data generated from Guassians
% Author : Wei Chu @ UCLA
% Date : 06/26/2008
%
% gmmdata(strModel, bplot, bsimu, strDir, trnum, tenum)
% bsimu    - 1 - simu data, 0 - do not
% bplot    - 1 - plot data, 0 - do not
% strMod   - contain multivariate Gaussian pdfs' name, mix weight, mean, covar
% trnum    - # of samples to be generated for each pdf for training set
% tenum    - # of samples to be generated for each pdf for test set
% bshow    - 1 - show the simu data, 0 - do not


gmm.mix = length(g.m);
gmm.dim = size(g.m{1},1);
gmm.cov = zeros(gmm.dim,gmm.dim,gmm.mix);
for i = 1:gmm.mix
    gmm.mu(:,i) = g.m{i};
    gmm.cov(:,:,i) = g.ir{i}^-1;
    gmm.wt(i) = g.a(i);
end

aname = ['te'; 'tr'];
acolor = ['b'; 'r'; 'k'];
astyle = ['+'; 'x'; '*'];



% plot pdf contour
x1 = linspace(-3, 3, 50);
x2 = linspace(-3, 3, 50);
[X1,X2] = meshgrid(x1, x2);
num = size(X1(:), 1);
X = [X1(:) X2(:)]; % N * dim
V = zeros(num, 1);

% p = \sum c x N (x, u, cov)
for k = 1:num
    V(k) = mgdprob(gmm, X(k,:));
end
V = reshape(V, length(x2), length(x1));
p = contour(x1, x2, V); hold on; grid on;
xlabel('x1'); ylabel ('x2');
smp = 100;
mell = getmellipse(gmm, smp);
for j = 1:gmm.mix
    plot(mell(1,:,j), mell(2,:,j), 'g', 'LineWidth', 2); hold on;
    plot(gmm.mu(1,j), gmm.mu(2,j), 'g+', 'LineWidth', 2);
end
axis([-3 3 -3 3]);

