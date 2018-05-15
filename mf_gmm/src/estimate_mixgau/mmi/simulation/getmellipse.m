% get the data for drawing the equal prob curve (the eclipse on 1 std)
% of each mixture of the multi Gaussians

function [y] = getmellipse(gstr, smp)

t = linspace(0, 2*pi, smp);
x = [cos(t); sin(t)];
y = zeros(gstr.dim, smp, gstr.mix);
for j = 1:gstr.mix
    [L, p] = chol(gstr.cov(:,:,j), 'lower');
    if p ~= 0
        error('cov - covariance should be PSD')
    end
    y(:,:,j) = gstr.mu(:,j)*ones(1, smp) + L*x;
end