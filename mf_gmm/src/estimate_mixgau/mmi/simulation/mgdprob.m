% p = \sum c x N (x, u, cov)
% X - 1 x dim

function [v] = mgdprob(gstr, X)

log2pi = log(2*pi);
LOGZERO = -1.0E20;
LSMALL = -0.5E20;
minLogExp = -log(-LOGZERO); % -20

if size(X,2) ~= gstr.dim
    s = ['feature dimesion', num2str(size(X,2)), not equal, 'model dimension', num2str(gstr.dim)]
    return;
end
v = LOGZERO;
for j = 1:gstr.mix
    %V = mvnpdf(X, mu(:, j)', cov(:,:,j));
    c = gstr.dim * log2pi + log( det(gstr.cov(:,:,j)) );
    z = c + (X - gstr.mu(:,j)')*inv(gstr.cov(:,:,j))*(X - gstr.mu(:,j)')';
    z = z * -0.5 + log(gstr.wt(j));

    % logv = log(v + z)
    if z < v
        temp = z; z = v; v = temp;
    end
    diff = v - z;
    if diff < minLogExp
        if z < LSMALL
            v = LOGZERO;
        else
            v = z;
        end
    else
        v = z + log(1.0 + exp(diff));
    end
end
