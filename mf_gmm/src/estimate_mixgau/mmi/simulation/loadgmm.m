function [gmm, iMod] = loadgmm(strModPath)
%[gmm, iMod] = loadgmm(strModPath)
% Read in Model

if nargin ~= 1
    printf('loadgmm(strModPath)');
    exit(-1);
end

fid = fopen(strModPath, 'rt');
if fid == -1
    s = sprintf('cannot open %s for read', strModPath);
end

%iMod = fscanf(fid, '<MODNUM> %d\n', 1);
iMod = fscanf(fid, '%d', 1);
gmm(iMod).name = 0;
gmm(iMod).dim = 0;
for i = 1:iMod
    gmm(i).name = ['model_' num2str(i)];
    tmp = fscanf(fid,'%d',2);
    gmm(i).mix = tmp(1);
    gmm(i).dim = tmp(2);
    gmm(i).ctype = '<DIAG>';
    
    gmm(i).wt  = zeros(1, gmm(i).mix);
    gmm(i).mu  = zeros(gmm(i).dim, gmm(i).mix);
    gmm(i).cov = zeros(gmm(i).dim, gmm(i).dim, gmm(i).mix);
    for j = 1:gmm(i).mix
        gmm(i).wt(j)  = fscanf(fid, '%f\n', 1);
        gmm(i).mu(:, j) = fscanf(fid, '%f', gmm(i).dim);
        if strcmp(gmm(i).ctype, '<DIAG>') == 1
            for k = 1:gmm(i).dim
                gmm(i).cov(k, k, j) = fscanf(fid, '%f', 1);
            end
        else
            gmm(i).cov(:,:,j) = fscanf(fid, '%f', [gmm(i).dim gmm(i).dim]);
            gmm(i).cov(:,:,j) = gmm(i).cov(:,:,j)';
        end
    end
    if abs(sum(gmm(i).wt) - 1) > 0.03
        s = sprintf('sum of mix weight - %.3f does not equal 1.\n', sum(gmm(i).wt))
    end
end