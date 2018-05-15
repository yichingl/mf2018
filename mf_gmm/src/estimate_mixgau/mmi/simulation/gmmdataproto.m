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

function [flag] = gmmdata(bsimu, bplot, strMod, trnum, tenum, bRecord, bshow)

if nargin < 3
    printf('gmmdata(bsimu, bplot, strMod, trnum, tenum, bshow)');
    exit(-1);
end

% cd '/net/homes/wchu/projects/GMM/matlab/mgd';
strRoot = '/net/homes/wchu/projects/GMM/';
% strDir   - output data dir
% strLst   - scriptt data dir
strModDir    = strcat(strRoot, 'model/');
strDataDir   = strcat(strRoot, 'data/');
strScpDir    = strcat(strRoot, 'script/');
strModPath   = strcat(strModDir, strMod);
sprintf('Model name: %s\n', strModPath)

log2pi = log(2*pi);
LOGZERO = -1.0E20;
LSMALL = -0.5E20;
minLogExp = -log(-LOGZERO); % -20

fid = fopen(strModPath, 'rt');
if fid == -1
    s = sprintf('cannot open %s for read', strModPath);
end

iMod = fscanf(fid, '<MODNUM> %d\n', 1);
for i = 1:iMod
    strName = fscanf(fid, '<NAME> %s USER\n', 1);
    dim = fscanf(fid, '<VECSIZE> %d', 1);
    strCovType = fscanf(fid, '%s\n', 1);
    mix = fscanf(fid, '<NUMMIXES> %d\n', 1);
    
    wt  = zeros(1, mix);
    mu  = zeros(dim, mix);
    cov = zeros(dim, dim, mix);
    for j = 1:mix
        mixnum = fscanf(fid, '<MIXTURE> %d', 1);
        wt(j)  = fscanf(fid, '%f\n', 1);
        s = fscanf(fid, '%s\n', 1); % <MEAN>
        mu(:, j) = fscanf(fid, '%f', dim);
        s = fscanf(fid, '%s\n', 1); % <VARIANCE>
        if strcmp(strCovType, '<DIAG>') == 1
            for k = 1:dim
                cov(k, k, j) = fscanf(fid, '%f', 1);
            end
        else
            cov(:,:,j) = fscanf(fid, '%f', [dim dim]);
            cov(:,:,j) = cov(:,:,j)';
        end
        fscanf(fid, '\n', 1);
    end
    s = fscanf(fid, '%s\n', 1); % <ENDGMM>
    
    if abs(sum(wt) - 1) > 1e-3
        s = sprintf('sum of mix weight - %.3f does not equal 1.\n', sum(wt))
    end
end

for i = 1:iMod
    % t = 1 -> simu data for testing
    % t = 2 -> simu data for training

    % simulate data
    if bsimu == 1
        anum = [tenum; trnum];
        aname = ['te'; 'tr'];
        acolor = ['b'; 'r'; 'k'];
        astyle = ['*'; 'x'; '+'];

        rand('twister', 0);
        randn('state', 0);
        for t = 1:2
            y = mgdistri(anum(t), dim, mix, wt, mu, cov);

            %write down in natural order
            if bRecord == 1
                strDataPath = strcat(strDataDir, strName, '.', aname(t,:), num2str(anum(t)) );
                fo = fopen(strDataPath, 'wb');
                if fo == -1
                    printf('cannot open %s for write', strDataPath);
                end
                %write htk header
                fwrite(fo, anum(t), 'int');       % # of samples
                fwrite(fo, 100000, 'int');  % sampling period
                fwrite(fo, 4*dim, 'short');     % sample size
                fwrite(fo, 9, 'short');     % feature type - 9:USER
                %write data y - n x dim
                fwrite(fo, y, 'float');
                fclose(fo);

                strScpPath = strcat(strScpDir, strName, '.', aname(t,:), num2str(anum(t)) );
                if i == 1 && exist(strScpPath, 'file')
                    delete(strScpPath);
                end
                fo = fopen(strScpPath, 'at+');        
                if fo == -1
                    printf('cannot open %s for write', strScpPath);
                end
                if i == 1
                    fprintf(fo, '%d\n', iMod);
                end
                fprintf(fo, '%d\n', 1);
                fprintf(fo, '%s\n', strDataPath);
                fclose(fo);
            end
            
            %some plot cmd
            if bshow == 1 && t == 2
                %figure;
                pt = strcat(acolor(i), astyle(i));
                plot(y(1,:), y(2,:), pt); hold on;
                xlabel('x1'); ylabel ('x2');
                title(['the data of ', num2str(i), 'th model: ', strName]);
            end
        end
    end
    
    % plot pdf contour
    if bplot == 2
        x1 = linspace(-3, 3, 50);
        x2 = linspace(-3, 3, 50);
        [X1,X2] = meshgrid(x1, x2);
        num = size(X1(:));
        X = [X1(:) X2(:)]; % N * dim
        V = zeros(num(1), 1);
        
        % p = \sum c x N (x, u, cov)
        for k = 1:num
            v = LOGZERO;
            for j = 1:mix
                %V = mvnpdf(X, mu(:, j)', cov(:,:,j));
                c = dim * log2pi + log( det(cov(:,:,j)) );
                z = c + (X(k,:) - mu(:,j)')*inv(cov(:,:,j))*(X(k,:) - mu(:,j)')';
                z = z * -0.5 + log(wt(j));
                                
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
            V(k) = v;
        end
        V = reshape(V, length(x2), length(x1));
        contour(x1, x2, V); hold on; grid on;hold on;
        xlabel('x1'); ylabel ('x2');
        title(['the ', num2str(i), ' th model: ', strName]);
    %elseif bplot == 2
        for j = 1:mix
            sn = 100;
            t = linspace(0, 2*pi, sn);
            x = [cos(t); sin(t)];
            [L, p] = chol(cov(:,:,j), 'lower');
            if p ~= 0
                error('cov - covariance should be PSD (actually symmetric)')
            end
            y = mu(:,j)*ones(1,sn) + L*x;
            plot(y(1,:), y(2,:), 'b', 'LineWidth', 2); hold on;
            plot(mu(1,j), mu(2,j), 'r+', 'LineWidth', 2);
            %axis([-3 3 -3 3]);
        end
    end

end
