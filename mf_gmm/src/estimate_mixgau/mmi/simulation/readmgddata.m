% Read Data

function [y] = readmgddata(str)
fo = fopen(str, 'rb');
if fo == -1
    printf('cannot open %s for read', str);
end
%write htk header
num = fread(fo, 1, 'int');       % # of samples
smp = fread(fo, 1, 'int');  % sampling period
dim = fread(fo, 1, 'short')/4;     % sample size, 4 - sizeof float
ft  = fread(fo, 1, 'short');     % feature type - 9:USER
%write data y - n x dim
y = zeros(num, dim);
for i = 1:num
    y(i,:) = fread(fo, dim, 'float');
end
fclose(fo);

figure;
plot(y(:,1), y(:,2), '*');
