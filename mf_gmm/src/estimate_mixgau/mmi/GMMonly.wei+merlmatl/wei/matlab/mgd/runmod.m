all = 1;
r = 2;
l = 4;
%figure;
for i = 1:all
    %subplot(r, l, i);hold on;
    %name = ['g.n1.d2.m', num2str(i), '.f.tr1000.def'];
    name = 'g.n2.d2.m2.f.smu';
    gmmdata(1, 2, name, 10000, 1000, 1, 1);
end