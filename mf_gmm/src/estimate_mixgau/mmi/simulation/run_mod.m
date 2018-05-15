for i = 1:8
    figure;
    col = 4;
    r = mod(i,col);
    l = i-(r-1)*col;
    subplot(r, l, i);
    name = ['g.n1.d2.m', num2str(i), '.f.tr100000.def'];
   %name = 'g.n1.d2.m4.f.smu';
    gmmdata(1, 2, name, 1000, 100, 0, 0);
end