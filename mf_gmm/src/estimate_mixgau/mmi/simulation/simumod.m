
%simu
num = 100;
while num <= 100000
    figure;
    name = 'g.n2.d2.m4.f';
    gmmdata(1, 1, name, num, 1000, 1, 1);
    num = num *10;
end
