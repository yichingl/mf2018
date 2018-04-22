boat_feat = load('cnn__boat_v55hojInS9c_011.csv');
heli_feat = load('cnn__heli_72qNna11Znw_045.csv');
N = length(boat_feat);
figure;
subplot(2,1,1);
plot(1:N,boat_feat);
title('CNN Proxy Features (Boat)');
xlabel('node');
ylabel('weight');

subplot(2,1,2);
plot(1:N,heli_feat);
title('CNN Proxy Features (Heli)');
xlabel('node');
ylabel('weight');