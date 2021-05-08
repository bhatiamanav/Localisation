x = [2.5,2.0,1.5,1.0,2.5,2.9,0.5,0.8,2.9,2.0];
y = [2.9,2.0,1.5,0.5,0.5,1.0,1.5,3.5,2.5,3.9];
z = [0.5,0.8,1.0,0.6,1.0,1.2,1.0,1.5,1.0,1.6];

x_pred = [3.5,1.0,1.5,1.9,3.6,1.9,1.4,1.8,2.1,1.0];
y_pred = [1.9,3.0,2.4,1.5,2.5,0.1,1.1,4.5,1.8,1.91];
z_pred = [0.3,1.9,2.0,1.6,1.4,2.2,2.0,3.5,2.5,0.9];

x_pred_better = [3.0,1.5,1.8,1.5,2.0,2.4,0.8,1.3,2.3,1.5];
y_pred_better = [2.3,2.5,1.1,1.0,1.0,0.6,1.7,3.8,2.1,2.5];
z_pred_better = [0.4,1.8,1.5,1.1,1.2,1.6,1.5,2.6,1.7,1.3];

figure(1)
scatter3(x, y, z, 'r', 'filled')
hold on
scatter3(x_pred, y_pred, z_pred, 'y', 'filled')
hold on
scatter3(x_pred_better, y_pred_better,z_pred_better,'b','filled')
hold off
grid on
legend('original','predicted','refined')