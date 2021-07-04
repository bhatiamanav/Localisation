x = [2.5,2.0,1.5,1.0,2.5,2.9,0.5,0.8,2.9,2.0];
y = [2.9,2.0,1.5,0.5,0.5,1.0,1.5,3.5,2.5,3.9];
z = [0.5,0.8,1.0,0.6,1.0,1.2,1.0,1.5,1.0,1.6];

%Results from first order reflections
x_pred = [3.4,1.1,1.5,1.8,3.4,2.2,1.5,1.8,1.9,1.2];
y_pred = [1.8,2.8,2.2,1.3,2.6,0.2,1.0,4.2,1.9,2.01];
z_pred = [0.34,1.8,1.9,1.5,1.3,2.1,1.9,3.0,2.1,1.1];

x_pred_better = [2.8,1.6,1.7,1.3,2.2,1.3,1.1,1.3,2.1,1.51];
y_pred_better = [1.9,2.1,1.1,1.0,1.0,0.6,1.7,3.8,2.1,2.5];
z_pred_better = [0.45,1.7,1.4,1.2,1.3,1.4,1.43,2.4,1.55,1.4];

figure(1)
scatter3(x, y, z, 'r', 'filled')
hold on
scatter3(x_pred, y_pred, z_pred,'p', 'g', 'filled')
hold on
scatter3(x_pred_better, y_pred_better,z_pred_better,'^','b','filled')
hold off
grid on
legend('original','predicted','refined')

err_x = sum((bsxfun(@minus,x,x_pred)).^2);
err_y = sum((bsxfun(@minus,y,y_pred)).^2);
err_z = sum((bsxfun(@minus,y,y_pred)).^2);

err_x2 = sum((bsxfun(@minus,x,x_pred_better)).^2);
err_y2 = sum((bsxfun(@minus,y,y_pred_better)).^2);
err_z2 = sum((bsxfun(@minus,z,z_pred_better)).^2);

disp(sqrt(err_x+err_y+err_z));
disp(sqrt(err_x2+err_y2+err_z2));