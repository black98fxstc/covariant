%% 
%% 
%% 

function [ data, nx, ny ] = read_joint(filename)
fid = fopen(filename, 'rb');
dims = fread(fid, 2, 'int32');
nx = dims(1);
ny = dims(2);
data = fread(fid, inf, 'single');
fclose(fid);
data = permute(reshape(data, nx, ny), [2 1]);
end

function data = read_marginal(filename)
fid = fopen(filename, 'rb');
points = fread(fid, 1, 'int32');
data = fread(fid, points, 'single');
fclose(fid);
end

[ F, nx, ny ] = read_joint('f.bin');
[ QC, nx, ny ] = read_joint('QC.bin');
[ W, nx, ny ] = read_joint('w.bin');
[ f1, nx, ny ] = read_joint('f1.bin');
[ f2, nx, ny ] = read_joint('f2.bin');
[ S1, nx, ny ] = read_joint('S1.bin');
[ S2, nx, ny ] = read_joint('S2.bin');
[ T1, nx, ny ] = read_joint('T1.bin');
[ T2, nx, ny ] = read_joint('T2.bin');
[ R, nx, ny ] = read_joint('R.bin');

P1 = read_marginal('P1.bin');
P2 = read_marginal('P2.bin');
Q1 = read_marginal('Q1.bin');
Q2 = read_marginal('Q2.bin');

% Create grid coordinates
x = linspace(0, 1, nx);
y = linspace(0, 1, ny);
[X, Y] = meshgrid(x, y);

% Plot
figure;
surf(X, Y, W);
shading interp;
view(2);
colorbar;
axis([0 1 0 1]);
axis square;
title('Weights w(x, y)');

figure;
surf(X, Y, F);
shading interp;
view(2); % Top-down view
colorbar;
axis([0 1 0 1]);
axis square;
xlabel('x');
ylabel('y');
title('Function f(x, y)');

figure;
hold on;
contour(X, Y, QC, 0.1:0.1:0.9, 'k');
contour(X, Y, QC, [0.001 0.001], 'b');
% contour(X, Y, S1', [0 0], 'b', 'LineWidth', 2);
% contour(X, Y, S2', [0 0], 'g', 'LineWidth', 2);
% contour(X, Y, L', [0 0], 'g', 'LineWidth', 2);
axis([0 1 0 1]);
axis square;
title('Quantile Contours');
hold off

figure;
surf(X, Y, f1);
shading interp;
view(2);
colorbar;
axis([0 1 0 1]);
axis square;
title('f1(x, y)');

figure;
surf(X, Y, f2);
shading interp;
view(2);
colorbar;
axis([0 1 0 1]);
axis square;
title('f2(x, y)');

figure;
surf(X, Y, S1);
shading interp;
view(2);
colorbar;
axis([0 1 0 1]);
axis square;
title('S1(x, y)');

figure;
surf(X, Y, S2);
shading interp;
view(2);
colorbar;
axis([0 1 0 1]);
axis square;
title('S2(x, y)');

figure;
surf(X, Y, T1);
shading interp;
view(2);
colorbar;
axis([0 1 0 1]);
axis square;
title('T1(x, y)');

figure;
surf(X, Y, T2);
shading interp;
view(2);
colorbar;
axis([0 1 0 1]);
axis square;
title('T2(x, y)');

figure;
surf(X, Y, R);
shading interp;
view(2);
colorbar;
axis([0 1 0 1]);
axis square;
title('R(x, y)');

% figure;
% surf(X, Y, L);
% shading interp;
% view(2);
% colorbar;
% axis([0 1 0 1]);
% axis square;
% title('Laplacian(x, y)');
% figure;
% surf(X, Y, S1');
% shading interp;
% view(2);
% hold on;
% contour(X, Y, P', 0.1:0.1:0.9, 'k');
% hold off;
% colorbar;
% axis([0 1 0 1]);
% axis square;
% title('S1(x, y) surface with P contours');
