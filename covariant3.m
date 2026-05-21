function [ data, nx, ny, nz ] = read_joint(filename)
fid = fopen(filename, 'rb');
dims = fread(fid, 3, 'int32');
nx = dims(1);
ny = dims(2);
nz = dims(3);
data = fread(fid, inf, 'single');
fclose(fid);
data = reshape(data, nx, ny, nz);
end

function data = read_marginal(filename)
fid = fopen(filename, 'rb');
points = fread(fid, 1, 'int32');
data = fread(fid, points, 'single');
fclose(fid);
end

function cmap = hueMap(n)
h = linspace(.6666, 0, n).';                 % column vector of hues
s = ones(n,1);                             % saturation = 1
v = 0.75 * ones(n,1);                      % value = 0.75
cmap = hsv2rgb([h s v]);                   % convert to RGB
end

% Generate a colormap for visualization
nColors = 256; % Number of colors in the colormap
cmap = hueMap(nColors);
qmap = hueMap(10);

[ quantile, nx, ny, nz ] = read_joint('QC.bin');
[ classes, nx, ny, nz ] = read_joint('classes.bin');

x = linspace(0, 1, nx);
y = linspace(0, 1, ny);
z = linspace(0, 1, nz);
[X, Y, Z] = meshgrid(x, y, z);

% Force displayed volume to be unit cube
ax = gca;
ax.XLim = [0 1];
ax.YLim = [0 1];
ax.ZLim = [0 1];
t = ax.CameraTarget;           % camera target
v = ax.CameraPosition - t;     % view vector from target to camera
%ax.CameraPosition = t + 2*v;   % move camera twice as far from target

% Make data units equal in all directions and prevent autoscaling on rotate
daspect([1 1 1])      % equal data-unit aspect ratio
axis vis3d           % freeze aspect ratio and limits during rotation
xlabel('x'); ylabel('y'); zlabel('z');

%figure;
% p = isosurface(X, Y, Z, quantile, 0.1);
% patch(p, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% p = isosurface(X, Y, Z, quantile, 0.3);
% patch(p, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% p = isosurface(X, Y, Z, quantile, 0.5);
% patch(p, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% p = isosurface(X, Y, Z, quantile, 0.7);
% patch(p, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% p = isosurface(X, Y, Z, quantile, 0.9);
% patch(p, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% figure;
p = isosurface(X, Y, Z, quantile, 0.5);
patch(p, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% % p = isosurface(X, Y, Z, classes, .5);
% % patch(p, 'FaceColor', 'blue');
% % % axis tight;
% % daspect([1 1 1]);
view(3);
grid(ax,"on");
% daspect([1 1 1])      % equal data-unit aspect ratio
% axis vis3d           % freeze aspect ratio and limits during rotation
% xlabel('x'); ylabel('y'); zlabel('z');
