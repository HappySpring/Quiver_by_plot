clear all
close all
clc


%% prepare data

x = -1:0.2:1;
y = -1:0.2:1;

[X,Y] = meshgrid(x,y);
X=X';
Y=Y';

theta = atan2d(Y,X);

r = sqrt( X.^2 + Y.^2);

u =  sind( theta ) .* r / 100;
v = -cosd( theta ) .* r / 100;

x = x*10 - 75;
y = y*10 + 35;


%%
figure
hold on

m_proj('lambert', 'lon', [-88 -62], 'lat', [22 48] );
m_grid

FUN_quiver_by_plotV2_cmap_patch_mmap( x, y, u, v, 0 );
m_coast('color','b')
colorbar


FUN_easy_export_fig('logo.jpg','-m2');

