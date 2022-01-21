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


%% Example: FUN_quiver_by_plotV2

figure
xlim([-88 -62]);
ylim([22 48]);

FUN_quiver_by_plotV2( x, y, u, v, 0 )

title('Example 3.1: FUN_quiver_by_plotV2','Interpreter','none');

FUN_easy_export_fig('Example_3.1.jpg','-m2');


% % %% Example: FUN_quiver_by_plotV2_mmap
% % figure
% % xlim([-88 -62]);
% % ylim([22 48]);
% % 
% % m_proj('lambert', 'lon', [-88 -62], 'lat', [22 48] );
% % FUN_quiver_by_plotV2_mmap( x, y, u, v, 0)
% % m_grid
% % m_coast('color','b')
% % 
% % title('Example: FUN_quiver_by_plotV2_mmap','Interpreter','none');
% % 
% % FUN_easy_export_fig('Example_1.2b.jpg','-m2');
% % 
% % 
% % %% Example: FUN_quiver_by_plotV2_cmap_patch
% % 
% % figure
% % xlim([-88 -62]);
% % ylim([22 48]);
% % 
% % FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, 0 )
% % colorbar
% % title('Example: FUN_quiver_by_plotV2_cmap_patch','Interpreter','none');
% % 
% % FUN_easy_export_fig('Example_1.2c.jpg','-m2');
% % 
% % 
% % %% Example: FUN_quiver_by_plotV2_cmap_patch_mmap
% % 
% % figure
% % xlim([-88 -62]);
% % ylim([22 48]);
% % 
% % m_proj('lambert', 'lon', [-88 -62], 'lat', [22 48] );
% % FUN_quiver_by_plotV2_cmap_patch_mmap( x, y, u, v, 0 )
% % m_grid
% % m_coast('color','b')
% % colorbar
% % 
% % title('Example: FUN_quiver_by_plotV2_cmap_patch','Interpreter','none');
% % 
% % FUN_easy_export_fig('Example_1.2d.jpg','-m2');
% % 
