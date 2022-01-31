clear all
close all
clc


%% prepare data

x = 1:40;
y = 1:30;

[X,Y] = meshgrid(x,y);
X=X';
Y=Y';

u = sin( 2*pi/30 .* X );
v = cos( 2*pi/30 .* Y );
%% Example: FUN_quiver_by_plotV2
figure
xlim([1 41]);
ylim([1 31]);

arrow_scale = 0; % 0: the scale is determined by the script automatically.
FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, arrow_scale );

grid on
box on
xlabel('x');
ylabel('y');
colorbar

title({'Example 2.1: a simple case for colorful arrows', 'Color is defined by the magnitude of arrows'});

FUN_easy_export_fig('Demo_2.1.jpg','-m2');


%% Example: FUN_quiver_by_plotV2_mmap
figure
xlim([1 41]);
ylim([1 31]);


arrow_scale = 0; %auto

value_for_color = sqrt(X.^2+Y.^2);
FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, arrow_scale, 'zval', value_for_color );

grid on
box on

xlabel('x');
ylabel('y');

colorbar

title({'Example 2.2: Color is defined by parameter "zval"'});

FUN_easy_export_fig('Example_2.2.jpg','-m2');


%% Example: FUN_quiver_by_plotV2_cmap_patch_mmap

x1 = x/1.5-88;
y1 = y + 22;

figure
xlim([-88 -62]);
ylim([22 48]);

m_proj('lambert', 'lon', [-88 -62], 'lat', [22 48] );
FUN_quiver_by_plotV2_cmap_patch_mmap( x1, y1, u, v, 0, 'interval', 2 );
m_grid
m_coast('color','b')
colorbar

title('Example: FUN_quiver_by_plotV2_cmap_patch','Interpreter','none');

FUN_easy_export_fig('Example_2.3.jpg','-m2');

