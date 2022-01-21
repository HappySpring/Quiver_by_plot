clear all
close all
clc

%% A simple case

x = 1:40;
y = 1:30;

[X,Y] = meshgrid(x,y);
X=X';
Y=Y';

u = sin( 2*pi/30 .* X );
v = cos( 2*pi/30 .* Y );


%%
figure

% scale for arrows.
%     0:         auto
%     otherwise: the scale will be applied to both u and v.
arrow_scale = 0; 

FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, arrow_scale );

grid on
box on

xlabel('x');
ylabel('y');

colorbar

title({'Example 2.1: a simple case', 'Color is defined by the magnitude of arrows'});

FUN_easy_export_fig('Example_2.1.jpg','-m2');


%%
value_for_color = sqrt(X.^2+Y.^2);

figure

arrow_scale = 0; 
FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, arrow_scale, 'zval', value_for_color );

grid on
box on

xlabel('x');
ylabel('y');

colorbar

title({'Example 2.2: Color is defined by parameter "zval"'});


FUN_easy_export_fig('Example_2.2_zval.jpg','-m2');


%% 


