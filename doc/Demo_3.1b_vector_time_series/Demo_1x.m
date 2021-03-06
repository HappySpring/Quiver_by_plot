clear all
close all
clc


%% prepare data

x = 0:5:270;
u = cosd(x);
v = sind(x);

%% Example 1.1 A simple case

figure('position',[  100 100  980  300])
hold on

% please always set xlim and ylim first
xlim([min(x) max(x)]);
ylim([-1.2 1.2])
vel_plot_scale = 1; % arrow sacle

% plot arrows
% 'is_correct_angle' is set to true to plot arrows according to its real
%   direction, ignorning the ratio between x & y axis. 
FUN_quiver_by_plotV2( x, zeros(size(x)), u, v, vel_plot_scale, 'is_correct_angle', true, 'is_plot_head', false);
box on

plot( x, zeros(size(x)), '-b');

title('Example: Vector time series');

% This is private function, please use saveas
FUN_easy_export_fig('Example_1.1.jpg','-m2');

