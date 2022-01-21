function [h1, h2, cbar] = FUN_quiver_by_plotV2_preset_cmap( x, y, u, v, vel_scale, zval, cmap, cmap_limit, varargin )
%  [h1,h2, cbar] = FUN_quiver_by_plotV2_preset_cmap( x, y, u, v, vel_scale, cmap, cmap_limit, varargin )
% Plot quiver with proset colormap
%
% *********************************************************************
% ***** This is replaced `FUN_quiver_by_plotV2_preset_cmap_patch` *****
% *********************************************************************
%
% Notice: Please note the the colormap should not be eidted after this function.
%         Currently, the color will not be updated with colormap.
%
%
% INPUT
%       x, y, u, v, vel_scale, varargin: same as FUN_quiver_by_plotV2
%       zval: value used to decide color
%       cmap [ N x 3 ]: colormap
%       cmap_limit [1 x 2]: limit for colormap, which is the input for caxis
%
%   INPUT-optional
%       interval: interval of quiver plotted
%                 e.g., if interval = 2, only the u(1:2:end, 1:2:end),
%                 v(1:2:end, 1:2:end) will be plotted.
%       {oters}: see FUN_quiver_by_plotV2
% OUTPUT
%       h1: array of handles for head body
%       h2: array of handles for head arrow
%       cbar: handle for colormap

% VEND: by L. Chi, This is replaced by `FUN_quiver_by_plotV2_preset_cmap_patch`. This function will not be updated any more
% V1.11 by L. Chi, 2021-04-01: replace `colormap(cmap)` by `set(gca,'colormap',cmap);`
% V1.10 By L. Chi, 2021-01-18: support a new input attribution: 'interval'
% V1.02 By L. Chi, 2020-11-01: fix a bug: An error may occur if all input
%                              values are NaN
% V1.01 By L. Chi, 2020-08-25: fix a bug
% V1.00 By L. Chi, 2020-08-25


% # load parameters =======================================================
 [interval, varargin] = FUN_codetools_read_from_varargin( varargin, 'interval', 1, true );
 [interval_x, varargin]       = FUN_codetools_read_from_varargin( varargin, 'interval_x', interval );
 [interval_y, varargin]       = FUN_codetools_read_from_varargin( varargin, 'interval_y', interval );

if isempty( zval )
    zval = sqrt( u.^2 + v.^2 );
end

% # Ini ===================================================================
if isvector( x ) && length(x) == size(u,1) && isvector(y) && length(y) == size(u,2) 
    
    [x,y]=meshgrid(x,y);
    x = x';
    y = y';
    
end

% ## reshape --------------------------------------------------------------

    % apply interval
    x = x(1:interval_x:end,1:interval_y:end);
    y = y(1:interval_x:end,1:interval_y:end);
    u = u(1:interval_x:end,1:interval_y:end);
    v = v(1:interval_x:end,1:interval_y:end);
    zval = zval(1:interval_x:end,1:interval_y:end);
    
    % reshape
    x = x(:).';
    y = y(:).';
    u = u(:).';
    v = v(:).';
    vel_scale = vel_scale(:)';
    zval = zval(:).';
    
% ## Remove nan -----------------------------------------------------------

    nanloc = isnan( x ) | isnan( y ) | isnan( u ) | isnan( v );
    x( nanloc ) = [];
    y( nanloc ) = [];
    u( nanloc ) = [];
    v( nanloc ) = [];
    if length( zval ) == length( nanloc )
       zval( nanloc ) = []; 
    end
    


% # Plot quivers ==========================================================
[h1,h2] = FUN_quiver_by_plotV2( x, y, u, v, vel_scale, 'is_plot_via_arrayfun', true, varargin{:});


nanloc = isnan( zval );
zval = zval(~nanloc);

if ~isempty(h1)
    h1 = h1(~nanloc);
end

if ~isempty(h2)
    h2 = h2(~nanloc);
end

% # prepare color for each  ==========================================================
cbar = colorbar;
caxis( cmap_limit );
%colormap(cmap); this is replaced by `set(gca,'colormap',cmap);` since the first one will apply to not only the current axes but all the other axes in the same figure
set(gca,'colormap',cmap);
%color_mat = cmapping( zval(~nanloc) , cmap );
color_mat = FUN_Plot_colormap_interp( cmap_limit, cmap, zval(~nanloc) );
% # Apply colormap ========================================================

if ~isempty(h1)
    arrayfun(@(x,y1,y2,y3)set(x,'color',[y1, y2, y3]),h1(:),color_mat(:,1),color_mat(:,2),color_mat(:,3));
end


if ~isempty(h2)
if isprop( h2(1), 'color')
    arrayfun(@(x,y1,y2,y3)set(x,'color',[y1, y2, y3]),h2(:),color_mat(:,1),color_mat(:,2),color_mat(:,3));
    
elseif isprop( h2(1), 'FaceColor')
    arrayfun(@(x,y1,y2,y3)set(x,'FaceColor',[y1, y2, y3]),h2(:),color_mat(:,1),color_mat(:,2),color_mat(:,3));
    
else
    error('Cannot find the property for color!');
    
end
end


