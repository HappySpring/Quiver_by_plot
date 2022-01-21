function [h1,h2, uu, vv, hu, hv] = FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, vel_scale, varargin )
% [h1,h2, uu, vv, hu, hv] = FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, vel_scale, varargin )
% [h1,h2, uu, vv, hu, hv] = FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, vel_scale, 'parameter','value',...)
% 
% Plot quivers with "patch" functions. This make it easier to control the final results.
% 
% **Attention**************************************************************
%    Please set xlim and ylim manually before calling this function and make sure that the 
%    values set for xlim and ylim will not be changed later
% *************************************************************************
%
% ** INSTASLL **
%     This function requires another matlab function "plotboxpos", which is available
%     from https://github.com/kakearney/plotboxpos-pkg
%
% =========================================================================
% INPUT:
%       x, y, u, v : input data
%       vel_scale [default: 1 ]
%                  0: auto scale,  
%                  1: No scale, the original u/v will be used
%                  otherwise: this scale will be applied to u and v.
%           
%   Optional parameters and their default values:
%
%       zval: value used to decide color. (default: sqrt(u.^2+v.^2) )
%
%       is_plot_head: [default: true ]
%                  Plot arrow head, or not
%       is_plot_body: [default: true ]
%                  Plot arrow body, which is the part except arrow head.
%       fill_head: [default: false]
%                  If this is true, a filled triangle will be used as the head of each arrow.
%                  In this case, patch, instead of plot, is used to plot arrow heads.
%
%       head_length: [default: 0]
%                   It shares the same unit with u and v. (Thus, the vel_scale 
%                            will also be applied to head_length)
%                    0: it is set according to the gca (width & length in pixel).
%                  ~=0: it is the absolute length of the arrow head in **pixel**
%
%       head_angle: half of arrow head's open angle (in degree)  [default: 10]
%                  half angle of the arrow head.
%
%       is_correct_angle [default: false]
%   
%   Remainings in varargin will be used as paramters for `plot`
%
% (NOTE: The meaning of parameters plotarrows, head_length and head_angleh are different 
%  from the ones in FUN_quiver_by_plot (or the V1) )
%
% =========================================================================
% OUTPUT:
%       h1: handle array for arrow body
%       h2: handle array for arrow head
%       
%       uu [ 3 x N ]: x coordinate for arrow body, the 3 values in each col are [ x_origin; x_head; Nan]
%       vv [ 3 x N ]: y coordinate for arrow body, the 3 values in each col are [ y_origin; y_head; Nan]
%       hu [ 4 x N ]: x coordinate for arrow head, the 4 values in each col are [ x_head_edge, x_head, x_the_other_edge_of_head; Nan]
%       hv [ 4 x N ]: y coordinate for arrow head, the 4 values in each col are [ y_head_edge, y_head, y_the_other_edge_of_head; Nan]
%            you can re-plot the arrows by uu,vv hu, hv:
%
%                     plot( uu(:), vv(:) ); % plot arrow body
%                     hold on
%                     plot( hu(:), hv(:) ); % plot arrow head
%
% =========================================================================
% Example:
%       x = 0;
%       y = 0;
%       u = 1;
%       v = 2;
%       vel_scale = 1;
%
%       h = FUN_quiver_by_plotV2( x, y, u, v, vel_scale )
%       h = FUN_quiver_by_plotV2( x, y, u, v, vel_scale, 'is_plot_body', true,'head_length', 0, 'head_angle', 15, 'is_correct_angle', false )
% =========================================================================


% =========================================================================
% V2.73 by L. Chi
%         + Add "interval"
% V2.72 by L. Chi
%         + update default head_length following V2.66 in FUN_quiver_by_plotV2.m
% V2.71 by L. Chi
%         + fix a bug in estimating arrow head length in some cases
%         + add more examples (see the end of this file)
% V2.70 by L. Chi (2021/06/18)
%         + Replace `plot` by `patch` to support colormap
%         + arrayfun is not necessary anymore
%         + update methods for estimating xlimit and ylimit if an axes
%         has not be created.
%         + This replaces `FUN_quiver_by_plotV2_preset_cmap`
%
% ------- This is edited based on FUN_quiver_by_plotV2 (V2.64) ------------
% V2.64 by L. Chi (03/31/2021)
%         - support x,y as [1x1] matrix.
% V2.63 By L. Chi (02/xx/2020)
%         - add output uu, vv, hu, hv.
% V2.62 By L. Chi (01/21/2020) 
%         - fix bugs in comments
%         - improve codes format (no effect on results)   
% V2.61 By L. Chi (10/30/2020)
%         Fix a (level:mid) bug in calculating p_ratio_y_x: 
%            if the axis is set to a mode other than normal, the gca position
%            from "get(gca,'position')" is not accurate. A thrid-part function 
%            "plotboxpos" is used to calculate the accurate position. 
%            See https://github.com/kakearney/plotboxpos-pkg for details.
% V2.60 By L. Chi (10/30/2020)
%         Add support to fill the arrow head by line_color 
%           This feature can be turned on by set optional parameter fill_head as true.
%         Edit comments slightly.
%         Fix a bug: an hold on is added in the function to make sure all plots are hold.
% V2.52 By L. Chi (10/08/2020)
%         Fix a bug: the positions of arrows may shift unexpectedly when their body is not plotted. 
%         More details can be found at https://gitlab.com/rs16/peach/P01_HF_Reproduce_Ben/-/issues/23
% V2.51 By L. Chi (08/28/2020):
%         Improve compatibility to matlab R2016b (linux).
% V2.50 by L. Chi (08/26/2020):
%         Use colorful arrows.
%         Stations with valid data will be noted at the northwest corner.
%         Other tunes
% V2.21 by L. Chi (08/25/2020)
%       It is possible to use arrayfun in plotting. This will return one handle 
%       for each arrow, making it possbile to edit each arrow (e.g., color) separately
% V2.20 by L. Chi
%       Input parameters are now given in {'param','value'} mode.
% V2.10 by L. Chi
%       Update the methods for calculating head length automatically.
% V2.00 by L. Chi
%
%       The old version (FUN_quiver_by_plot.m) is not replaced by this one
%           to keep compatiblility to some old codes
% V1.00 
% =========================================================================

%%
% =========================================================================
% # default values
% =========================================================================
%                                                                    |input  | paramter       | default value | 
    [zval, varargin] = FUN_codetools_read_from_varargin(             varargin, 'zval',         sqrt(u.^2+v.^2) ); % zval defines color 

    [is_plot_head, varargin]     = FUN_codetools_read_from_varargin( varargin, 'is_plot_head', true ); % show head of arrows or not
    [is_plot_body, varargin]     = FUN_codetools_read_from_varargin( varargin, 'is_plot_body', true ); % show body of arrows or not
    [fill_head,    varargin]     = FUN_codetools_read_from_varargin( varargin, 'fill_head',    false );% fill heead of arrows or not

    [head_length, varargin]      = FUN_codetools_read_from_varargin( varargin, 'head_length',  0 );    % length of head.  
    [head_length_min, varargin]  = FUN_codetools_read_from_varargin( varargin, 'head_length_min', 0 ); % minimal length of head
    [head_angle, varargin]       = FUN_codetools_read_from_varargin( varargin, 'head_angle',   15 );   
    [is_correct_angle, varargin] = FUN_codetools_read_from_varargin( varargin, 'is_correct_angle', false );
    [is_correct_angle_by_edit_v, varargin] = FUN_codetools_read_from_varargin( varargin, 'is_correct_angle_by_edit_v', true );

    [interval,   varargin]       = FUN_codetools_read_from_varargin( varargin, 'interval', 1 );
    [interval_x, varargin]       = FUN_codetools_read_from_varargin( varargin, 'interval_x', interval );
    [interval_y, varargin]       = FUN_codetools_read_from_varargin( varargin, 'interval_y', interval );

    %[is_plot_arrow_body, varargin] = FUN_codetools_read_from_varargin( varargin, 'is_plot_arrow_body', true );
    %[is_plot_via_arrayfun, varargin] = FUN_codetools_read_from_varargin( varargin, 'is_plot_via_arrayfun', false );

    % vel_scale must be given explicitly
    if isempty( vel_scale )
        vel_scale = 1; 
    end
    
    % check
    if isnumeric( vel_scale )
        %PASS
    else
        error('vel_scale must be a number!')
    end
    
%%
% =========================================================================
% # initialization
% =========================================================================
%     set(gca,'xlimMode','manual')
%     set(gca,'ylimMode','manual')

    if length(x) == 1 && length(u)>1
       x = ones(size(u)).*x; 
    end
    
    if length(y) == 1 && length(u)>1
       y = ones(size(u)).*y; 
    end

    if isvector( x ) && length(x) == size(u,1) && isvector(y) && length(y) == size(u,2) 

        [x,y]=meshgrid(x,y);
        x = x';
        y = y';

    end
    
    if interval_x ~= 1 || interval_y ~= 1
        x = x(1:interval_x:end,1:interval_y:end);
        y = y(1:interval_x:end,1:interval_y:end);
        u = u(1:interval_x:end,1:interval_y:end);
        v = v(1:interval_x:end,1:interval_y:end);
        zval = zval(1:interval_x:end,1:interval_y:end);
    end

% ---- ## Make velocity vectors -------------------------------------------
    %size0 = size(x);
    x = x(:).'; 
    y = y(:).';
    u = u(:).'; 
    v = v(:).';
    zval = zval(:).';
    
% ---- ## Check -----------------------------------------------------------
    plot_size = size( u );

    if     isequal( size( x ), plot_size ) ...
        && isequal( size( y ), plot_size ) ...
        && isequal( size( u ), plot_size ) ...
        && isequal( size( v ), plot_size ) ...
        && isequal( size( zval ), plot_size )
        % Pass: x, y, u, v must share same size
    else
        error('The input variable x/y/u/v should share same size')
    end
    
% ---- ## Remove NaNs -----------------------------------------------------
    nanloc = isnan( x ) | isnan( y ) | isnan( u ) | isnan( v );
    x( nanloc ) = [];
    y( nanloc ) = [];
    u( nanloc ) = [];
    v( nanloc ) = [];
    zval(nanloc)= [];

% ---- ## return empty for empty results ----------------------------------
    if isempty( u ) || all( isnan(u) )
        h1 = [];
        h2 = [];
        uu = [];
        vv = [];
        hu = [];
        hv = [];
        return
    end

% ---- ## get figure and gca info -----------------------------------------
    % ### create axis if it does not exist
    tem = get(gcf,'children');

    if isempty(tem)
       % set xlimit/ylimit automatically

       if vel_scale == 0
           t_vel_scale = 1;
       else
           t_vel_scale  = vel_scale ;
       end
       tx = [ x; x+u*t_vel_scale*1.2 ];
       ty = [ y; y+v*t_vel_scale*1.2 ];
       xlimit = [ nanmin(tx(:)) nanmax(tx(:)) ];
       ylimit = [ nanmin(ty(:)) nanmax(ty(:)) ];

    else
       % read xlimit/ylimit from the current figure
       xlimit = get(gca,'xlim');
       ylimit = get(gca,'ylim');
    end
    
    if xlimit(1) == xlimit(2)
        xlimit = [xlimit(1) xlimit(1)];
    end
    
    if ylimit(1) == ylimit(2)
        ylimit = [ylimit(1) ylimit(1)];
    end   

   % ### check axis position/size
   gcf_pos_now = get( gcf, 'position');
   
   % gca_pos_now = get( gca, 'position');
   % the gca position from "get( gca, 'position' )" is not accurate in certain cases. 
   % see https://github.com/kakearney/plotboxpos-pkg for details.
   gca_pos_now = plotboxpos( gca );
   
   pix_x = gcf_pos_now(3) * gca_pos_now(3);
   pix_y = gcf_pos_now(4) * gca_pos_now(4);
   
   xlim_diff = xlimit(2) - xlimit(1);
   ylim_diff = ylimit(2) - ylimit(1);
   
   % ### prepare for angle correction
   pix_per_x = pix_x ./ xlim_diff;
   pix_per_y = pix_y ./ ylim_diff;
   
   p_ratio_y_x = pix_per_x ./ pix_per_y; % This is "pix_per_x ./ pix_per_y", not "pix_per_y ./ pix_per_x"
   
   if is_correct_angle 
       r_correct_v    = p_ratio_y_x;
       r_correct_head = p_ratio_y_x;
       
   else
       r_correct_v    = 1 ; % No correction for arrows. 
       r_correct_head = p_ratio_y_x; % The head should still be corrected.       
   end
%%
% =========================================================================
% # Plot velocity
% =========================================================================
    
% ---- ## preparation - save current settings -----------------------------
    cax  = gca;
    next = lower(get(cax,'NextPlot'));
    hold_state = ishold(cax);
    
    hold on

% ---- ## Velocity Scale --------------------------------------------------

    if vel_scale == 0
        % Base vel_scale value on average spacing in the x and y
        % directions.  Estimate number of points in each direction as
        % either the size of the input arrays or the effective square
        % spacing if x and y are vectors.
        if min(size(x))==1
            n=sqrt(numel(x)); 
            m=n;
        else
            [m,n]=size(x); 
        end
        
        if all( size(x) == 1 )
            delx = xlim_diff/3;
            dely = ylim_diff/3;
        else
            delx = diff([min(x(:)) max(x(:))])/n;
            dely = diff([min(y(:)) max(y(:))])/m;
        end
        del = delx.^2 + dely.^2;
        if del>0
            len = sqrt((u.^2 + v.^2)/del);
            len2 = median(len(:));
        else
            len2 = 0;
        end

        if len2>0
            vel_scale = 0.9 / len2 ;
        else
            vel_scale = 1;
        end
    end
    
    u = u*vel_scale;
    v = v*vel_scale;
    
% ---- ## Ploit quiver body -----------------------------------------------
    if is_plot_body
        if is_correct_angle_by_edit_v
            up = u;
            vp = v.*r_correct_v;
        else
            up = u./r_correct_v;
            vp = v;
        end
    else 
        if  isempty( head_length ) || head_length == 0
            % head_length = 1;
            warning('Please set head_length manually when is_plot_arrow_body = false')
        end
        up = zeros(size(u));
        vp = zeros(size(v));
    end

    if is_plot_body
        %uu = [x ; x+up; nan(size(u))];
        %vv = [y ; y+vp; nan(size(v))];
        uu = [ x; x+up ];
        vv = [ y; y+vp ];
        h1 = patch( 'XData', uu, 'YData', vv, 'CData', repmat(zval(:).',2,1), 'FaceColor','none', 'edgecolor','flat', varargin{:} );
        
    else
        uu = [];
        vv = [];
        h1 = [];
    end
%%
% =========================================================================
% # Plot quiver head
% =========================================================================

    if is_plot_head

        % ### arrow angle on the screen -----------------------------------
        if is_correct_angle
            tem_ang1   = atan2d( v, u );
        else
            tem_ang1   = atan2d( v/r_correct_head, u );
        end


        % ### Calculate the length of arrow head automatically ------------
        if  isempty( head_length ) || head_length == 0
            
            r_head_L = 0.4 ; % 40% of the arrow length
            
            abs_velp = sqrt( up.^2 + vp.^2 ) ;             
            head_length = abs_velp .* r_head_L ;
            
            % velocity scale has been applied to u and v (also up and vp).
            % Thus, it is should not applied used here.
            % % head_length = head_length * vel_scale;
            
        elseif head_length < 0 && head_length >= -1
            
            r_head_L = abs(head_length); 
            abs_velp = sqrt( up.^2 + vp.^2 ) ;             
            head_length = abs_velp .* r_head_L ;
            
        elseif head_length < -1
            error('Unexpected head_length!');
            
        else
            head_length = head_length * vel_scale;
        end
        
        % set min head length 
        if head_length_min > 0 
            head_length_min = head_length_min * vel_scale;
            head_length( head_length < head_length_min ) = head_length_min;            
        end

        % ### Generate arrays for arrow heads -----------------------------
        if is_plot_body
            % head is placed at the end of the arrow
            hu = [x+up-head_length.*cosd( tem_ang1+head_angle );  x+up; ...
                  x+up-head_length.*cosd( tem_ang1-head_angle ); ];

            hv = [y+vp-head_length.*sind( tem_ang1+head_angle )*r_correct_head;  y+vp; ...
                  y+vp-head_length.*sind( tem_ang1-head_angle )*r_correct_head;  ];
        else
            % the middle point of the head triangle is placed at (x,y)
            %
            % up, vp are kept in calculating tem_hu, tem_hv to keep a correct size.
            tem_hu = [up-head_length.*cosd( tem_ang1+head_angle );  up; ...
                      up-head_length.*cosd( tem_ang1-head_angle );  ];

            tem_hv = [vp-head_length.*sind( tem_ang1+head_angle )*r_correct_head;  vp; ...
                      vp-head_length.*sind( tem_ang1-head_angle )*r_correct_head; ];
               
            hu = x + tem_hu - nanmean(tem_hu);
            hv = y + tem_hv - nanmean(tem_hv);
        end
        
        % ### Plot head of arrows
        hold(cax,'on')
        
        % ### plot
        if fill_head
            h2 = patch( 'XData', hu, 'YData', hv, 'CData', repmat(zval(:).',3,1), 'FaceColor','flat', 'edgecolor','none', varargin{:} );
        else
            h21 = patch( 'XData', hu(1:2,:), 'YData', hv(1:2,:), 'CData', repmat(zval(:).',2,1), 'FaceColor','none', 'edgecolor','flat', varargin{:} );
            h22 = patch( 'XData', hu(2:3,:), 'YData', hv(2:3,:), 'CData', repmat(zval(:).',2,1), 'FaceColor','none', 'edgecolor','flat', varargin{:} );
            h2 = [h21; h22];
        end

        
    else
        hu = [];
        hv = [];
        h2 = [];
    end
    
%%
% =========================================================================
% # Postprocess
% =========================================================================

    if ~hold_state
        hold(cax,'off');
        view(cax,2); 
        set(cax,'NextPlot',next); 
    end

    %h = [h1;h2];

    % keep the original xlimit
    % The angle is corrected based on fixed the xlimit and ylimit.
    xlim(xlimit);
    ylim(ylimit);

%%
% =========================================================================
% # return
% =========================================================================

%     head_x = hu(2,:);
%     head_y = hv(2,:);

% if any( nanloc );
%     
%     if ~isempty( h1 )
%         h1_0 = nan( size0);
%         h1_0(nanloc) = h1;
%     end
%     
%     if ~isempty( h2 )
%         h2_0 = nan( size0);
%         h2_0(nanloc) = h2;
%     end
%     
% end
%%
% =========================================================================
% # TEST
% =========================================================================
% figure
% 
% xlim([-1 1]*1.5)
% ylim([-5 5])
% 
%       x = 0;
%       y = 0;
%       u = 1;
%       v = 2;
%       vel_scale = 1;
% 
%       %h = FUN_quiver_by_plotV2( x, y, u, v, vel_scale )
%       h = FUN_quiver_by_plotV2( x, y, u, v, 0, 'is_plot_head', 'true','head_length', 0, 'head_angle', 15, 'is_correct_angle', false )

%%
% =========================================================================
% # DEMO
% =========================================================================
% % 
% % % Here is a set of examples for FUN_quiver_by_plotV2_cmap_patch
% % 
% % % prepare data
% % 
% % theta = 0:10:360;
% % r = linspace(1,2,37);
% % 
% % x = r.*sind(theta);
% % y = r.*cosd(theta);
% % 
% % u = r.*cosd(theta);
% % v =-r.*sind(theta);
% % 
% % %% regular plot
% % figure
% % xlim([-3 3])
% % ylim([-3 3])
% % 
% % plot_scale = 1;
% % FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, plot_scale )
% % colorbar
% % 
% % %% Do not plot body of arrows
% % figure
% % xlim([-3 3])
% % ylim([-3 3])
% % 
% % plot_scale = 1;
% % FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, plot_scale, 'is_plot_body', false, 'is_plot_head', true, 'head_length', 0.5 )
% % colorbar
% % 
% % %% Do not plot head of arrows
% % figure
% % xlim([-3 3])
% % ylim([-3 3])
% % 
% % plot_scale = 1;
% % FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, plot_scale, 'is_plot_body', false, 'is_plot_head', true, 'head_length', 0.5, 'head_angle', 10 )
% % colorbar
% % 
% % 
% % %% Custom head angle
% % figure
% % xlim([-3 3])
% % ylim([-3 3])
% % 
% % plot_scale = 1;
% % FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, plot_scale, 'is_plot_body', true, 'is_plot_head', false, 'head_length', 0.5 )
% % colorbar
% % 
% % %% fill arrow head
% % figure
% % xlim([-3 3])
% % ylim([-3 3])
% % 
% % plot_scale = 1;
% % FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, plot_scale, 'is_plot_body', true, 'is_plot_head', true, 'head_length', 0.5, 'fill_head', true )
% % colorbar
% % 
% % %% custom values used to set colors
% % figure
% % xlim([-3 3])
% % ylim([-3 3])
% % 
% % plot_scale = 1;
% % FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, plot_scale, 'is_plot_body', true, 'is_plot_head', true, 'head_length', 0.5, 'fill_head', true, 'zval', rand(size(u)) )
% % colorbar
% % 
% % 
% % %% this function is based on patch. Any parameters works with patch can also be used here.
% % figure
% % xlim([-3 3])
% % ylim([-3 3])
% % 
% % plot_scale = 1;
% % FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, plot_scale, 'is_plot_body', true, 'is_plot_head', false, 'head_length', 0.5,  'linewidth', 2 )
% % colorbar



