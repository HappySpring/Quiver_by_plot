function [h1,h2, uu, vv, hu, hv] = FUN_quiver_by_plotV2( x, y, u, v, vel_scale, varargin )
% [h1, h2] = FUN_quiver_by_plotV2( x, y, u, v, vel_scale )
% [h1, h2] = FUN_quiver_by_plotV2( x, y, u, v, vel_scale, 'parameter')
% [h1, h2, uu, vv, hu, hv] = FUN_quiver_by_plotV2(  ... )
% 
% Plot quivers with "plot" functions. This make it easier to control the final results.
% 
% **Attention**************************************************************
%    Please set xlim and ylim before using this function and make sure that the 
%    values used for xlim and ylim will not change after this function is called.
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
%       is_plot_head: [default: true ]
%                  Plot arrow head, or not
%       is_plot_body: [default: true ]
%                  Plot arrow body, which is the part except arrow head.
%       fill_head: [default: false]
%                  If this is true, a filled triangle will be used as the head of each arrow.
%                  In this case, patch, instead of plot, is used to plot arrow heads.
%
%       head_length_pixel: [default: 0]
%                   It shares the same unit with u and v. (Thus, the vel_scale 
%                            will also be applied to head_length)
%                    0  : it is set according to the gca (width & length in pixel).
%                   >0  : it is the absolute length of the arrow head in unix of pixel
%                   -1~0: it is the percent of total vector length at each point.
%
%       head_length_min_pixel: [default: 0]
%                   minimal head length in unix of pixel
%
%       head_angle: 0.5* angles between the two wings of the arrow head (in degree)  [default: 10]
%
%       is_correct_angle [default: false]
%   
%       .... Remainings in varargin will be used as paramters for `plot` in general cases and `patch` for colorful arrows.
%
%
%   Parameters not recommended anymore.
%
%       head_length: [default: 0] (This is same as head_length_pixel but follows the same unit of the x axis.)
%                   It shares the same unit with u and v. (Thus, the vel_scale 
%                            will also be applied to head_length)
%                    0  : it is set according to the gca (width & length in pixel).
%                   >0  : it is the absolute length of the arrow head
%                   -1~0: it is the percent of total vector length at each point.
%
%       head_length_min: [default: 0] (This is same as head_length_min_pixel but follows the same unit of the x axis.)
%                   minimal head length following the same unit of the x axis.
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
%       h = FUN_quiver_by_plotV2( x, y, u, v, vel_scale, 'is_plot_body', 'true','head_length', 0, 'head_angle', 15, 'is_correct_angle', false )
% =========================================================================


% =========================================================================
% V3.10 by L. Chi: add a new arrow style. 
%                  Fix a minor bug: a warning message may appear when there is nothing to worry about.
% V3.01 by L. Chi: update comments. Some minor edits.
% V3.00 by L. Chi: 
%         - fix a bug in calculating head length defined by the parameter "head_length"
%           when the ratio of pixels per distance between x and y distance
%           are very large, the head length be longer than what expected.
%
%         - [New] Support colorful arrows following the colormap.
%
%         - [important] "is_correct_angle_by_edit_v" is set to false by default.
%         - add a warning if 'XLimMode' or "YLimMode' is set to 'auto'.
%         - move "hold on" after the first time "plot" is called to make
%         the behavior matches exception.
%
% V2.67 by L. Chi: Add is_arrow_origin_at_xy
% V2.66 by L. Chi: update head length parameters
%         - Update default behaiver. The default head length is
%                    40% of the vector length at each point.
%         - If the input 'head_length' is at [-1, 0), then it is 
%                    the percent of total vector length at each point.
%         - Support minimal head length (parameter: 'head_length_min')
% V2.65b by L. Chi: correct an error in estimating `head_length`
% V2.65 by L. Chi (06/18/2021)
%         - update methods for estimating xlimit and ylimit if an axes
%         has not be created.
%         - clean some useless codes
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
% # Optional Paramters
% =========================================================================
    
    % Parameter                                                                     Parameter,    Default value 
    [zval,         varargin]          = FUN_codetools_read_from_varargin( varargin, 'zval',             [] ); % for colorful arrows 
    [is_plot_head, varargin]          = FUN_codetools_read_from_varargin( varargin, 'is_plot_head',    true  );
    [is_plot_body, varargin]          = FUN_codetools_read_from_varargin( varargin, 'is_plot_body',    true  );
    [fill_head,    varargin]          = FUN_codetools_read_from_varargin( varargin, 'fill_head',       false );

    [line_color, varargin]            = FUN_codetools_read_from_varargin( varargin, 'color',           'k' );
    [head_length, varargin]           = FUN_codetools_read_from_varargin( varargin, 'head_length',      0 );
    [head_length_min, varargin]       = FUN_codetools_read_from_varargin( varargin, 'head_length_min',  0 );
    [head_length_pixel,varargin]      = FUN_codetools_read_from_varargin( varargin, 'head_length_pixel', 0 );
    [head_length_pixel_min, varargin] = FUN_codetools_read_from_varargin( varargin, 'head_length_pixel_min', 0 );
    % style 
    % style - 1 (default): a simple arrow
    % style - 2: see readme.md section 3.5 Arrow style
    [head_style, varargin]            = FUN_codetools_read_from_varargin( varargin, 'head_style',      1  );
    [head_length_2, varargin]         = FUN_codetools_read_from_varargin( varargin, 'head_length_2',   0.5);
    
    [head_angle, varargin]            = FUN_codetools_read_from_varargin( varargin, 'head_angle',      15 );
    [is_correct_angle, varargin]      = FUN_codetools_read_from_varargin( varargin, 'is_correct_angle', false );
    [is_correct_angle_by_edit_v, varargin] = FUN_codetools_read_from_varargin( varargin, 'is_correct_angle_by_edit_v', false );
    
    % velocity lower than this value will be ignored.
    [min_zval_displayed, varargin]       = FUN_codetools_read_from_varargin( varargin, 'min_zval_displayed', nan );
    
    [interval,   varargin]            = FUN_codetools_read_from_varargin( varargin, 'interval',          1 );
    [interval_x, varargin]            = FUN_codetools_read_from_varargin( varargin, 'interval_x',        interval );
    [interval_y, varargin]            = FUN_codetools_read_from_varargin( varargin, 'interval_y',        interval );
    
    % gca_pos_now = plotboxpos(gca);
    % This can be given as an input paramter instead of finding it every
    % time in this function. This helps speedup the function.
    [gca_pos_now, varargin]           = FUN_codetools_read_from_varargin( varargin, 'gca_pos_now',        [] );
    
    
    % if it is true, arrows origin from input (x,y). Otherwise, the center of arrow bodies locate at (x,y) 
    [is_arrow_origin_at_xy, varargin] = FUN_codetools_read_from_varargin( varargin, 'is_arrow_origin_at_xy', true );

    
    %[is_plot_arrow_body, varargin] = FUN_codetools_read_from_varargin( varargin, 'is_plot_arrow_body', true );
    [is_plot_via_arrayfun, varargin] = FUN_codetools_read_from_varargin( varargin, 'is_plot_via_arrayfun', false );

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
    
    if ~isempty( zval ) 
       is_colorful_arrows = true;
    else
       is_colorful_arrows = false;
    end
    
%%
% =========================================================================
% # initialization
% =========================================================================

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
        if is_colorful_arrows
            zval = zval(1:interval_x:end,1:interval_y:end);
        end
    end

% ---- ## Make velocity vectors -------------------------------------------
    %size0 = size(x);
    x = x(:).'; 
    y = y(:).';
    u = u(:).'; 
    v = v(:).';
    if is_colorful_arrows
        zval = zval(:).';
    end
    
% ---- ## Check -----------------------------------------------------------
    plot_size = size( u );

    if     isequal( size( x ), plot_size ) ...
        && isequal( size( y ), plot_size ) ...
        && isequal( size( u ), plot_size ) ...
        && isequal( size( v ), plot_size )
        % Pass: x, y, u, v must share same size
    else
        error('The input variable x/y/u/v should share same size')
    end
    
    if is_colorful_arrows && ~isequal( size( zval ), plot_size )
        error('The input variable x/y/u/v should share same size')
    end
    
    
% ---- ## Remove NaNs -----------------------------------------------------
    nanloc = isnan( x ) | isnan( y ) | isnan( u ) | isnan( v );
    
    if ~isnan( min_zval_displayed )
        if isempty( zval ) 
            zval = sqrt( u.^2 + v.^2 );
        end
        nanloc = nanloc | zval < min_zval_displayed;
    end
    
    x( nanloc ) = [];
    y( nanloc ) = [];
    u( nanloc ) = [];
    v( nanloc ) = [];
    
    if is_colorful_arrows
        zval(nanloc)= [];
    end

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
            tem = get(gcf,'children');
            
            if isempty(tem)
                xlimit = [0, 1];
                ylimit = [0, 1];
            else
                xlimit = get(gca,'xlim');
                ylimit = get(gca,'ylim');
            end

            if xlimit(1) == xlimit(2)
                xlimit = [xlimit(1) xlimit(1)];
            end
            
            if ylimit(1) == ylimit(2)
                ylimit = [ylimit(1) ylimit(1)];
            end
            
            xlim_diff = xlimit(2) - xlimit(1);
            ylim_diff = ylimit(2) - ylimit(1);
   
       
            delx = xlim_diff/3;
            dely = ylim_diff/3;
        else
            delx = diff([min(x(:)) max(x(:))])/n;
            dely = diff([min(y(:)) max(y(:))])/m;
        end
        del = delx.^2 + dely.^2;
        
        % the variables cleared in the next line are only used for auto size and my be modifed
        % later.
        clear delx dely xlim_diff ylim_diff xlimit ylimit
         
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
    
% ---- ## get figure and gca info -----------------------------------------
    % ### create axis if it does not exist
    tem = get(gcf,'children');

    if isempty(tem)
       % set xlimit/ylimit automatically
       tx = [ x; x+u*1.2 ];
       ty = [ y; y+v*1.2 ];
       xlimit = [ nanmin(tx(:)) nanmax(tx(:)) ];
       ylimit = [ nanmin(ty(:)) nanmax(ty(:)) ];
       
       xlim( xlimit );
       ylim( ylimit );
       
    else
       % read xlimit/ylimit from the current figure
       xlimit = get(gca,'xlim');
       ylimit = get(gca,'ylim');
    end
    
    if strcmpi(get(gca,'XLimMode'), 'auto')
       warning('XLimMode is auto!It is high recommended to set a xlim manually before using this function!');
    end
    
    if strcmpi(get(gca,'YLimMode'), 'auto')
       warning('YLimMode is auto! It is high recommended to set a ylim manually before using this function!');
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
   if isempty( gca_pos_now )
        gca_pos_now = plotboxpos( gca );
   else
      % this is given as an input parameter to speedup the function.
   end
   
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
    
% ---- ## Ploit quiver body -----------------------------------------------
    if is_plot_body
        if is_correct_angle_by_edit_v
            up = u;
            vp = v.*r_correct_v; % r_correct_v is 1 if is_correct_angle is set to flase.
        else
            up = u./r_correct_v; % r_correct_v is 1 if is_correct_angle is set to flase.
            vp = v;
        end
    else 
        if  ( isempty( head_length ) || head_length == 0 ) && isempty(head_length_pixel)
            % head_length = 1;
            warning('Please set head_length manually when is_plot_arrow_body = false')
        end
        up = zeros(size(u));
        vp = zeros(size(v));
    end

    if is_arrow_origin_at_xy == true
        % origins of arrows locate at (x,y)
        x_offset = 0;
        y_offset = 0;
    else
        % centers of arrows locate at (x,y)
        x_offset = -up./2;
        y_offset = -vp./2;    
    end
    
    if is_plot_body
        
        if is_arrow_origin_at_xy == true
            uu = [x ; x+up; nan(size(u))];
            vv = [y ; y+vp; nan(size(v))];
        else
            uu = [x+x_offset ; x+up+x_offset; nan(size(u))];
            vv = [y+y_offset ; y+vp+y_offset; nan(size(v))];
        end
        
        if is_colorful_arrows            
            % plot colorful arrows
            h1 = patch( 'XData', uu(1:2,:), 'YData', vv(1:2,:), 'CData', repmat(zval(:).',2,1), 'FaceColor','none', 'edgecolor','flat', varargin{:} );
        else
            % plot reguar arrows
            if is_plot_via_arrayfun
                % Using arrayfun will return h1 as a matrix. This make it possbile to chagne the color
                % of each arrows.
                % h1 is a 1 x N matrix
                h1 = arrayfun( @(x01,x02,y01,y02)plot( [x01 x02], [y01 y02],'-', 'parent', cax, 'color', line_color, varargin{:}), uu(1,:), uu(2,:), vv(1,:), vv(2,:)  );
            else
                % h1 is in size of 1 x 1
                h1 = plot( uu(:), vv(:), '-', 'parent', cax, 'color', line_color, varargin{:} );
            end
        end
    else
        uu = [];
        vv = [];
        h1 = [];
    end
    
    hold on
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
        
        % set default value
        if isempty( head_length_pixel ) || head_length_pixel ==0
            
            % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            % head_length only works if head_length_pixel is not defined!
            %
            % This block is kept here for compatibility purpose only. 
            % It won't be updated anymore 
            
            if  isempty( head_length ) || head_length == 0 
                head_length = -0.4 ; % 40% of the arrow length
            end

            if head_length < 0 && head_length >= -1

                r_head_L = abs(head_length); 
                head_length_pixel = sqrt( (up.*pix_per_x).^2 + (vp.*pix_per_y).^2 ) * r_head_L;       

            elseif head_length < -1
                error('Unexpected head_length!');

            else
                % velocity scale has been applied to u and v (also up and vp).
                % Thus, it is should not applied used here.
                % % head_length = head_length * vel_scale;
                head_length = head_length * vel_scale;
                head_length_pixel = head_length .* pix_per_x; 
                                
            end
            
            % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            
        elseif head_length_pixel > 0
            % head_length_pixel = head_length_pixel ;
            
        elseif head_length_pixel < 0 && head_length_pixel >= -1
            r_head_L = abs(head_length_pixel);
            head_length_pixel = sqrt( (up.*pix_per_x).^2 + (vp.*pix_per_y).^2 ) * r_head_L;
        else
            error('Unexpected head_length_pixel!');
        end
        
        % set min head length 
        if head_length_min > 0
            % Please do not set head_length_min. 
            % It is left here for compatiblility purpose. 
            head_length_min = head_length_min * vel_scale;
            head_length_pixel_min = head_length_min * pix_per_x;
        end
        
            
        if head_length_pixel_min > 0
            head_length_pixel( head_length_pixel < head_length_pixel_min ) = head_length_pixel_min;    
        end

        % ### Generate arrays for arrow heads -----------------------------
        if head_style == 1 % default, simple arrow like ->
            
            npa = 3; % number of point for each arrow.
            
            if is_plot_body
                % head is placed at the end of the arrow
                hu = [x+up-head_length_pixel./pix_per_x.*cosd( tem_ang1+head_angle );  x+up; ...
                      x+up-head_length_pixel./pix_per_x.*cosd( tem_ang1-head_angle );  NaN(size(x))];

                hv = [y+vp-head_length_pixel./pix_per_y.*sind( tem_ang1+head_angle );  y+vp; ...
                      y+vp-head_length_pixel./pix_per_y.*sind( tem_ang1-head_angle );  NaN(size(x))];
            else
                % the middle point of the head triangle is placed at (x,y)
                %
                % up, vp are kept in calculating tem_hu, tem_hv to keep a correct size.
                tem_hu = [up-head_length_pixel./pix_per_x.*cosd( tem_ang1+head_angle );  up; ...
                          up-head_length_pixel./pix_per_x.*cosd( tem_ang1-head_angle );  NaN(size(x))];

                tem_hv = [vp-head_length_pixel./pix_per_y.*sind( tem_ang1+head_angle );  vp; ...
                          vp-head_length_pixel./pix_per_y.*sind( tem_ang1-head_angle );  NaN(size(x))];

                hu = x + tem_hu - nanmean(tem_hu(1:3,:));
                hv = y + tem_hv - nanmean(tem_hv(1:3,:));
            end

        elseif head_style == 2
            npa = 5; % number of point for each arrow.
            
            if is_plot_body
                % head is placed at the end of the arrow
                hu = [x+up-head_length_pixel./pix_per_x.*cosd( tem_ang1+head_angle );  
                      x+up; ...
                      x+up-head_length_pixel./pix_per_x.*cosd( tem_ang1-head_angle );   ...
                      x+up-head_length_pixel./pix_per_x.*cosd( tem_ang1 )*head_length_2;...
                      x+up-head_length_pixel./pix_per_x.*cosd( tem_ang1+head_angle );  ...
                      NaN(size(x))];

                hv = [y+vp-head_length_pixel./pix_per_y.*sind( tem_ang1+head_angle );        % p1
                      y+vp; ...                                                              % p2
                      y+vp-head_length_pixel./pix_per_y.*sind( tem_ang1-head_angle ); ...    % p3
                      y+vp-head_length_pixel./pix_per_y.*sind( tem_ang1 )*head_length_2; ... % p4
                      y+vp-head_length_pixel./pix_per_y.*sind( tem_ang1+head_angle );  ...   % p1
                      NaN(size(x))];
            else
                % the middle point of the head triangle is placed at (x,y)
                %
                % up, vp are kept in calculating tem_hu, tem_hv to keep a correct size.
                tem_hu = [up-head_length_pixel./pix_per_x.*cosd( tem_ang1+head_angle );  
                          up; ...
                          up-head_length_pixel./pix_per_x.*cosd( tem_ang1-head_angle ); ... 
                          up-head_length_pixel./pix_per_x.*cosd( tem_ang1 )*head_length_2; ... 
                          up-head_length_pixel./pix_per_x.*cosd( tem_ang1+head_angle );  ...
                          NaN(size(x))];

                tem_hv = [vp-head_length_pixel./pix_per_y.*sind( tem_ang1+head_angle );  
                          vp; ...
                          vp-head_length_pixel./pix_per_y.*sind( tem_ang1-head_angle );  ...
                          vp-head_length_pixel./pix_per_y.*sind( tem_ang1 )*head_length_2;  ...
                          vp-head_length_pixel./pix_per_y.*sind( tem_ang1+head_angle );...
                          NaN(size(x))];

                hu = x + tem_hu - nanmean(tem_hu(1:3,:)); % tem_hu(4,:) is added for tuning arrow style. It should not be counted here.
                hv = y + tem_hv - nanmean(tem_hv(1:3,:));
            end
        else
            error('Unknown head_style');            
        end
        
        if is_arrow_origin_at_xy == false
           hu = bsxfun(@plus, hu, x_offset);
           hv = bsxfun(@plus, hv, y_offset);
           %hu(1,:) = hu(1,:) + x_offset;
           %hu(2,:) = hu(2,:) + x_offset;
           %hu(3,:) = hu(3,:) + x_offset;
           % 
           %hv(1,:) = hv(1,:) + y_offset;
           %hv(2,:) = hv(2,:) + y_offset;
           %hv(3,:) = hv(3,:) + y_offset;
        end

        
        % ### Plot head of arrows
        hold(cax,'on')
        
        % ### plot
        % Using arrayfun will return h1 as a matrix. This make it possbile to chagne the color
        % of each arrows.        
        
        if is_colorful_arrows
            
            if fill_head
                h2  = patch( 'XData', hu(1:npa,:), 'YData', hv(1:npa,:), 'CData', repmat(zval(:).',npa,1), 'FaceColor','flat', 'edgecolor','none', varargin{:} );
            else
                h21 = patch( 'XData', hu(1:2,:), 'YData', hv(1:2,:), 'CData', repmat(zval(:).',2,1), 'FaceColor','none', 'edgecolor','flat', varargin{:} );
                h22 = patch( 'XData', hu(2:3,:), 'YData', hv(2:3,:), 'CData', repmat(zval(:).',2,1), 'FaceColor','none', 'edgecolor','flat', varargin{:} );
                if head_style == 1                   
                   h2 = [h21; h22];
                elseif head_style == 2
                   h23 = patch( 'XData', hu(3:4,:), 'YData', hv(3:4,:), 'CData', repmat(zval(:).',2,1), 'FaceColor','none', 'edgecolor','flat', varargin{:} );
                   h24 = patch( 'XData', hu([1 4],:), 'YData', hv([1 4],:), 'CData', repmat(zval(:).',2,1), 'FaceColor','none', 'edgecolor','flat', varargin{:} );
                   h2 = [h21; h22; h23; h24];
                else
                   error('Unknown head_style!');
                end
            end
            
        else
        
            if is_plot_via_arrayfun
            % h1 is a 1 x N matrix
                if head_style == 1
                    if fill_head
                        h2 = arrayfun( @(x01,x02,x03,y01,y02,y03)patch([x01;x02;x03],[y01;y02;y03], line_color, 'parent', cax, 'edgecolor','none', varargin{:}), hu(1,:), hu(2,:), hu(3,:), hv(1,:), hv(2,:), hv(3,:), 'UniformOutput',false);
                        h2 = cat(1, h2{:});
                    else    
                        %h2 = arrayfun( @(x01,x02,x03,y01,y02,y03)plot([x01,x02,x03],[y01,y02,y03],'-','parent',cax, 'color', line_color, varargin{:}), hu(1,:), hu(2,:), hu(3,:), hv(1,:), hv(2,:), hv(3,:) );
                        h2 = arrayfun( @(x01,x02,x03,y01,y02,y03)plot([x01,x02,x03],[y01,y02,y03],'-','parent',cax, 'color', line_color, varargin{:}), hu(1,:), hu(2,:), hu(3,:), hv(1,:), hv(2,:), hv(3,:), 'UniformOutput',false);
                        h2 = cat(1, h2{:});
                    end
                elseif head_style == 2
                    if fill_head
                        h2 = arrayfun( @(x01,x02,x03,x04,x05,y01,y02,y03,y04,y05)patch([x01;x02;x03;x04;x05],[y01;y02;y03;y04;y05], line_color, 'parent', cax, 'edgecolor','none', varargin{:}), hu(1,:), hu(2,:), hu(3,:), hu(4,:), hu(5,:), hv(1,:), hv(2,:), hv(3,:), hv(4,:), hv(5,:), 'UniformOutput',false);
                        h2 = cat(1, h2{:});
                    else    
                        %h2 = arrayfun( @(x01,x02,x03,y01,y02,y03)plot([x01,x02,x03],[y01,y02,y03],'-','parent',cax, 'color', line_color, varargin{:}), hu(1,:), hu(2,:), hu(3,:), hv(1,:), hv(2,:), hv(3,:) );
                        h2 = arrayfun( @(x01,x02,x03,x04,x05,y01,y02,y03,y04,y05)plot([x01,x02,x03,x04,x05],[y01,y02,y03,y04,y05],'-','parent',cax, 'color', line_color, varargin{:}), hu(1,:), hu(2,:), hu(3,:), hu(4,:), hu(5,:), hv(1,:), hv(2,:), hv(3,:),hv(4,:), hv(5,:),'UniformOutput',false);
                        h2 = cat(1, h2{:});
                    end
                else
                    error('unknown head_style!');
                end
            else
            % h1 is in size of 1 x 1
                if fill_head
                    
                    h2 = patch(hu(1:npa,:),hv(1:npa,:), line_color, 'parent', cax, 'edgecolor','none', varargin{:});
                else
                    h2 = plot(hu(:),hv(:),'-','parent',cax, 'color', line_color, varargin{:});
                end
            end
        
        end
        
    else
        hu = [];
        hv = [];
        h2 = [];
    end
    
%%
% =========================================================================
% # Finalization
% =========================================================================

    if ~hold_state
        hold(cax,'off');
        view(cax,2); 
        set(cax,'NextPlot',next); 
    end

    %h = [h1;h2];
    h1 = h1(:);
    h2 = h2(:);

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
% # TESTING
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

