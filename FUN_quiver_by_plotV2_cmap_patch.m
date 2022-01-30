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
% V3.00 by L. Chi
%         + This function has been integrated with FUN_quiver_by_plotV2. 
%
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

    is_rm_loadedd_param = true;
%                                                                    |input  | paramter       | default value | 
    [zval, varargin] = FUN_codetools_read_from_varargin(             varargin, 'zval',         sqrt(u.^2+v.^2), is_rm_loadedd_param ); % zval defines color 

% =========================================================================
% # call FUN_quiver_by_plotV2
% =========================================================================

    [h1,h2, uu, vv, hu, hv] = FUN_quiver_by_plotV2( x, y, u, v, vel_scale, 'zval', zval, varargin );


