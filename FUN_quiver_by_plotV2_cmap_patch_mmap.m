function h=FUN_quiver_by_plotV2_cmap_patch_mmap(long,lat,u,v,vel_scale,varargin)
% This is an edited version of m_quiver.
% It is modified to be comptible with FUN_quiver_by_plotV2
% V1.00 07/09/2021: first version
%
% =========================================================================
% Notes from m_quiver
% =========================================================================
% M_QUIVER Makes a quiverplot on a map (QUIVER-style)
%    M_QUIVER(LONG,LAT,U,V) plots velocity vectors as arrows with components 
%    (U,V) at the points (LONG,LAT) on the currently defined map.  The 
%    matrices LONG,LAT,U,V must all be the same size. U and V contain the 
%    eastward and northward components of velocity (in m/s or equivalent, NOT
%    degrees lat/long  per sec or equivalent). Arrow scaling is automatic.
% 
%    M_QUIVER(X,Y,U,V,S) automatically scales the arrows to fit within the 
%    grid and then stretches them by S.  Use S=0 to plot the arrows without 
%    the automatic scaling; In this case the scaling is 1 unit/degree 
%    latitude. Note that we do not scale arrows with respect to map 
%    coordinates! Instead, the arrows will correspond better to actual motions
%    over some time step. The tradeoff is that a single scale arrow cannot
%    be accurate for the entire map (M_VEC scales arrows according to
%    map coordinates).
% 
%    M_QUIVER(...,LINESPEC) uses the plot linestyle specified for
%    the velocity vectors.  Any marker in LINESPEC is drawn at the base
%    instead of an arrow on the tip.  Use a marker of '.' to specify
%    no marker at all.  See PLOT for other possibilities. M_QUIVER is a wrapper
%    for QUIVER - for fancier arrows it is possible to replace the call to 
%    QUIVER with one to another routine that draws fancy arrows, e.g. 
%    ARROW (from TMW user-contrib software archive), or to use M_VEC.
%
%    M_QUIVER(...,'filled') fills any markers specified.
% 
%    H = M_QUIVER(...) returns a vector of line handles.
% 
%    See also QUIVER, M_VEC

% Rich Pawlowicz (rich@ocgy.ubc.ca) 20/Jan/97
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%

% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)
% 7/jul/06 - changed angle calc to work correctly very near boundaries.
% 12/jul/06 - fixed a factor of 10 error that crept into the length of unscaled
%             arrows between version 1.3f and 1.4a (pointed out D. Kaplan).

global MAP_PROJECTION MAP_VAR_LIST

% Have to have initialized a map first

if isempty(MAP_PROJECTION)
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end

%%
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
% ** Please note that zval must not be passed to FUN_quiver_by_plotV2_cmap_patch directly. **
% It must be handled before both (u,v) is converted to mmap coordinate.
%
is_rm_loadedd_param = true; % this must be true for zval.
[zval, varargin]       = FUN_codetools_read_from_varargin( varargin, 'zval', sqrt(u.^2+v.^2), is_rm_loadedd_param );


if isvector( long ) && length(long) == size(u,1) && isvector(lat) && length(lat) == size(u,2)
    [long,lat]=meshgrid(long,lat);
    long = long';
    lat = lat';
end
    

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%% 
[X,Y]=m_ll2xy(long,lat,'clip','point');

% This is the old way, now replaced  - RP 7/jun/06
%[XN,YN]=m_ll2xy(long,lat+.001,'clip','point');
%[XE,YE]=m_ll2xy(long+(.001)./cos(lat*pi/180),lat,'clip','point');

%mU=u.*(XE-X)*100 + v.*(XN-X)*100;
%mV=u.*(YE-Y)*100 + v.*(YN-Y)*100;

[XN ,YN ]=m_ll2xy([long(:) long(:)]',[lat(:) lat(:)+.001]','clip','off');
[XE ,YE ]=m_ll2xy([long(:) long(:)+(.001)./cos(lat(:)*pi/180)]',[lat(:) lat(:)]','clip','off');
mU=u.*reshape(diff(XE),size(lat))*1000 + v.*reshape(diff(XN),size(lat))*1000;
mV=u.*reshape(diff(YE),size(lat))*1000 + v.*reshape(diff(YN),size(lat))*1000;

% Only pass over mapped stuff (otherwise auto-scaling doesn't work)

% ii=isfinite(X(:)) & isfinite(u(:));
 
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Replace the following command
%      h=quiver(X(ii),Y(ii),mU(ii),mV(ii),varargin{:});
% by this one
h=FUN_quiver_by_plotV2_cmap_patch(X,Y,mU,mV,vel_scale, 'zval', zval, varargin{:});
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

set(h,'tag','m_quiver');

if nargout==0
 clear h
end
