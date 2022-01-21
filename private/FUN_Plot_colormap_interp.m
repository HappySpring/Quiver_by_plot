function res = FUN_Plot_colormap_interp( caxis_limit, cmap, val )

%% # check
if isvector( val )
else
   error('input must be vector'); 
end


%% # get x for colormap
x = linspace( caxis_limit(1), caxis_limit(2), size(cmap,1) );

%% # Linear interpolation
res = nan( length(val), 3 );

res(:,1) = interp1( x, cmap(:,1), val );
res(:,2) = interp1( x, cmap(:,2), val );
res(:,3) = interp1( x, cmap(:,3), val );

%% # Handle extrapolation 

if x(end) > x(1)
    
    
    loc_ex_hi = val >= x( end ) ;
    res(loc_ex_hi,1) =  cmap(end,1);
    res(loc_ex_hi,2) =  cmap(end,2);
    res(loc_ex_hi,3) =  cmap(end,3);
    
    loc_ex_lo = val <= x( 1 ) ;
    res(loc_ex_lo,1) =  cmap(1,1);
    res(loc_ex_lo,2) =  cmap(1,2);
    res(loc_ex_lo,3) =  cmap(1,3);
    
else
    
    loc_ex_hi = val >= x( 1 ) ;
    res(loc_ex_hi,1) =  cmap(1,1);
    res(loc_ex_hi,2) =  cmap(1,2);
    res(loc_ex_hi,3) =  cmap(1,3);
    
    loc_ex_lo = val <= x( end ) ;
    res(loc_ex_lo,1) =  cmap(end,1);
    res(loc_ex_lo,2) =  cmap(end,2);
    res(loc_ex_lo,3) =  cmap(end,3);
    
end


