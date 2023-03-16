


figure

% ----------------------------------------------------------------------

subplot(2,2,1)
xlim([0 2])
ylim([0 2])
FUN_quiver_by_plotV2( .5 , .5, .6, .8, 1, ...
                      'head_length_pixel', 50, 'fill_head', false, ...
                      'head_style', 1);
title('head_style=1, fill_head=false','Interpreter','none')
% ----------------------------------------------------------------------
             
subplot(2,2,2)
xlim([0 2])
ylim([0 2])
FUN_quiver_by_plotV2( .5 , .5, .6, .8, 1, ...
                      'head_length_pixel', 50, 'fill_head', true, 'is_plot_head',true, ...
                      'head_style', 1);
title('head_style=1, fill_head=true','Interpreter','none')
% ----------------------------------------------------------------------

subplot(2,2,3)
xlim([0 2])
ylim([0 2])     
FUN_quiver_by_plotV2( .5 , .5, .6, .8, 1, ...
                      'head_length_pixel', 50, 'fill_head', false, ...
                      'head_style', 2);
title('head_style=2, fill_head=false','Interpreter','none')

% ----------------------------------------------------------------------
subplot(2,2,4)
xlim([0 2])
ylim([0 2])   
FUN_quiver_by_plotV2( .5 , .5, .6, .8, 1, ...
                      'head_length_pixel', 50, 'fill_head', true, ...
                      'head_style', 2);
title('head_style=2, fill_head=true','Interpreter','none')

% ----------------------------------------------------------------------
FUN_easy_export_fig('Example_3.4.jpg','-m2');

                                           