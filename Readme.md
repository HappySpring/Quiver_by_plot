[toc]

## 1. Introduction

### 1.1 Features/Why I wrote this

+ More flexibility than the built-in function `quiver`.

  The arrows are generated by `plot`/`patch`, which provides the possibility that everything can be adjusted. Here are some of the examples:

  + The length of head (fixed length, or a ratio to the arrow length, etc.)
  + Angles between two edges of head. 
  + Show/hide the head of arrows.

  + Show/hide the body of arrows. 
  + Other properties in `plot` or `patch`.
  + ...

+ Colorful arrows adopted to current colormap

+ Support creating vector time series

+ "Perfect" arrows

  

### 1.2 Functions included

+ `FUN_quiver_by_plotV2.m`

  ​	Plot single color arrows by `plot`

+ `FUN_quiver_by_plotV2_mmap.m`

  ​	This is same as `FUN_quiver_by_plotV2.m` but compatible with [`mmap`](https://www.eoas.ubc.ca/~rich/mapug.html)

+ `FUN_quiver_by_plotV2_cmap_patch.m`

  ​	Plot colorful arrows based on `patch`. The color of arrows depends on colormap.

+ `FUN_quiver_by_plotV2_cmap_patch_mmap.m`

  ​	This is same as `FUN_quiver_by_plotV2_cmap_patch.m` but compatible with [`mmap`](https://www.eoas.ubc.ca/~rich/mapug.html)

+ `FUN_quiver_by_plotV2_preset_cmap.m`

  ​	[not recommended!] This is similar to `FUN_quiver_by_plotV2_cmap_patch` but the color is defined based on a fixed colormap, which must be specified 



### Limitation

**To create a *perfect* arrow, the `xlim` and `ylim` must be defined explicitly and it should not be changed after calling this toolbox.**



## 2. Install

You need to add this to the searching path of your Matlab environment before using them. The `private` folder must be put in the same folder as `FUN_quiver_by_plotV2*.m`.  The subfolders (`private`, `doc`) should not be added to the searching path. It can be done by two ways:

### If you have GUI access to Matlab : 

Click "Home tab> Set Path". It will open a dialog for setting the path. Then, click "Add Folder...", add the root path of this package (the folder contains functions `FUN_quiver_by_plotV2*.m`), then click "Save" near the bottom of the dialog.  

### If you are in a command line environment:

+ Method 1 (recommended):

  ```
  addpath('/path/to/Quiver_by_plotV2/');
  savepath
  ```

+ Method 2: 
  Matlab will run `startup.m` during boot automatically if the file exists. Thus, you can add `addpath('/path/to/Quiver_by_plotV2/');` to the `startup.m` file and make sure that the `startup.m` is put in existing searching path. This provide more flexibility.



## 3. How to use this toolbox

### 3.1 Plot arrows without colormap (`FUN_quiver_by_plotV2.m`)

`FUN_quiver_by_plotV2`


```matlab
[h1,h2, uu, vv, hu, hv] = FUN_quiver_by_plotV2( x, y, u, v)
[h1,h2, uu, vv, hu, hv] = FUN_quiver_by_plotV2( x, y, u, v, vel_scale)
[h1,h2, uu, vv, hu, hv] = FUN_quiver_by_plotV2( x, y, u, v, vel_scale, ...)
```

#### INPUT:
​       x, y, u, v : input data
​       vel_scale [optional, default: 1 ]
​                  0: auto scale  
​                  1: No scale, the original u/v will be used
​                  otherwise: this scale will be applied to u and v.

##### Other parameters

| Parameter                  | Default value | Notes                                                        |
| -------------------------- | ------------- | ------------------------------------------------------------ |
| is_plot_head               | true          | show/hide head of arrows                                     |
| is_plot_body               | true          | show/hide body of arrows                                     |
| fill_head                  | false         | If this is true, a filled triangle will be used as the head of each arrow.<br/>In this case, patch, instead of plot, is used to plot arrow heads. |
| color                      | 'k'           | color of arrows                                              |
| head_length_pixel          | 0             | length of (the wing of) arrow head in units of pixel<br /><br />It shares the same unit with u and v. (Thus, the vel_scale will also be applied to head_length)<br/>        0  : it is set according to the gca (width & length in pixel).<br/>        >0  : it is the absolute length of the arrow head<br/>        -1~0: it is the percent of total vector length at each point. |
| head_length_min_pixel      | 0             | minimal length of arrow head length in units of pixel.       |
| head_angle                 | 15            | Angle between arrow head and arrow body<br />`2*head_angle` is the angle between the two wings of the arrow head. |
| is_correct_angle           | false         | This is for vector time series and some special cases.       |
| is_correct_angle_by_edit_v | false         |                                                              |
| interval                   | 1             | interval                                                     |
| interval_x                 | interval      | interval for the first dim (usually x)                       |
| interval_y                 | interval      | interval for the second dim (usually y)                      |
| is_plot_via_arrayfun       | false         | Call `plot` via  `arrayfun`. This is useful if you want to change properties of a subsample of arrows.<br />**false**: The returned handles h1 & h2 are 1 x 1 matrix (one for array body and one for arrow head)<br />**true**: The returned handles h1 & h2 are  [m x 1] matrix (one for array body and one for arrow head)<br />Note: this does not support colorful plots yet. |
| .... [varargin]            | ....          | Other parameters (e.g., 'linewidth') applicable to `plot` for regular plots or `patch` for colorful arrows. |

##### Parameters not recommended

| Parameter       | Default value | Notes                                                        |
| --------------- | ------------- | ------------------------------------------------------------ |
| head_length     | 0             | length of (the wing of) arrow head in the same unit as x-axis<br /><br />It shares the same unit with u and v. (Thus, the vel_scale will also be applied to head_length)<br/>        0  : it is set according to the gca (width & length in pixel).<br/>        >0  : it is the absolute length of the arrow head<br/>        -1~0: it is the percent of total vector length at each point. |
| head_length_min | 0             | minimal length of arrow head<br />                           |



#### Output

+ h1: handle array for arrow body
+ h2: handle array for arrow head
+ uu [ 3 x N ]: x coordinate for arrow body, the 3 values in each col are [ x_origin; x_head; Nan]​
+ vv [ 3 x N ]: y coordinate for arrow body, the 3 values in each col are [ y_origin; y_head; Nan]
+ hu [ 4 x N ]: x coordinate for arrow head, the 4 values in each col are [ x_head_wing1, x_head_center, x_head_wing2; Nan]​ 
+ hv [ 4 x N ]: y coordinate for arrow head, the 4 values in each col are [ y_head_wing1, y_head_center, y_head_wing2; Nan]

​     you can re-plot the arrows by uu,vv hu, hv:

   ```matlab
      plot( uu(:), vv(:) ); % plot arrow body
      hold on
      plot( hu(:), hv(:) ); % plot arrow head
   ```



#### Demos

Data used in this demo is generated by the following codes

```matlab
x = -1:0.2:1;
y = -1:0.2:1;

[X,Y] = meshgrid(x,y);
X=X';
Y=Y';

theta = atan2d(Y,X);

r = sqrt( X.^2 + Y.^2);

u =  sind( theta ) .* r / 100;
v = -cosd( theta ) .* r / 100;

x = x*10 - 75;
y = y*10 + 35;
```

This is also available here: [Data_demo_1.mat](doc\Data_demo_1.mat).

##### Demo - plot arrows. 

```matlab
figure
xlim([-88 -62]);
ylim([22 48]);

FUN_quiver_by_plotV2( x, y, u, v, 0 )

title('Example 3.1: FUN_quiver_by_plotV2','Interpreter','none');

FUN_easy_export_fig('Example_3.1.jpg','-m2');
```

<img src="doc/Demo_3.1_arrows/Example_1.1.jpg" alt="Example_1.2a" style="zoom: 50%;" />

##### Demo - set absolute scale 

```matlab
vel_scale1 = 200;
vel_scale2 = 300;

figure
xlim([-88 -62]);
ylim([22 48]);

hold on

hd1 = FUN_quiver_by_plotV2( x, y, u, v, vel_scale1, 'color','m' );
hd2 = FUN_quiver_by_plotV2( x, y, u, v, vel_scale2, 'color','k' );

lid = legend( [hd1, hd2], 'scale=200', 'scale=300');
set(lid,'Interpreter','none');
box on

title('Example 1.2 "scale"');
```



<img src="doc\Demo_3.1_arrows\Example_1.2.jpg" alt="Example_1.2" style="zoom:50%;" />

##### Demo - set head length in pixel (parameter: 'head_length_pixel'). 

```matlab
head_length1 = 20;
head_length2 = 10;

figure
xlim([-88 -62]);
ylim([22 48]);

hold on

hd1 = FUN_quiver_by_plotV2( x, y, u, v, 0, 'head_length_pixel', head_length1, 'color','m' );
hd2 = FUN_quiver_by_plotV2( x, y, u, v, 0, 'head_length_pixel', head_length2, 'color','k' );

lid = legend( [hd1, hd2], 'head_length_pixel=20', 'head_length_pixel=10');
set(lid,'Interpreter','none');
box on

title('Example 1.3: "head_length (positive)"','Interpreter','none');		
```



<img src="E:\working_dir_sync\2022-01-30_quiver_by_plot\Quiver_by_plot\doc\Demo_3.1_arrows\Example_1.3a.jpg" alt="Example_1.3a" style="zoom:50%;" />

```matlab
head_length1 = -0.6;
head_length2 = -0.3;

figure
xlim([-88 -62]);
ylim([22 48]);

hold on

hd1 = FUN_quiver_by_plotV2( x, y, u, v, 0, 'head_length_pixel', head_length1, 'color','m' );
hd2 = FUN_quiver_by_plotV2( x, y, u, v, 0, 'head_length_pixel', head_length2, 'color','k' );

lid = legend( [hd1, hd2], 'head_length_pixel=-0.6', 'head_length_pixel=-0.3');
set(lid,'Interpreter','none');
box on

title('Example 1.3: "head_length (negative)"','Interpreter','none');
```

<img src="doc\Demo_3.1_arrows\Example_1.3b.jpg" alt="Example_1.3b" style="zoom:50%;" />

##### Demo - set head angle (parameter: 'head_angle'). 

```matlab
head_angle1 = 20;
head_angle2 = 40;

figure
xlim([-88 -62]);
ylim([22 48]);

hold on

hd1 = FUN_quiver_by_plotV2( x, y, u, v, 0, 'head_angle', head_angle1, 'color','m' );
hd2 = FUN_quiver_by_plotV2( x, y, u, v, 0, 'head_angle', head_angle2, 'color','k' );

lid = legend( [hd1, hd2], 'head_angle=20', 'head_angle=40');
set(lid,'Interpreter','none');
box on

title('Example 1.4: "head_angle"','Interpreter','none');
```

<img src="doc\Demo_3.1_arrows\Example_1.4.jpg" alt="Example_1.4" style="zoom:50%;" />

##### Demo - hide arrow head (parameter: 'is_plot_head'). 

```matlab
figure
xlim([-88 -62]);
ylim([22 48]);

hold on

hd1 = FUN_quiver_by_plotV2( x, y, u, v, 0, 'is_plot_head', false );
box on

title('Example 1.5: "is_plot_head=false"','Interpreter','none');
```

<img src="doc\Demo_3.1_arrows\Example_1.5.jpg" alt="Example_1.5" style="zoom:50%;" />

##### Demo - hide arrow body (parameter: 'is_plot_body'). 

```matlab

figure
xlim([-88 -62]);
ylim([22 48]);

hold on

hd1 = FUN_quiver_by_plotV2( x, y, u, v, 0, 'is_plot_body', false, 'head_length_pixel', 30, 'color','k' );
box on

title({'Example 1.6: "is_plot_body=false"','In this case, length of the arrow head must be','provided in parameter "head_length"'},'Interpreter','none');
```

<img src="doc\Demo_3.1_arrows\Example_1.6.jpg" alt="Example_1.6" style="zoom:50%;" />

##### Demo - perfect arrow head 

```matlab
figure
xlim([-0.2 1]*0.4);
ylim([-0.2 1]*1);
hold on

% this toolbox
hd1 = FUN_quiver_by_plotV2( 0, 0, 5, 5,0.050 );

% built-in matlab function
hd2 = quiver(0, 0, 5, 5 , 0.05, 'color','r', 'AutoScale', 'off','MaxHeadSize',15);

lid = legend( [hd1, hd2], 'FUN_quiver_by_plotV2', 'quiver');
set(lid,'Interpreter','none');
box on

title('Example 1.8: a "perfect" arrow','Interpreter','none');
```

<img src="doc\Demo_3.1_arrows\Example_1.8.jpg" alt="Example_1.8" style="zoom:50%;" />



### 3.2 Plot colorful arrows (`FUN_quiver_by_plotV2_cmap_patch`)



```matlab
[h1,h2, uu, vv, hu, hv] = FUN_quiver_by_plotV2( x, y, u, v)
[h1,h2, uu, vv, hu, hv] = FUN_quiver_by_plotV2( x, y, u, v, vel_scale)
[h1,h2, uu, vv, hu, hv] = FUN_quiver_by_plotV2( x, y, u, v, vel_scale, 'zval', zval)
[h1,h2, uu, vv, hu, hv] = FUN_quiver_by_plotV2( x, y, u, v, vel_scale, 'zval', zval, ...)
```

#### INPUT

* paramter 'zval' is used to determine the color. By default, `zval = sqrt( u.^2 + v.^2)`.
* Other parameters are identical to `FUN_quiver_by_plotV2`.



#### output 

* Identical to `FUN_quiver_by_plotV2`.

  

#### Demo

##### Demo - a simple case

In this case, the color is defined by the length of arrows, which is usually the magnitude of velocities.

```matlab
%% prepare data
x = 1:40;
y = 1:30;

[X,Y] = meshgrid(x,y);
X=X';
Y=Y';

u = sin( 2*pi/30 .* X );
v = cos( 2*pi/30 .* Y );

% plot
figure

arrow_scale = 0; % 0: the scale is determined by the script automatically.
FUN_quiver_by_plotV2_cmap_patch( x, y, u, v, arrow_scale );

grid on
box on
xlabel('x');
ylabel('y');
colorbar

title({'Example 1: a simple case', 'Color is defined by the magnitude of arrows'});


```

<img src="doc\Demo_3.2_colorful_arrows\Demo_2.1.jpg" alt="Demo_2.1" style="zoom:50%;" />

##### Demo - the color is defined by a given matrix

In this example, the arrow color is defined in parameter `zval`.

```matlab
%% prepare data
x = 1:40;
y = 1:30;

[X,Y] = meshgrid(x,y);
X=X';
Y=Y';

u = sin( 2*pi/30 .* X );
v = cos( 2*pi/30 .* Y );

%% plot
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
```

<img src="doc\Demo_3.2_colorful_arrows\Example_2.2.jpg" alt="Example_2.2" style="zoom:50%;" />

### 3.3 Work with mmap

#### 3.3.1 Simple arrows (`FUN_quiver_by_plotV2_mmap`)

Inputs are identical to `FUN_quiver_by_plotV2`.

It is recommended to call m_grid before `FUN_quiver_by_plotV2_mmap`

```matlab
%% prepare data

x = -1:0.2:1;
y = -1:0.2:1;

[X,Y] = meshgrid(x,y);
X=X';
Y=Y';

theta = atan2d(Y,X);

r = sqrt( X.^2 + Y.^2);

u =  sind( theta ) .* r / 100;
v = -cosd( theta ) .* r / 100;

x = x*10 - 75;
y = y*10 + 35;

%% plot
figure
hold on

m_proj('lambert', 'lon', [-88 -62], 'lat', [22 48] );
m_grid;

hd1 = FUN_quiver_by_plotV2_mmap( x, y, u, v, 0 );

m_coast('color','b')

title({'Example 3.1: Arrows on m_map'},'Interpreter','none');
```

<img src="doc/Demo_3.3_mmap/Example_3.1.jpg" alt="Example_3.1" style="zoom:50%;" />



### 3.3.2 Colorful arrows (`FUN_quiver_by_plotV2_cmap_patch_mmap`)

Inputs are identical to `FUN_quiver_by_plotV2_cmap_patch`.

It is recommended to call m_grid before `FUN_quiver_by_plotV2_mmap`

```matlab
%% prepare data

x = -1:0.2:1;
y = -1:0.2:1;

[X,Y] = meshgrid(x,y);
X=X';
Y=Y';

theta = atan2d(Y,X);

r = sqrt( X.^2 + Y.^2);

u =  sind( theta ) .* r / 100;
v = -cosd( theta ) .* r / 100;

x = x*10 - 75;
y = y*10 + 35;

%% plot
figure
hold on

m_proj('lambert', 'lon', [-88 -62], 'lat', [22 48] );
m_grid

FUN_quiver_by_plotV2_cmap_patch_mmap( x, y, u, v, 0 );
m_coast('color','b')
colorbar

title({'Example 3.2: Colorful arrows on m_map'},'Interpreter','none');
```

<img src="doc\Demo_3.3_mmap\Example_3.2.jpg" alt="Example_3.2" style="zoom:50%;" />

## 4. Vector Time Series

To plot vector time series, plot set `is_correct_angle` to `true`.

```matlab
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
```

<img src="doc/Demo_3.1b_vector_time_series/Example_1.1.jpg" alt="Example_1.1" style="zoom:50%;" />





## 5. Acknowledgements

Matlab script `plotboxpos`, which is available from https://github.com/kakearney/plotboxpos-pkg, is adopted to calculate the accurate axis ratio. 

