# VisAna
Scripts and functions for reading, visualizing and analyzing simulation data.
This package contains the reader for all kinds of output data from SWMF.
It also contains the series of analysis on Ganymede simulations.
A new version with more functionalities has been implemented with Julia.
Checkout [Batsrus.jl](https://github.com/henry2004y/Batsrus.jl/) and
[VisAnaJulia](https://github.com/henry2004y/VisAnaJulia) for more.

## Getting Started

### Prerequisites

Recent versions of MATLAB.

### Examples

Log data
```
filename = 'logTest.log';
[filehead,logdata] = read_data(filename,'filetype','log');

% Test of reading flyby trajectory data & compare with simulation outputs
[filehead,data] = read_log_data('Galileo_G8_flyby_MAG.dat');

time = datetime(data(:,1:6));
xyz  = data(:,7:9);
B_obs= data(:,10:12);
```

2D single snapshot
```
filename = '~/SWMF/SWMF/GM/BATSRUS/run_test/GM/IO2/y*';
[filehead,data] = read_data(filename);
plot_data(data.file1,filehead,'p','plotmode','contbar')
```

2D multiple snapshots
```
npict = 61; % index of snapshot to be read
filename = 'z.outs';
[filehead,data] = read_data(filename,'npict',npict);
plot_data(data.file1,filehead(1),'p','plotmode','grid');
func = 'p';
plotmode='contbar';
plot_data(data.file1,filehead(1),func,'plotmode',plotmode,'plotinterval',0.2);
plotrange = [-4 4 -4 4];
plot_data(data.file1,filehead(1),func...
   ,'plotmode','trimesh','plotrange',plotrange)

plot_data(data.file1,filehead.file1,func...
   ,'plotmode','trisurf','plotrange',plotrange)

plot_data(data.file1,filehead.file1,func...
   ,'plotmode','meshbar','plotrange',plotrange,'plotinterval',0.05)


filename='y=0*.out';
[filehead,data] = read_data(filename);

func = 'jy ux;uz';
plotmode = 'contbar streamover';
plotrange = [-8 8 -8 8];

plot_data(data.file1,filehead(1),func...
   ,'plotmode',plotmode,'plotrange',plotrange,'plotinterval',0.05)

rectangle('Position',[-.5,-.5,1,1],'Curvature',[1,1]...
   ,'FaceColor',[.6 .6 .6])
viscircles([0 0],1,'color',[.4 .2 1]);

filename='y=0_var_1_n00060000.out';
[filehead,data] = read_data(filename,'npict',1,'verbose',false);

func='jx ux;uy';
plotmode='contbar streamover';
plotrange = [-4 4 -4 4];
plot_data(data.file1,filehead(1),func,'plotmode',plotmode,...
   'plotrange',plotrange,'plotinterval',0.05)
hold on

rectangle('Position',[-2.8 -3 (-1.125+2.8) 6],'EdgeColor','r',...
  'LineWidth',2)

rectangle('Position',[-.5,-.5,1,1],'Curvature',[1,1]...
   ,'FaceColor',[.6 .6 .6]);
viscircles([0 0],1,'color',[.4 .2 1]);
```

Test of animation (not fully implemented)
```
func = 'P';
plotmode = 'contbar';
plotrange = [-4 4 -4 4];

animate_data(filename,func,'plotmode',plotmode,'plotrange',plotrange,...
   'firstpict',1,'dpict',1,'npictmax',61,...
   'plotinterval',0.05,'savemovie',true)
```

3D box output, Cartesian coordinates
```
filename = 'box*.outs';
[filehead,data] = read_data(filename);
```

3D IPIC3D output
```
filename = 'y=0_fluid_region0_1_n0_35781.outs';
[filehead,data] = read_data(filename,'npict',60);
plot_data(data.file1,filehead,'ps0','plotmode','contbar')
```

## Authors

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* All the nice guys who share their codes


