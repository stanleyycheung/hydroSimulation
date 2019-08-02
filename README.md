# hydroSimulation  
Written in Python 3  
Modules needed: ```numpy```, ```matplotlib```  
UROP project (July - Aug 2019) under Dr. Haworth in Imperial College Astrophysics Department  
## Description  
1D pseudo disc-wind model  
Also includes a general hydro solver  
Features: includes rotational motion,
### Instructions:  
Running ```solver.py``` will produce .dat files in corresponding output  
Running the corresponding ```movie.py``` will produce .png in corresponding movie folder
Use ffmpeg command in the ```movie.py``` file to produce .mp4  
### Test cases:  
* ```sodshock.py``` for sodshock case  
    *  includes energy equation and advection  
* ``` simulation_isothermal.py``` includes atmosphere test  
    * includes gravity, but not rotation  
    * ```atmosphere.py``` for initial atmosphere - doesn't work  

## Mainfiles (in order of importance)  
* **```simulation_isothermal.py```** isothermal hydro simulation - no energy advection  
    * fully documented  
* ```simulation.py``` full hydro simulation - includes rotation and energy advection and gravity  
* ```sodshock.py``` sod shock test - energy advection but no rotation or gravity  
* various ```plot.py``` files - plotting specific output.dat files  
* various ```movie.py``` files - taking .dat files and producing graphs for movie  
* /debug folder - if ```debug = True```, will output log in this folder
* various /output folders - where solver outputs .dat raw data files  
* various /movie folders - where ```movie.py``` output .png files  
* /FRIED_grid - examples of FRIED_grid plots  
* /examples - a bunch of example solvers from notes by V. Springel and C.P. Dullemond  

## Multimedia    
* ```moviesod.mp4``` and ```moviesod2.mp4``` - results of sod shock test  
* ```moviedisc.mp4``` - initial disc model (incorrect)  
* ```movieatmosphere.mp4``` - atmosphere test (incorrect)  
* ```moviedisc_iso.mp4``` - attempt at a psuedo wind model (unsure if correct)  
* ```moviedisc_iso_grav.mp4``` - more accurate testing of atmosphere case  


## Changelog  
* 19/07  
    * created repo  
* 22/07  
    * ```advection.py```, ```advection2.py```: works  
    * working on ```hydro1.py```: works(??) see movie.mp4  
    * started on ```hydro2.py```  
* 23/07  
    * ```hydro2.py``` works!!!  
    * testing sod shock case  
    * implimented energy equation (section 5.5)  
* 24/07  
    * implimented gravity (section 5.8)  
    * atmosphere test case in ```atmosphere.py``` - initial uniform density let to evolve  
    * implimented rotational force (rotational velocity is not being advected!!!)   
* 25/07  
    * working on ```discmodel.py```  
    * implimented disc BC (need to double check if right or not)  
    * rotational velocity is still not advecting  
    * time step dt seems to go smaller and smaller - freezing the simulation  
        * cause: disc BC implimentation  
* 26/07  
    * rescaled ```discmodel.py``` and ```moviedisc.py``` for disc problem  
    * freezing problem from ~~ghost cells~~ sharp boundary conditions  
        * need BC on velocity  
    * fixed dt problem by revising inner BC of u  
        * can maybe put back in ghost cells?  
    * created ```simulation.py``` - cleaner version of ```discmodel.py```  
* 29/07  
    * completed ```simulation.py```  
    * created ```simulation_isothermal.py``` - isothermal version of ```simulation.py```  
        * no energy equation  
* 30/07  
    * nothing works :'(  
    * recreated atmosphere test in disc model (```/~_grav.py``` files)  
* 31/07  
    * working on atmosphere test in disc  
* 01/08  
    * still working on atmosphere test  
    * updated documentation  

## Todo  
- [ ] use ```simulation_isothermal.py``` or ```simulation.py``` as a class to automate plotting and making movie into one file.  
Example (can be used in ```moviedisc_iso_grav.py```):  
```python
import simulation_isothermal  
simulation = discSimulation(500, 50000, 0, 2e4, 0.3, 'exp', 'flat', 'kep', 'mirror/free', 'van Leer', rotation=False)  
simulation.run(time_interval=50, debug=False)  
# rest of code to plot movie  
```  

- [ ] tidy up ```self._rotation``` implimentation (it is just a bunch of if states right now)  
- [ ] note that m_gas is hard coded into ```simulation_isothermal.py```- can be implimented in the class structure  
- [ ] physical constants can be included in class  
