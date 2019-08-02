# hydroSimulation  
1D pseudo disc-wind model  
Also includes a general hydro solver  

## Changelog  
* 19/07  
    * created repo  
* 22/07  
    * '''advection.py''', '''advection2.py''': works  
    * working on '''hydro1.py''': works(??) see movie.mp4  
    * started on '''hydro2.py'''  
* 23/07  
    * '''hydro2.py''' works!!!  
    * testing sod shock case  
    * implimented energy equation (section 5.5)  
* 24/07  
    * implimented gravity (section 5.8)  
    * atmosphere test case in '''atmosphere.py''' - initial uniform density let to evolve  
    * implimented rotational force (rotational velocity is not being advected!!!)   
* 25/07  
    * working on '''discmodel.py'''  
    * implimented disc BC (need to double check if right or not)  
    * rotational velocity is still not advecting  
    * time step dt seems to go smaller and smaller - freezing the simulation  
        * cause: disc BC implimentation  
* 26/07  
    * rescaled '''discmodel.py''' and '''moviedisc.py''' for disc problem  
    * freezing problem from ~~ghost cells~~ sharp boundary conditions  
        * need BC on velocity  
    * fixed dt problem by revising inner BC of u  
        * can maybe put back in ghost cells?  
    * created '''simulation.py''' - cleaner version of '''discmodel.py'''  
* 29/07  
    * completed '''simulation.py'''  
    * created '''simulation_isothermal.py''' - isothermal version of '''simulation.py'''  
        * no energy equation  
* 30/07  
    * nothing works :'(  
    * recreated atmosphere test in disc model (_grav files)  
* 31/07  
    * working on atmosphere test in disc  
* 01/08  
    * still working on atmosphere test  
