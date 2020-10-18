#include "planet.h"
#include "solver.h"

//TODO:
/*
Implement unit test.
Simplify code
Solve all tasks
Remove surplus code
Structure better
Kinetic/potential energy?
*/

int main()
{
    /*
    NASA data is distance in AU, velocity in AU/day.
    I am using G = 39.478, which assumes velocity 
    in AU/year, so convert velocity by multiplying by 365. This is done elsewhere in code
    Also have to measure mass in terms of solar mass, also done elsewhere
    */

    Solver solver;
    Planet sun;
    Planet earth;
    Planet jupiter;
    Planet saturn;
    Planet venus;
    Planet mars;
    Planet mercury;
    Planet uranus;
    Planet neptune;
    Planet pluto;

    solver.addPlanet(sun);
    solver.addPlanet(earth);
    solver.addPlanet(jupiter);
    solver.addPlanet(saturn);
    solver.addPlanet(venus);
    solver.addPlanet(mars);
    solver.addPlanet(mercury);
    solver.addPlanet(uranus);
    solver.addPlanet(neptune);
    solver.addPlanet(pluto);

    //just enough time for pluto to complete on orbit, ie 248 years
    solver.verletMethod(248.0, 10000.0);

    //Sun and earth system with earth having a perfectly circular orbit.
    Solver solve2;
    Planet sun2(1.989e30, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    Planet earth2(5.972e24, 1.0, 0.0, 0.0, 0.0, 2 * M_PI / sun.YEARS, 0.0);
    solve2.addPlanet(sun2);
    solve2.addPlanet(earth2);
    solve2.setResultsFolder("data/sun_earth");
    solve2.verletMethod(5, 10000.0, false);

    return 0;
}