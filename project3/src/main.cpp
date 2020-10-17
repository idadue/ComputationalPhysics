#include "planet.h"
#include "solver.h"

//TODO:
/*
Implement unit test.
Simplify code
Solve all tasks
Add timing of algorithm/method
Remove surplus code
Structure better
Solve other todo's
Kinetic/potential energy?
*/

int main()
{
    /*
    NASA data is distance in AU, velocity in AU/day.
    I am using G = 4 * pow(M_PI, 2), which assumes velocity 
    in AU/year, so convert by multiplying by 365. This is done elsewhere in code
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

    //solver.readData();

    //just enough time for pluto to complete on orbit, ie 248 years
    double years = 248;
    solver.solarSystem(years, years * 365, 0);

    return 0;
}
