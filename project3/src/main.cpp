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

void all_planets()
{

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
    solver.setResultsFolder("solar_system");

    //just enough time for pluto to complete on orbit, ie 248 years
    solver.verletMethod(248, 10000.0);
    //solver.setResultsFolder("f_euler");
    //solver.forwardEulerMethod(5, 10000.0);
}

void earth_sun_system()
{
    //Sun and earth system with earth having a perfectly circular orbit.
    //Sun and earth system with earth having escape velocity.
    //circ orbit : v = 2 * M_PI
    //v_e = (sqrt(2) * 2 * M_PI)
    Solver solve;
    Planet sun(1.989e30, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    Planet earth(5.972e24, 1.0, 0.0, 0.0, 0.0, (2 * M_PI) / sun.YEARS, 0.0);
    solve.addPlanet(sun);
    solve.addPlanet(earth);
    solve.setReadFile("data/sun_earth.txt");

    solve.setResultsFolder("sun_earth");
    //solve.verletMethod(1, 10000.0, false);
    solve.forwardEulerMethod(50, 10000.0, false);
    solve.setResultsFolder("sun_earth_e_cromer");
    solve.eulerCromerMethod(50, 10000.0, false);
}

int main()
{
    /*
    NASA data is distance in AU, velocity in AU/day. days as earth days.
    I am using G = 4*pi^2, which assumes velocity 
    in AU/year, so convert velocity by multiplying by 365. This is done elsewhere in code
    Also have to measure mass in terms of solar mass, also done elsewhere
    */
    earth_sun_system();
    //all_planets();

    return 0;
}