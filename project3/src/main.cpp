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

/*
Sun and earth system:
    Choose sun position to be at origin ie (0, 0, 0)
    Task 3c:
        Find initial value for velocity that gives a circular orbit
        Test stability of algorithm as function of different time steps
        Make plot of results
        Check that kinetic and potential energy is conserved and explain
            why they should be
        Discuss differences between Verlet and euler algorithm
        Compare timing of both methods
        Consider number of flops involved for both 
    
    Task 3d:
        Discuss your results with
            circular and elliptical orbits for the Earth-Sun system.
        Make plots?
    
    Task 3e:
        Implement force calculation with variable beta, ie pow(r, beta)
            where beta = [2,3]
        What happens to system when beta -> 3?
        Is the rest purely analytic?
    
    Task 3f:
        Find by trial and error the escape velocity.
        Analytical escape velocity is sqrt(2)*2*pi

Three body problem:
    Task 3g:
        How much does Jupiter alter Earth's motion?
        Make a plot(include sun?)
        Discuss the stability of the solutions using the verlet method
        Repeat above subtasks with mass 
            of Jupiter changed by a factor of 10 and 1000

Misc:    
    Task 3h:
        Simulate three body system with sun not at origin and a fixed center
            of mass and compare with results from previous task.
        Extend to all planets \\
        Optional: Add moons?
    
    Task 3i:
        Add a general relativistic correction to the force
        Run simulation over one century(100 yr) with just Mercury and Sun
            with Mercury at perihilion(closest to Sun) on x-axis.
        Find perihilion angle theta_p = y_p/x_p 
        Speed at perihelion is 12.44 AU/yr
        Distance is 0.3075 AU
        Make sure that the time resolution used in simulation is sufficient
            by checking that the perihelion precession is at least a few 
                orders of magnitude smaller than the observed 
                    perihelion precession of Mercury.
        Can the observed perihelion precession of Mercury be explained 
            by the general theory of relativity?



*/

//TODO: Turn N into double

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

    //Stability tests
    double t[4] = {100, 1000, 10000, 100000};
    double energy;
    for (int i = 0; i < 4; ++i)
    {
        std::string n = std::to_string(t[i]);
        solve.setResultsFolder("sun_earth/stability/verlet/step_" + n);
        solve.verletMethod(3, t[i], false);
        energy = earth.kineticEnergy() + earth.potentialEnergy(sun, 4 * pow(M_PI, 2), 0.0);
        printf("Energy = %f \n", energy);
        solve.setResultsFolder("sun_earth/stability/euler/step_" + n);
        solve.forwardEulerMethod(3, t[i], false);
        energy = earth.kineticEnergy() + earth.potentialEnergy(sun, 4 * pow(M_PI, 2), 0.0);
        printf("Energy = %f \n", energy);
        printf("pos = %f \n", earth.getPosition(0));
    }

    //Timing test
    /*for (int i = 0; i < 10; ++i)
    {
        solve.setResultsFolder("verlet_timing");
        solve.verletMethod(3, 1000000, false);
        solve.setResultsFolder("euler_timing");
        solve.forwardEulerMethod(3, 1000000, false);
    }*/

    /*
    double beta[7] = {3, 3.2, 3.4, 3.6, 3.8, 3.9, 4};
    for (int i = 0; i < 7; ++i)
    {
        std::string n = std::to_string(beta[i]);
        solve.setResultsFolder("sun_earth/variable_beta/" + n);
        solve.variableBeta(1, 10000.0, beta[i]);
    }*/
    //solve.setResultsFolder("sun_earth");
    //solve.verletMethod(1, 10000.0, false);
    //solve.forwardEulerMethod(50, 10000.0, false);
    /*solve.setResultsFolder("sun_earth_e_cromer");
    for (int i = 0; i < 10; i++)
    {
        solve.eulerCromerMethod(5, 100000.0, false);
    }*/
}

void sun_merkur()
{
    Solver solve;
    Planet sun;
    Planet mercury;
    solve.addPlanet(sun);
    solve.addPlanet(mercury);
    solve.setReadFile("data/sun_mercury.txt");

    solve.setResultsFolder("sun_mercury");
    solve.verletMethod(3, 10000.0);
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
    //sun_merkur();
    //all_planets();

    return 0;
}