#include "planet.h"
#include "solver.h"

/*
Check if a planet is bound to the system or not
*/
void check_if_bound(double energy, double kinetic, double potential, std::string method)
{
    if (energy < 0.0)
    {
        printf("%f %f %f\n", energy, kinetic, potential);
        if (abs(2 * kinetic) == abs(potential))
        {
            printf("Circular orbit! \n");
        }
        std::cout << method;
        printf(": System is bound!\n");
    }
    else
    {
        printf("System is not bound. \n");
    }
}

/*
Run a simulaton of all planets in the solar system, including pluto
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
}

/*
Run various test of the program using just the Earth and the Sun
*/
void earth_sun_system(bool stability = true, bool timing = false, bool var_beta = false, bool escape = false)
{
    Solver solve;
    Planet sun(1.989e30, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    Planet earth(5.972e24, 1.0, 0.0, 0.0, 0.0, (2 * M_PI) / sun.YEARS, 0.0);
    solve.addPlanet(sun);
    solve.addPlanet(earth);

    if (stability)
    {
        //Stability tests
        double t[4] = {100, 1000, 10000, 100000};
        double kinetic, potential, energy;
        double init_kinetic, init_potential, init_energy;
        init_kinetic = earth.kineticEnergy();
        init_potential = earth.potentialEnergy(sun, solve.G, 0.0);
        init_energy = init_kinetic + init_potential;
        printf("Initial values: \n Kinetic: %f, Potential: %f, Total = %f \n", init_kinetic, init_potential, init_energy);

        for (int i = 0; i < 4; ++i)
        {
            std::string n = std::to_string(t[i]);
            solve.setResultsFolder("sun_earth/stability/verlet/step_" + n);
            solve.verletMethod(3, t[i], false);

            Planet temp = solve.getPlanet(1);
            kinetic = temp.kineticEnergy();
            potential = temp.potentialEnergy(solve.getPlanet(0), solve.G, 0.0);
            energy = kinetic + potential;

            check_if_bound(energy, kinetic, potential, "Verlet");

            solve.setResultsFolder("sun_earth/stability/euler/step_" + n);
            solve.forwardEulerMethod(3, t[i], false);

            Planet temp2 = solve.getPlanet(1);
            kinetic = temp2.kineticEnergy();
            potential = temp2.potentialEnergy(solve.getPlanet(0), solve.G, 0.0);
            energy = kinetic + potential;

            check_if_bound(energy, kinetic, potential, "Euler");

            solve.setResultsFolder("sun_earth/stability/euler_cromer/step_" + n);
            solve.eulerCromerMethod(3, t[i], false);
        }
    }

    if (timing)
    {
        //Timing test
        for (int i = 0; i < 10; ++i)
        {
            solve.setResultsFolder("verlet_timing");
            solve.verletMethod(3, 1000000, false);
            solve.setResultsFolder("euler_timing");
            solve.forwardEulerMethod(3, 1000000, false);
        }
    }

    if (var_beta)
    {
        Solver solve2;
        Planet sun(1.989e30, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        Planet earth(5.972e24, 1.0, 0.0, 0.0, 0.0, (5) / sun.YEARS, 0.0);
        solve2.addPlanet(sun);
        solve2.addPlanet(earth);
        double beta[7] = {3, 3.2, 3.4, 3.6, 3.8, 3.9, 4};
        for (int i = 0; i < 7; ++i)
        {
            std::string n = std::to_string(beta[i]);
            solve.setResultsFolder("sun_earth/variable_beta/" + n);
            solve.verletMethod(3, 10000.0, false, beta[i]);
            solve2.setResultsFolder("sun_earth/variable_beta_ellipse/" + n);
            solve2.verletMethod(3, 10000.0, false, beta[i]);
        }
    }

    if (escape)
    {
        Solver solve;
        Planet sun(1.989e30, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        Planet earth(5.972e24, 1.0, 0.0, 0.0, 0.0, (sqrt(2) * 2 * M_PI) / sun.YEARS, 0.0);
        solve.addPlanet(sun);
        solve.addPlanet(earth);
        solve.setResultsFolder("sun_earth/escape");
        double kinetic, potential, energy;

        Planet temp = solve.getPlanet(1);
        kinetic = temp.kineticEnergy();
        potential = temp.potentialEnergy(solve.getPlanet(0), solve.G, 0.0);
        energy = kinetic + potential;

        solve.verletMethod(10, 10000, false);
        Planet temp2 = solve.getPlanet(1);
        kinetic = temp2.kineticEnergy();
        potential = temp2.potentialEnergy(solve.getPlanet(0), solve.G, 0.0);
        energy = kinetic + potential;

        check_if_bound(energy, kinetic, potential, "Verlet");
    }
}

/*
Run simulation of the three body system, checking what happens when the mass of Jupiter increases by a 
factor of 10 and 1000
*/
void three_body_problem()
{
    Solver solver;
    Planet sun;
    Planet earth;
    Planet jupiter;

    solver.addPlanet(sun);
    solver.addPlanet(earth);
    solver.addPlanet(jupiter);

    double years = 10;

    solver.setReadFile("data/sun_earth_jupiter.txt");
    solver.setResultsFolder("three_body_problem/normal_mass");
    solver.verletMethod(years, 10000);

    solver.setReadFile("data/sun_earth_jupiter2.txt");
    solver.setResultsFolder("three_body_problem/normal_mass_x10");
    solver.verletMethod(years, 10000);

    solver.setReadFile("data/sun_earth_jupiter3.txt");
    solver.setResultsFolder("three_body_problem/normal_mass_x1000");
    solver.verletMethod(years, 10000);
}

/*
Simulate a system of just the Sun and mercury, using both a classical force of gravity, and a force
containing a relativistic correction term. 
*/
void perihelion()
{
    //Sun and Mercury system with mercury beginning at perihelion.
    //Using classical and relativistic corrected gravity
    //Perihelion at 0.3075 AU
    //v_m = 12.44 AU/yr = 0.03408219178082192 AU/day
    //Should be run with more than 1e7 values for good precision
    Solver solve;
    Planet sun;
    Planet mercury;
    solve.addPlanet(sun);
    solve.addPlanet(mercury);
    solve.setReadFile("data/sun_mercury.txt");

    solve.setResultsFolder("sun_mercury_cla");
    solve.verletMethod(100, 10000000.0);
    solve.setResultsFolder("sun_mercury_rel");
    solve.relVerletMethod(100, 10000000.0);
}

int main()
{
    //Uncomment to run desired simulation

    //earth_sun_system(false, false, false, true);
    //all_planets();
    //three_body_problem();
    //perihelion();

    return 0;
}