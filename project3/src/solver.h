#pragma once
#include "planet.h"
#include <iostream>
#include <filesystem>
#include <fstream>
#include <time.h>
#include <iomanip>

namespace fs = std::filesystem;

/*
Class for solving the particular differential equation of particles(planets)
interacting in only a conservative gravitational field, using either the verlet velocity
method, Euler's forward method or the Euler-Cromer method.
*/

class Solver
{
public:
    Solver();

    const double G = 4 * pow(M_PI, 2);
    void addPlanet(const Planet &planet);

    //Solver methods

    void forwardEulerMethod(int endingTime, double N, bool read = true);
    void eulerCromerMethod(int endingTime, double N, bool read = true);
    void verletMethod(double endingTime, double N, bool read = true, double beta = 3);
    void relVerletMethod(double endingTime, double N, bool read = true);

    void readData(const std::string &readFile);
    void setResultsFolder(std::string folder);
    void setReadFile(const std::string &readFile);
    Planet getPlanet(int index);

private:
    double systemMass = 0;
    std::vector<Planet> planets;

    double cmp[3], cmv[3];
    double pos[3], vel[3], acc[3];
    double force[3];
    double ang[3];

    double lightSpeed = 173.1;

    time_t start, finish;

    void initialisePlanets(bool read = true);

    void gravitationalForce(const Planet &current, double &Fx, double &Fy, double &Fz, double beta = 3);

    void gravitationalForceRel(Planet current, double &Fx, double &Fy, double &Fz);

    double centerMassPosition(uint16_t index);
    double centerMassVelocity(uint16_t index);

    void printer(std::ofstream &out, double x, double y, double z);
    std::ofstream outputStream(int i);

    void prepareFolder();
    std::string dir_path;
    std::string folder = "results/default";
    std::string readFile = dir_path + "data/solar_system.txt";
};
