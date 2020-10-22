#pragma once
#include "planet.h"
#include <iostream>
#include <filesystem>
#include <fstream>
#include <time.h>

namespace fs = std::filesystem;

/*
Class for solving the particular differential equation of particles(planets)
interacting in only a conservative gravitational field, using either the verlet velocity
method, or Euler's forward method(not implemented yet).
*/

class Solver
{
public:
    Solver();

    void addPlanet(const Planet &planet);
    void forwardEulerMethod(int endingTime, unsigned int N, bool read = true);
    void verletMethod(double endingTime, unsigned int N, bool read = true);
    //void specialForce();
    void readData(const std::string &readFile);
    void setResultsFolder(std::string folder);
    void setReadFile(const std::string &readFile);

private:
    double systemMass = 0;
    std::vector<Planet> planets;
    const double G = 4 * pow(M_PI, 2);
    //const double G = 39.478;
    double cmp[3];
    double cmv[3];
    double pos[3];
    double vel[3];
    double acc[3];
    double force[3];

    void initialisePlanets(bool read = true);
    void gravitationalForce(Planet current, double &Fx, double &Fy, double &Fz);
    double centerMassPosition(uint16_t index);
    double centerMassVelocity(uint16_t index);
    void printer(std::ofstream &out, double x, double y, double z);
    void prepareFolder();
    void init();

    void delete_dir_content(const fs::path &dir_path);
    bool checkIfEmpty(const std::string &path);
    bool checkIfExists(const std::string &path);
    std::string path;
    std::string dir_path;
    std::string folder = "results";
    std::string readFile = dir_path + "data/solar_system.txt";
};