#pragma once
#include "planet.h"
#include <cmath>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
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
    void forwardEulerMethod(int endingTime, unsigned int N);
    void verletMethod(int endingTime, unsigned int N, bool read = true);
    void readData();

private:
    double systemMass = 0;
    std::vector<Planet> planets;
    const double G = 4 * pow(M_PI, 2);

    void gravitationalForce(Planet current, double &Fx, double &Fy, double &Fz);
    double centerMassPosition(uint16_t index);
    double centerMassVelocity(uint16_t index);
    void printer(std::ofstream &out, double x, double y, double z);
    void prepareFolder(const std::string &folder);

    void delete_dir_content(const fs::path &dir_path);
    bool checkIfEmpty(const std::string &path);
    bool checkIfExists(const std::string &path);
    std::string path;
    std::string dir_path;
};