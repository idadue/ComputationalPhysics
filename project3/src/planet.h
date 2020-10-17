#pragma once
#include <cmath>
#include <random>

//CONSIDER: IF G DOES NOT NEED TO BE CHANGED, THEN MAKE IT A PART OF THE PRIVATE VARS.

class Planet
{
public:
    Planet();
    Planet(double mass);
    Planet(double mass, double x, double y, double z, double vx, double vy, double vz);
    void assignID();

    void setMass(double mass);
    double distance(const Planet &otherPlanet);
    double gravitationalForce(const Planet &otherPlanet, const double G);
    double acceleration(const Planet &otherPlanet, const double G);
    double kineticEnergy();
    double potentialEnergy(const Planet &otherPlanet, const double G, double epsilon);

    double getMass() const;
    double getID() const;

    bool operator==(const Planet &rhs)
    {
        return (id == rhs.id) ? true : false;
    }

    void setPosition(double x, double y, double z);
    void setVelocity(double vx, double vy, double vz);

    double getPosition(int i) const;
    double getVelocity(int i) const;

    const double YEARS = 365.242199;
    const double solarMass = 1.989e30;

private:
    double id;
    double mass;
    double position[3];
    double velocity[3];
    double potential;
    double kinetic;
};