#pragma once
#include <cmath>
#include <random>

/*
Class containing the information associated with a Planet, such as mass, position and velocity.
*/

class Planet
{
public:
    Planet();
    Planet(double mass);
    Planet(double mass, double x, double y, double z, double vx, double vy, double vz);
    void assignID();

    void setMass(double mass);
    double distance(const Planet &otherPlanet) const;
    double kineticEnergy();
    double potentialEnergy(const Planet &otherPlanet, const double G, double epsilon);

    double getMass() const;
    double getID() const;

    void setPosition(double x, double y, double z);
    void setVelocity(double vx, double vy, double vz);
    void setInitialPosition(double x, double y, double z);
    void setInitialVelocity(double vx, double vy, double vz);

    void setPos(int index, double pos);
    void setVel(int index, double vel);

    double getPosition(int i) const;
    double getVelocity(int i) const;

    double getInitPos(int i) const;
    double getInitVel(int i) const;

    const double YEARS = 365.242199;
    const double solarMass = 1.989e30;

private:
    double id;
    double mass;
    double position[3];
    double velocity[3];
    double potential;
    double kinetic;
    double initPos[3];
    double initVel[3];
};