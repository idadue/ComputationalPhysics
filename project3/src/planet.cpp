#include "planet.h"

Planet::Planet() : mass(0)
{
    for (int i = 0; i < 3; i++)
    {
        position[i] = 0;
        velocity[i] = 0;
        initPos[i] = 0;
        initVel[i] = 0;
    }
    assignID();
}

Planet::Planet(double mass) : mass{mass / solarMass}
{
    assignID();
}

Planet::Planet(double mass, double x, double y, double z, double vx, double vy, double vz) : mass{mass / solarMass}
{
    initPos[0] = x, initVel[0] = YEARS * vx;
    initPos[1] = y, initVel[1] = YEARS * vy;
    initPos[2] = z, initVel[2] = YEARS * vz;

    for (int i = 0; i < 3; i++)
    {
        position[i] = initPos[i];
        velocity[i] = initVel[i];
    }

    assignID();
}

/*Assign an unqiue id to instance of class. This is overkill for planets with many differences, but
adds the option for the program to work with several particles of same size, mass, pos/vel or all. */
void Planet::assignID()
{
    std::random_device r;
    std::uniform_real_distribution<double> number(0, 100000);
    std::default_random_engine engine{r()};

    id = number(engine);
}

void Planet::setMass(double mass)
{
    this->mass = mass / solarMass;
}

double Planet::distance(const Planet &otherPlanet) const
{
    double x = position[0] - otherPlanet.position[0];
    double y = position[1] - otherPlanet.position[1];
    double z = position[2] - otherPlanet.position[2];

    return sqrt(x * x + y * y + z * z);
}

double Planet::gravitationalForce(const Planet &otherPlanet, const double G, const double beta)
{
    double r = distance(otherPlanet);
    double otherMass = otherPlanet.mass;
    //return (r != 0) ? G * mass * otherMass / pow(r, 2) : 0;
    return G * mass * otherMass / pow(r, beta);
}

void Planet::componentGravitationalForce(const Planet &other, double &Fx, double &Fy, double &Fz, const double G, const double beta)
{
    double x, y, z;

    x = this->getPosition(0) - other.getPosition(0);
    y = this->getPosition(1) - other.getPosition(1);
    z = this->getPosition(2) - other.getPosition(2);

    Fx = -G * this->getMass() * other.getMass() * x / pow(distance(other), beta);
    Fy = -G * this->getMass() * other.getMass() * y / pow(distance(other), beta);
    Fz = -G * this->getMass() * other.getMass() * z / pow(distance(other), beta);
}

double Planet::acceleration(const Planet &otherPlanet, const double G, const double beta)
{
    double force = gravitationalForce(otherPlanet, G, beta);
    return (force != 0) ? force / mass : 0;
}

double Planet::kineticEnergy()
{
    double vel2 = (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);
    return 0.5 * mass * vel2;
}

double Planet::potentialEnergy(const Planet &otherPlanet, const double G, double epsilon)
{
    double r = distance(otherPlanet);
    if (epsilon == 0.0)
        return -G * mass * otherPlanet.mass / r;
    else
        return (G * mass * otherPlanet.mass / epsilon) * (atan(r / epsilon) - (0.5 * M_PI));
}

double Planet::getMass() const
{
    return mass;
}

double Planet::getID() const
{
    return id;
}

void Planet::setPosition(double x, double y, double z)
{
    position[0] = x;
    position[1] = y;
    position[2] = z;
}

void Planet::setVelocity(double vx, double vy, double vz)
{
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
}

void Planet::setInitialPosition(double x, double y, double z)
{
    initPos[0] = x;
    initPos[1] = y;
    initPos[2] = z;
}

void Planet::setInitialVelocity(double vx, double vy, double vz)
{
    initVel[0] = vx;
    initVel[1] = vy;
    initVel[2] = vz;
}

void Planet::setPos(int index, double pos)
{
    position[index] = pos;
}

void Planet::setVel(int index, double vel)
{
    velocity[index] = vel;
}

double Planet::getPosition(int i) const
{
    return position[i];
}
double Planet::getVelocity(int i) const
{
    return velocity[i];
}

double Planet::getInitPos(int i) const
{
    return initPos[i];
}

double Planet::getInitVel(int i) const
{
    return initVel[i];
}