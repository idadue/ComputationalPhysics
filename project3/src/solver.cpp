#include "solver.h"

Solver::Solver()
{
    //Finds current path, but deletes lasts 3 chars, which in this case means
    //returning to the ./project3 folder
    //Note: This only works if current folder is src or other 3 letter length word
    path = fs::current_path();
    int str_len = path.length();
    path = path.erase(str_len - 3, 3);
    dir_path = path;
}

void Solver::forwardEulerMethod(int endingTime, unsigned int N)
{
    /*
    Pass
    */
    double dt = double(endingTime) / double(N);
    double elapsedTime = 0;
}

/*
Solve diff equation using the Verlet velocity method.
*/
void Solver::verletMethod(int endingTime, unsigned int N, bool read)
{
    double dt = endingTime / double(N);
    double elapsedTime = 0;

    double x, y, z, vx, vy, vz;
    double Fx, Fy, Fz;
    double ax, ay, az, ax_dt, ay_dt, az_dt;

    if (read)
    {
        readData();
    }
    prepareFolder("data");

    while (elapsedTime <= endingTime)
    {
        for (long unsigned int i = 0; i < planets.size(); ++i)
        {
            std::ofstream out;
            std::string filename = path + "/planet_" + std::to_string(i) + ".txt";
            out.open(filename, std::ios::app);

            if (elapsedTime == 0.0)
            {
                x = planets[i].getPosition(0) - centerMassPosition(0);
                y = planets[i].getPosition(1) - centerMassPosition(1);
                z = planets[i].getPosition(2) - centerMassPosition(2);

                vx = planets[i].getVelocity(0) - centerMassPosition(0);
                vy = planets[i].getVelocity(1) - centerMassPosition(1);
                vz = planets[i].getVelocity(2) - centerMassPosition(2);
                printer(out, x, y, z);
            }
            else
            {
                x = planets[i].getPosition(0);
                y = planets[i].getPosition(1);
                z = planets[i].getPosition(2);

                vx = planets[i].getVelocity(0);
                vy = planets[i].getVelocity(1);
                vz = planets[i].getVelocity(2);
            }

            gravitationalForce(planets[i], Fx, Fy, Fz);

            ax = Fx / planets[i].getMass();
            ay = Fy / planets[i].getMass();
            az = Fz / planets[i].getMass();

            x = x + vx * dt + ax * pow(dt, 2) * 0.5;
            y = y + vy * dt + ay * pow(dt, 2) * 0.5;
            z = z + vz * dt + az * pow(dt, 2) * 0.5;
            planets[i].setPosition(x, y, z);

            gravitationalForce(planets[i], Fx, Fy, Fz);
            printer(out, x, y, z);

            ax_dt = Fx / planets[i].getMass();
            ay_dt = Fy / planets[i].getMass();
            az_dt = Fz / planets[i].getMass();

            vx = vx + (ax_dt + ax) * (0.5 * dt);
            vy = vy + (ay_dt + ay) * (0.5 * dt);
            vz = vz + (az_dt + az) * (0.5 * dt);
            planets[i].setVelocity(vx, vy, vz);

            x = y = z = 0;
            vx = vy = vz = 0;
            Fx = Fy = Fz = 0;
        }
        elapsedTime += dt;
    }
}

void Solver::addPlanet(const Planet &planet)
{
    planets.push_back(planet);
    systemMass += planet.getMass();
}

void Solver::gravitationalForce(Planet current, double &Fx, double &Fy, double &Fz)
{
    double x, y, z;
    x = y = z = 0;
    for (auto it = planets.begin(); it != planets.end(); ++it)
    {
        if (it->getID() != current.getID())
        {
            x += (it->getMass() * (current.getPosition(0) - it->getPosition(0)) / pow(current.distance(*it), 3));
            y += (it->getMass() * (current.getPosition(1) - it->getPosition(1)) / pow(current.distance(*it), 3));
            z += (it->getMass() * (current.getPosition(2) - it->getPosition(2)) / pow(current.distance(*it), 3));
        }
    }
    Fx = -G * current.getMass() * x;
    Fy = -G * current.getMass() * y;
    Fz = -G * current.getMass() * z;
}

void Solver::printer(std::ofstream &out, double x, double y, double z)
{
    out << x << ", " << y << ", " << z << std::endl;
}

double Solver::centerMassPosition(uint16_t index)
{
    double res = 0;
    for (auto it = planets.begin(); it != planets.end(); ++it)
    {
        res += it->getMass() * it->getPosition(index);
    }
    return res / systemMass;
}

double Solver::centerMassVelocity(uint16_t index)
{
    double res = 0;
    for (auto it = planets.begin(); it != planets.end(); ++it)
    {
        res += it->getMass() * it->getVelocity(index);
    }
    return res / systemMass;
}

void Solver::delete_dir_content(const fs::path &dir_path)
{
    for (auto &path : fs::directory_iterator(dir_path))
    {
        fs::remove_all(path);
    }
}

bool Solver::checkIfEmpty(const std::string &path)
{
    return fs::is_empty(path);
}

bool Solver::checkIfExists(const std::string &path)
{
    return fs::exists(path);
}

void Solver::prepareFolder(const std::string &folder)
{
    path = path + folder;
    if (!checkIfExists(path))
    {
        fs::create_directory(path);
    }
    if (!checkIfEmpty(path))
    {
        fs::path p = path;
        delete_dir_content(p);
    }
}

void Solver::readData()
{
    std::ifstream inData;
    inData.open(dir_path + "NASA/indata.txt");
    double mass, x, y, z, vx, vy, vz;
    for (auto &it : planets)
    {
        inData >> mass;
        inData >> x >> y >> z;
        inData >> vx >> vy >> vz;
        it.setMass(mass);
        systemMass += mass;
        it.setPosition(x, y, z);
        it.setVelocity(it.YEARS * vx, it.YEARS * vy, it.YEARS * vz);
    }
}