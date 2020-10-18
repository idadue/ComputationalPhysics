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
    time_t start, finish;

    double cmp[3];
    double cmv[3];
    double pos[3];
    double vel[3];
    double acc[3];
    double force[3];

    if (read)
    {
        readData(readFile);
    }
    for (int i = 0; i < 3; i++)
    {
        cmp[i] = centerMassPosition(i);
        cmv[i] = centerMassVelocity(i);
    }
    prepareFolder();

    start = clock();
    while (elapsedTime <= endingTime)
    {
        for (long unsigned int i = 0; i < planets.size(); ++i)
        {
            std::ofstream out;
            std::string filename = path + "/planet_" + std::to_string(i) + ".txt";
            out.open(filename, std::ios::app);

            if (elapsedTime == 0.0)
            {
                for (int j = 0; j < 3; ++j)
                {
                    pos[j] = planets[i].getPosition(j) - cmp[j];
                    vel[j] = planets[i].getVelocity(j) - cmv[j];
                    planets[i].setPos(j, pos[j]);
                }
            }
            else
            {
                for (int j = 0; j < 3; ++j)
                {
                    pos[j] = planets[i].getPosition(j);
                    vel[j] = planets[i].getVelocity(j);
                }
            }
            gravitationalForce(planets[i], force[0], force[1], force[2]);

            for (int j = 0; j < 3; ++j)
            {
                acc[j] = force[j] / planets[i].getMass();
                pos[j] += vel[j] * dt + acc[j] * pow(dt, 2) * 0.5;
                planets[i].setPos(j, pos[j]);
            }

            gravitationalForce(planets[i], force[0], force[1], force[2]);
            printer(out, pos[0], pos[1], pos[2]);

            for (int j = 0; j < 3; ++j)
            {
                acc[j] += force[j] / planets[i].getMass();
                vel[j] += acc[j] * (0.5 * dt);
                planets[i].setVel(j, vel[j]);
            }
        }
        elapsedTime += dt;
    }
    finish = clock();
    printf("Execution time is : %f seconds. \n", double(finish - start) / double(CLOCKS_PER_SEC));
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

void Solver::prepareFolder()
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

void Solver::setResultsFolder(std::string folder)
{
    this->folder = folder;
}

void Solver::setReadFile(const std::string &readFile)
{
    this->readFile = readFile;
}

void Solver::readData(const std::string &readFile)
{
    std::ifstream inData;
    inData.open(dir_path + readFile);
    double mass, x, y, z, vx, vy, vz;
    for (auto &it : planets)
    {
        inData >> mass;
        inData >> x >> y >> z;
        inData >> vx >> vy >> vz;
        it.setMass(mass);
        systemMass += it.getMass();
        it.setPosition(x, y, z);
        it.setVelocity(it.YEARS * vx, it.YEARS * vy, it.YEARS * vz);
    }
}