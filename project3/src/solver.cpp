#include "solver.h"

Solver::Solver()
{
    //Finds current path, but deletes lasts 3 chars, which in this case means
    //returning to the ./project3 folder
    //Note: This only works if current folder is src or other 3 letter length word
    std::string path = fs::current_path();
    int str_len = path.length();
    path = path.erase(str_len - 3, 3);
    dir_path = path;
}

/*
Add a Planet object to the planets vector
*/
void Solver::addPlanet(const Planet &planet)
{
    planets.push_back(planet);
    //systemMass += planet.getMass();
}

/*
Initialise the system with the option to read in data from file
Sets all planet inital positions and calculates the center of mass position and
center of mass velocity.
*/
void Solver::initialisePlanets(bool read)
{
    if (read)
        readData(readFile);

    prepareFolder();

    for (int i = 0; i < 3; i++)
    {
        cmp[i] = centerMassPosition(i);
        cmv[i] = centerMassVelocity(i);
        std::cout << "cmp before = " << centerMassPosition(i) << std::endl;
        std::cout << "cmv before = " << centerMassVelocity(i) << std::endl;

    }

    int i = 0;
    for (auto &it : planets)
    {
        std::ofstream out = outputStream(i);
        it.setPosition(it.getInitPos(0) - cmp[0], it.getInitPos(1) - cmp[1], it.getInitPos(2) - cmp[2]);
        it.setVelocity(it.getInitVel(0) - cmv[0], it.getInitVel(1) - cmv[1], it.getInitVel(2) - cmv[2]);
        printer(out, it.getPosition(0), it.getPosition(1), it.getPosition(2));
        ++i;
    }
    for (int i = 0; i < 3; i++)
    {
        std::cout << "cmp after = " << centerMassPosition(i) << std::endl;
        std::cout << "cmv after = " << centerMassVelocity(i) << std::endl;
    }
}

/*
Solve diff equation using forward euler method
*/
void Solver::forwardEulerMethod(int endingTime, double N, bool read)
{
    double dt = double(endingTime) / double(N);
    double elapsedTime = 0;

    initialisePlanets(read);

    start = clock();
    while (elapsedTime <= endingTime)
    {
        for (long unsigned int i = 0; i < planets.size(); ++i)
        {
            std::ofstream out = outputStream(i);
            gravitationalForce(planets[i], force[0], force[1], force[2]);

            for (int j = 0; j < 3; ++j)
            {
                acc[j] = force[j] / planets[i].getMass();
                pos[j] = planets[i].getPosition(j) + planets[i].getVelocity(j) * dt;
                vel[j] = planets[i].getVelocity(j) + acc[j] * dt;

                planets[i].setPos(j, pos[j]);
                planets[i].setVel(j, vel[j]);
            }
            printer(out, pos[0], pos[1], pos[2]);
        }
        elapsedTime += dt;
    }
    finish = clock();
    printf("Execution time for forward euler is: %f seconds. \n", double(finish - start) / double(CLOCKS_PER_SEC));
}

/*
Solve diff equation using euler-cromer method
*/
void Solver::eulerCromerMethod(int endingTime, double N, bool read)
{
    double dt = double(endingTime) / double(N);
    double elapsedTime = 0;

    initialisePlanets(read);

    start = clock();
    while (elapsedTime <= endingTime)
    {
        for (long unsigned int i = 0; i < planets.size(); ++i)
        {
            std::ofstream out = outputStream(i);
            gravitationalForce(planets[i], force[0], force[1], force[2]);

            for (int j = 0; j < 3; ++j)
            {
                acc[j] = force[j] / planets[i].getMass();
                vel[j] = planets[i].getVelocity(j) + acc[j] * dt;
                pos[j] = planets[i].getPosition(j) + vel[j] * dt;

                planets[i].setPos(j, pos[j]);
                planets[i].setVel(j, vel[j]);
            }
            printer(out, pos[0], pos[1], pos[2]);
        }
        elapsedTime += dt;
    }
    finish = clock();
    printf("Execution time for euler-cromer is: %f seconds. \n", double(finish - start) / double(CLOCKS_PER_SEC));
}

/*
Solve diff equation using the Verlet velocity method.
*/

void Solver::verletMethod(double endingTime, double N, bool read, double beta)

{
    double dt = endingTime / double(N);
    double elapsedTime = 0.0;

    initialisePlanets(read);
    std::cout << centerMassVelocity(0) << ", " << centerMassVelocity(1) << ", " << centerMassVelocity(2) << std::endl;
    std::cout << systemMass << std::endl;
    std::cout << planets[0].getVelocity(1) << std::endl;
    std::cout << planets[1].getVelocity(1) << std::endl;
    start = clock();
    while (elapsedTime <= endingTime)
    {
        int i = 0;
        for (auto &it : planets)
        {
            std::ofstream out = outputStream(i);

            gravitationalForce(it, force[0], force[1], force[2], beta);


            for (int j = 0; j < 3; ++j)
            {
                acc[j] = force[j] / it.getMass();
                pos[j] = it.getPosition(j) + it.getVelocity(j) * dt + acc[j] * pow(dt, 2) * 0.5;
                it.setPos(j, pos[j]);
            }

            gravitationalForce(it, force[0], force[1], force[2], beta);

            for (int j = 0; j < 3; ++j)
            {
                acc[j] += force[j] / it.getMass();
                vel[j] = it.getVelocity(j) + acc[j] * (0.5 * dt);
                it.setVel(j, vel[j]);
            }
            printer(out, pos[0], pos[1], pos[2]);
            ++i;
        }
        elapsedTime += dt;
    }
    finish = clock();
    printf("Execution time is : %f seconds. \n", double(finish - start) / double(CLOCKS_PER_SEC));
    systemMass = 0;
}


/*
Solve diff equation using the Verlet velocity method and a relativistic force model
*/
void Solver::relVerletMethod(double endingTime, double N, bool read)
{
    double dt = endingTime / double(N);
    double elapsedTime = 0.0;

    initialisePlanets(read);
    start = clock();
    while (elapsedTime <= endingTime)
    {
        int i = 0;
        for (auto &it : planets)
        {
            std::ofstream out = outputStream(i);
            gravitationalForceRel(it, force[0], force[1], force[2]);

            for (int j = 0; j < 3; ++j)
            {
                acc[j] = force[j] / it.getMass();
                pos[j] = it.getPosition(j) + it.getVelocity(j) * dt + acc[j] * pow(dt, 2) * 0.5;
                it.setPos(j, pos[j]);
            }
            gravitationalForceRel(it, force[0], force[1], force[2]);

            for (int j = 0; j < 3; ++j)
            {
                acc[j] += force[j] / it.getMass();
                vel[j] = it.getVelocity(j) + acc[j] * (0.5 * dt);
                it.setVel(j, vel[j]);
            }
            printer(out, pos[0], pos[1], pos[2]);
            ++i;
        }
        elapsedTime += dt;
    }
    finish = clock();
    printf("Execution time for Verlet with relativistic force is : %f seconds. \n", double(finish - start) / double(CLOCKS_PER_SEC));
}

/*
Find the gravitational force acting on a Planet in the system, in 3 dimensions
*/
void Solver::gravitationalForce(const Planet &current, double &Fx, double &Fy, double &Fz, double beta)
{
    double x, y, z;
    x = y = z = 0;
    for (auto it = planets.begin(); it != planets.end(); ++it)
    {
        if (it->getID() != current.getID())
        {
            x += (it->getMass() * (current.getPosition(0) - it->getPosition(0)) / (pow(current.distance(*it), beta)));
            y += (it->getMass() * (current.getPosition(1) - it->getPosition(1)) / (pow(current.distance(*it), beta)));
            z += (it->getMass() * (current.getPosition(2) - it->getPosition(2)) / (pow(current.distance(*it), beta)));
        }
    }
    Fx = -G * current.getMass() * x;
    Fy = -G * current.getMass() * y;
    Fz = -G * current.getMass() * z;
}

/*
Find the gravitational force acting on a Planet in the system, in 3 dimensions, with a relativistic correction
*/
void Solver::gravitationalForceRel(Planet current, double &Fx, double &Fy, double &Fz)

{
    double x, y, z;
    x = y = z = 0;
    double relCorrect;
    for (auto it = planets.begin(); it != planets.end(); ++it)
    {
        if (it->getID() != current.getID())
        {
            double xpos, ypos, zpos;
            xpos = current.getPosition(0) - it->getPosition(0);
            ypos = current.getPosition(1) - it->getPosition(1);
            zpos = current.getPosition(2) - it->getPosition(2);

            ang[0] = ypos * current.getVelocity(2) - zpos * current.getVelocity(1);
            ang[1] = zpos * current.getVelocity(0) - xpos * current.getVelocity(2);
            ang[2] = xpos * current.getVelocity(1) - ypos * current.getVelocity(0);

            double lSquare = 0;
            for (int i = 0; i < 3; ++i)
            {
                lSquare += pow(ang[i], 2);
            }

            relCorrect = 1 + ((3 * lSquare) / (pow(current.distance(*it), 2) * pow(lightSpeed, 2)));

            x += (it->getMass() * (current.getPosition(0) - it->getPosition(0)) / (pow(current.distance(*it), 3))) * relCorrect;
            y += (it->getMass() * (current.getPosition(1) - it->getPosition(1)) / (pow(current.distance(*it), 3))) * relCorrect;
            z += (it->getMass() * (current.getPosition(2) - it->getPosition(2)) / (pow(current.distance(*it), 3))) * relCorrect;
        }
    }
    Fx = -G * current.getMass() * x;
    Fy = -G * current.getMass() * y;
    Fz = -G * current.getMass() * z;
}


/*
Method that writes to a file by appending
*/
void Solver::printer(std::ofstream &out, double x, double y, double z)
{
    out << std::fixed << std::setprecision(15) << x << ", " << y << ", " << z << std::endl;
}

/*
Create an output stream for a file
*/
std::ofstream Solver::outputStream(int i)
{
    std::ofstream out;
    std::string filename = dir_path + folder + "/planet_" + std::to_string(i) + ".txt";
    out.open(filename, std::ios::app);
    return out;
}

/*
Calculate the center of mass postition of the system
*/
double Solver::centerMassPosition(uint16_t index)
{
    double res = 0;
    for (auto it = planets.begin(); it != planets.end(); ++it)
    {
        res += it->getMass() * it->getInitPos(index);
    }
    return res / systemMass;
}

/*
Calculate the center of mass velocity of the system
*/
double Solver::centerMassVelocity(uint16_t index)
{
    double res = 0;
    for (auto it = planets.begin(); it != planets.end(); ++it)
    {
        res += it->getMass() * it->getInitVel(index);
    }
    return res / systemMass;
}

/*
Prepars a folder for the addition of files. If folder dpes not exist, it will create one, and if the
folder already does exists, it will delete the contents. 
*/
void Solver::prepareFolder()
{
    std::string path = dir_path + folder;

    if (!fs::exists(path))
    {
        fs::create_directory(path);
    }
    if (!fs::is_empty(path))
    {
        fs::path dir = path;
        for (auto &p : fs::directory_iterator(dir))
        {
            fs::remove_all(p);
        }
    }
}

/*
Sets the folder to store results
*/
void Solver::setResultsFolder(std::string folder)
{
    this->folder = "results/" + folder;
}

/*
Sets the file to read data from
*/
void Solver::setReadFile(const std::string &readFile)
{
    this->readFile = "/" + readFile;
}

/*
Read data from readFile and insert relevant values into each planet object., whilst also increasing the system mass
*/
void Solver::readData(const std::string &readFile)
{
    std::ifstream inData;
    inData.open(dir_path + readFile);
    double mass, x, y, z, vx, vy, vz;
    systemMass = 0;
    for (auto &it : planets)
    {
        inData >> mass;
        inData >> x >> y >> z;
        inData >> vx >> vy >> vz;
        it.setMass(mass);
        systemMass += it.getMass();
        it.setInitialPosition(x, y, z);
        it.setInitialVelocity(it.YEARS * vx, it.YEARS * vy, it.YEARS * vz);
    }
}


/*
Return a planet object
*/
Planet Solver::getPlanet(int index)
{
    return planets[index];
}

