#pragma once
#include <fstream>
#include <string>
#include <armadillo>

/*Class for writing to files, creating folders and possibly reading from files(future)*/

/*
* This part checks during pre proccessing what system the user is using. If the user is using windows, they should have
* direct.h and io.h located on their system, and if the user has a unix based system, they should have unistd.h and sys/stat.h
* on their system. We only want to include libraries that the user already has, and ignore the others, such that we do not get any compile errors.
*/

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#include <io.h>
#define access _access_s
#define GetCurrentDir _getcwd
#define createFolder _mkdir
#define F_OK 0
std::string slash = "//";
#else
#include <unistd.h>
#include <sys/stat.h>
#define createFolder mkdir
#define GetCurrentDir getcwd
std::string slash = "/";
#endif

/*Creates a folder in current path*/
void create_directory(std::string filename)
{
	char *dir = const_cast<char *>(filename.c_str());
	//mkdir and _mkdir take different arguments, so we need to make sure we are passing them correctly
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	createFolder(dir);
#else
	createFolder(dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
}

/*Returns path to current directory*/
std::string get_current_dir()
{
	char buff[FILENAME_MAX];
	GetCurrentDir(buff, FILENAME_MAX);
	std::string current_working_dir(buff);

	for (int i = 0; i < int(current_working_dir.length()); i++)
	{
		if (current_working_dir[i] == '\\')
		{
			current_working_dir.replace(i, 1, slash);
		}
	}
	return current_working_dir;
}

class FileHandler
{
private:
	std::ofstream outFile;
	std::ifstream inFile;
	std::string folder;
	std::string currentPath;
	bool append = false;

public:
	FileHandler(std::string folder)
	{
		create_directory("results");
		create_directory("results" + slash + folder);
		currentPath = get_current_dir();
		currentPath += slash + "results" + slash + folder + slash;
	}

	void writeToFile(std::string &filename, const arma::vec &v)
	{
		filename = currentPath + filename + ".csv";
		bool does_exist = exists(filename);
		fileOpener(filename, does_exist);

		if (does_exist && outFile.is_open())
		{
			outFile << "----------------- \n";
			outFile << 0 << "\n";
			unsigned int n = v.size();
			for (unsigned int i = 0; i < n; i++)
			{
				outFile << v(i) << "\n";
			}
		}
		else if (outFile.is_open())
		{
			outFile << 0 << "\n";
			unsigned int n = v.size();
			for (unsigned int i = 0; i < n; i++)
			{
				outFile << v(i) << "\n";
			}
		}
		outFile << 0 << "\n";
		outFile.close();
	}

	/*Dont add 0 to beginning and end. Needs a better implementation*/
	void writeRawToFile(std::string &filename, const arma::vec &v)
	{
		filename = currentPath + filename + ".csv";
		bool does_exist = exists(filename);
		fileOpener(filename, does_exist);

		if (does_exist && outFile.is_open())
		{
			outFile << "----------------- \n";
			unsigned int n = v.size();
			for (unsigned int i = 0; i < n; i++)
			{
				outFile << v(i) << "\n";
			}
		}
		else if (outFile.is_open())
		{
			unsigned int n = v.size();
			for (unsigned int i = 0; i < n; i++)
			{
				outFile << v(i) << "\n";
			}
		}
		outFile.close();
	}

	inline bool exists(const std::string &name)
	{
		std::ifstream f(name.c_str());
		return f.good();
	}

	void fileOpener(std::string filename, const bool file)
	{
		(file) ? outFile.open(filename, std::ofstream::app) : outFile.open(filename, std::ofstream::out);
	}
};