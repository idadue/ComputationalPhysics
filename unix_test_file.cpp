#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#if __has_include(<direct.h>)
#include <direct.h>
#define GetCurrentDir _getcwd
#define createFolder _mkdir

#else
#include <unistd.h>
#include <sys/stat.h>
#define createFolder mkdir
#define GetCurrentDir getcwd
#endif

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
std::string slash = "//";
#else
std::string slash = "/";
#endif

void create_directory(std::string filename) {
	char* dir = const_cast<char*>(filename.c_str());    //convert string to char*
	#if defined(__CYGWIN__)
	createFolder(dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	#else
	createFolder(dir);
	#endif
}

std::string get_current_dir() {

	char buff[FILENAME_MAX];
	GetCurrentDir(buff, FILENAME_MAX);
	std::string current_working_dir(buff);
	for (int i = 0; i < int(current_working_dir.length()); i++) {
		if (current_working_dir[i] == '\\') {
			current_working_dir.replace(i,1, slash);
		}
	}
	return current_working_dir;
}

void createFile(std::string filename, int n) {
	create_directory(filename);

	std::ofstream ofile;
	std::string file = filename;

	file.append("_" + std::to_string(n) + ".csv");

	std::string dir = get_current_dir();
	ofile.open(dir + slash + filename + slash + file, std::fstream::out);

	if (ofile.is_open()) {
		std::cout << "File opened successfully." << std::endl;

		ofile << "Testing" << std::endl;
		for (int i = 0; i < 10; i++) {
			ofile << i << std::endl;
		}
	}
	else {
		std::cout << "Could not create file!" << std::endl;
	}

	ofile.close();
}

int main() {
	createFile("data", 10);

	return 0;
}