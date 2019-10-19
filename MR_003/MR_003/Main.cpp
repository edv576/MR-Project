#include <string>
#include <filesystem>
#include <iostream>
#include <Eigen/Dense>
#include <string>
#include "OFF_PLYConverter.h"


namespace fs = std::filesystem;

using namespace Eigen;
using namespace std;

int main() {

	printf("Please decide what to do:\n");
	printf("1 to process models\n");
	printf("2 to analize model\n");
	printf("Please proceed with your selection ... ");
	

	int option;

	scanf("%d", &option);

	if (option == 1) {

		std::string path = "DATA3";
		for (const auto& entry : fs::directory_iterator(path)) {
			std::cout << entry.path() << std::endl;
			string s = entry.path().string();
			string s2 = s.substr(s.find("\\") + 1);
			string sName = s2.substr(0, s2.find("."));
			string s3 = "DATA3/" + s2;
			char const* c = s3.data();
			FILE* fo = fopen(c, "r");

			string s3_resized = "DATA4/" + sName + ".ply";

			char const* c2 = s3_resized.data();

			FILE* fd = fopen(c2, "w");

			OFF_PLYConverter* converter = new OFF_PLYConverter();

			converter->Process_Post_Norm_PLY(fo, fd);

			fclose(fo);
			fclose(fd);
		}
	}
	else
	{
		string s = "DATA4/3.ply";

		char const* c = s.data();
		FILE* f = fopen(c, "r");

		OFF_PLYConverter* converter = new OFF_PLYConverter();

		converter->Read_PLY_File(f);


	}

	getchar();

		


	
}