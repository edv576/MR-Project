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
		string s1 = "DATA4/1.ply";

		char const* c1 = s1.data();
		FILE* f1 = fopen(c1, "r");

		OFF_PLYConverter* converter1 = new OFF_PLYConverter();

		converter1->Read_PLY_File(f1);

		string s2 = "DATA4/2.ply";

		char const* c2 = s2.data();
		FILE* f2 = fopen(c2, "r");

		OFF_PLYConverter* converter2 = new OFF_PLYConverter();

		converter2->Read_PLY_File(f2);

		float d1 = converter1->CalculateDiameter();
		float d2 = converter2->CalculateDiameter();

		float comp1 = converter1->CalculateCompactness();
		float comp2 = converter2->CalculateCompactness();

		VectorXi histogram_Bary_RandVert1 = converter1->CalculateHistogram_Bary_RandVert(1000, 10);
		VectorXi histogram_Bary_RandVert2 = converter2->CalculateHistogram_Bary_RandVert(1000, 10);

		//float comp1 = converter1->CalculateCompactness()

		printf((s1+ " - Diameter: ").data());
		printf(std::to_string(d1).data());

		printf("\n");

		printf((s2 + " - Diameter: ").data());
		printf(std::to_string(d2).data());

		printf("\n");

		printf((s1 + " - Compactness: ").data());
		printf(std::to_string(comp1).data());

		printf("\n");

		printf((s2 + " - Compactness: ").data());
		printf(std::to_string(comp2).data());



	}

	getchar();

		


	
}