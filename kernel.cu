#include <iostream>
#include "SVG.cuh"

using namespace std;

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		cout << "[.exe] [in.obj]" << endl;
		return -1;
	}

	Mesh mesh, d_mesh;
	mesh.LoadFromFile(argv[1]);
	mesh.copyToGPU(&d_mesh);


	SVG svg;

	mesh.clear();
	d_mesh.clearGPU();
    return 0;
}