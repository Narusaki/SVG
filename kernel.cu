#include <iostream>
#include "SVG.cuh"
#include "book.cuh"

using namespace std;

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		cout << "[.exe] [in.obj]" << endl;
		return -1;
	}
	HANDLE_ERROR(cudaSetDevice(0));
	Mesh mesh, d_mesh;
	mesh.LoadFromFile(argv[1]);
	cout << "Mesh is loaded." << endl;
	mesh.copyToGPU(&d_mesh);

	cout << "Mesh is copied into GPU." << endl;

	SVG svg;
	svg.AssignMesh(&mesh, &d_mesh);
	svg.SetParameters(500);
	svg.Allocation();

	cout << "Please check initial memory." << endl;
	system("pause");

	svg.ConstructSVG();

	cout << "SVG test ended." << endl;
	cout << "Please check final memory." << endl;
	system("pause");

	svg.Free();
	mesh.clear();
	d_mesh.clearGPU();

	cout << "Please check SVG-structure + Mesh memory." << endl;
	system("pause");

	HANDLE_ERROR(cudaDeviceReset());
    return 0;
}