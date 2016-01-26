#ifndef INITIALVALUEGEODESIC_H
#define INITIALVALUEGEODESIC_H

#include "Mesh.cuh"

class InitialValueGeodesic
{
public:

	struct GeodesicKeyPoint
	{
		unsigned edgeIndex;
		double pos;
		bool isInterior;
		unsigned faceIndex;
		Vector3D facePos3D;

		__host__ __device__ GeodesicKeyPoint() { isInterior = false; edgeIndex = -1; faceIndex = -1; }
	};

	__device__ InitialValueGeodesic();
	__device__ ~InitialValueGeodesic();

	__device__ void AssignMesh(Mesh *mesh_);

	__device__ void AssignStartPoint(unsigned startPointIndex_);
	__device__ void AssignStartPoint(unsigned startPointFaceIndex_, Vector3D startPointPos_);

	__device__ void AssignStartDirection(Vector3D direction_);
	__device__ void AssignFirstKeyPoint(unsigned iEdge, double pos);

	__device__ void AssignLength(double geodesicLength_);
	
	__device__ GeodesicKeyPoint BuildGeodesicPath();

	/*__device__ void OutputGeodesicPathWithMesh(const char* fileName);*/

private:
// 	__device__ GeodesicKeyPoint ProjectDirectionOnto1RingNeighOfSource();
// 	__device__ GeodesicKeyPoint ProjectDirectionOntoFace();

private:
	Mesh *mesh;

	unsigned startPointIndex;
	unsigned startPointFaceIndex; Vector3D startPointPos;

	Vector3D direction;
	double geodesicLength;
	GeodesicKeyPoint firstKeyPoint;

};

#endif