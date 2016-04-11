#ifndef INITIALVALUEGEODESIC_H
#define INITIALVALUEGEODESIC_H

#include "Mesh.cuh"

class InitialValueGeodesic
{
public:

	struct GeodesicKeyPoint	// 41 bytes
	{
		unsigned edgeIndex;
		double pos;
		bool isInterior;
		unsigned faceIndex;
		Vector3D facePos3D;

		__host__  GeodesicKeyPoint() { isInterior = false; edgeIndex = -1; faceIndex = -1; }
	};

	 InitialValueGeodesic();
	 ~InitialValueGeodesic();

	 void AssignMesh(Mesh *mesh_);

	 void AssignStartPoint(unsigned startPointIndex_);
	 void AssignStartPoint(unsigned startPointFaceIndex_, Vector3D startPointPos_);

	 void AssignStartDirection(Vector3D direction_);
	 void AssignFirstKeyPoint(unsigned iEdge, double pos);

	 void AssignLength(double geodesicLength_);

	 GeodesicKeyPoint BuildGeodesicPath();

	/* void OutputGeodesicPathWithMesh(const char* fileName);*/

private:
	// 	 GeodesicKeyPoint ProjectDirectionOnto1RingNeighOfSource();
	// 	 GeodesicKeyPoint ProjectDirectionOntoFace();

private:
	Mesh *mesh;

	unsigned startPointIndex;
	unsigned startPointFaceIndex; Vector3D startPointPos;

	Vector3D direction;
	double geodesicLength;
	GeodesicKeyPoint firstKeyPoint;

};

#endif