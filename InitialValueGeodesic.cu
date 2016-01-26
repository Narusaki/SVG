#include "InitialValueGeodesic.cuh"

__device__ InitialValueGeodesic::InitialValueGeodesic()
{
	mesh = NULL;
	startPointIndex = -1;
	startPointFaceIndex = -1;
	firstKeyPoint.edgeIndex = -1;
	firstKeyPoint.isInterior = false;
}

__device__ InitialValueGeodesic::~InitialValueGeodesic()
{

}

__device__ void InitialValueGeodesic::AssignMesh(Mesh *mesh_)
{
	mesh = mesh_;
}

__device__ void InitialValueGeodesic::AssignStartPoint(unsigned startPointIndex_)
{
	startPointIndex = startPointIndex_;
	startPointFaceIndex = -1;
}

__device__ void InitialValueGeodesic::AssignStartPoint(unsigned startPointFaceIndex_, Vector3D startPointPos_)
{
	startPointFaceIndex = startPointFaceIndex_;
	startPointPos = startPointPos_;
	startPointIndex = -1;
}

__device__ void InitialValueGeodesic::AssignStartDirection(Vector3D direction_)
{
	direction = direction_;
	direction.normalize();
}

__device__ void InitialValueGeodesic::AssignFirstKeyPoint(unsigned iEdge, double pos)
{
	firstKeyPoint.edgeIndex = iEdge;
	firstKeyPoint.pos = pos;
	firstKeyPoint.isInterior = false;
}


__device__ void InitialValueGeodesic::AssignLength(double geodesicLength_)
{
	geodesicLength = geodesicLength_;
}

__device__ InitialValueGeodesic::GeodesicKeyPoint InitialValueGeodesic::BuildGeodesicPath()
{
	// collect adjacent faces
	GeodesicKeyPoint keyPoint;
	Vector2D src2D;
	
	/*if (firstKeyPoint.edgeIndex != -1)*/
		keyPoint = firstKeyPoint;
// 	else if (startPointIndex != -1)
// 		keyPoint = ProjectDirectionOnto1RingNeighOfSource();
// 	else
// 		keyPoint = ProjectDirectionOntoFace();
	Vector2D keyPoint2D = Vector2D(keyPoint.pos, 0.0);

	// make sure the 2d position of the source
	double l0 = mesh->edges[keyPoint.edgeIndex].edgeLen;
	double l1 = mesh->edges[mesh->edges[keyPoint.edgeIndex].nextEdge].edgeLen;
	double l2 = mesh->edges[mesh->edges[mesh->edges[keyPoint.edgeIndex].nextEdge].nextEdge].edgeLen;
	src2D.x = (l2*l2 + l0*l0 - l1*l1) / (2.0 * l0);
	src2D.y = sqrt(fabs(l2*l2 - src2D.x*src2D.x));

	// refine src2D for the case where source in the interior of a face
	if (startPointIndex == -1)
	{
		double signedAreas[3], lambda[3], totalArea = 0.0;
		unsigned IStart = 0;
		if (mesh->faces[startPointFaceIndex].verts[0] == mesh->edges[keyPoint.edgeIndex].verts[0])
			IStart = 0;
		else if (mesh->faces[startPointFaceIndex].verts[1] == mesh->edges[keyPoint.edgeIndex].verts[0])
			IStart = 1;
		else IStart = 2;

		unsigned i = IStart;
		do
		{
			Vector3D p0 = mesh->verts[mesh->faces[startPointFaceIndex].verts[i]].pos;
			Vector3D p1 = mesh->verts[mesh->faces[startPointFaceIndex].verts[(i+1)%3]].pos;
			Vector3D signedAreaVector = (p0 - startPointPos) ^ (p1 - startPointPos);
			signedAreas[(i - IStart + 2) % 3] = signedAreaVector.length();
			totalArea += signedAreas[(i-IStart+2)%3];
			i = (i + 1) % 3;
		} while (i != IStart);
		for (unsigned i = 0; i < 3; ++i) lambda[i] = signedAreas[i] / totalArea;

		src2D = lambda[1] * Vector2D(l0, 0.0) + lambda[2] * src2D;
	}
	// start length
	double curLen = (src2D - keyPoint2D).length();
	GeodesicKeyPoint prevKeyPoint = keyPoint;
	int keyPointsNum = 1;

	while (curLen < geodesicLength)
	{
		unsigned e0 = mesh->edges[keyPoint.edgeIndex].twinEdge;
		// if e0 == -1, reach the boundary
		if (e0 == -1) break;
		unsigned e1 = mesh->edges[e0].nextEdge;
		unsigned e2 = mesh->edges[e1].nextEdge;
		double l3 = mesh->edges[e1].edgeLen, l4 = mesh->edges[e2].edgeLen;
		double l0 = mesh->edges[keyPoint.edgeIndex].edgeLen;
		// opposite vertex position
		Vector2D opPos;
		opPos.x = (l3*l3 + l0*l0 - l4*l4) / (2.0*l0);
		opPos.y = -sqrt(fabs(l3*l3 - opPos.x*opPos.x));
		// two endpoints position
		Vector2D endPoint0(0.0, 0.0), endPoint1(l0, 0.0);

		// the next geodesic key point
		GeodesicKeyPoint nextKeyPoint;

		if (keyPoint.pos != 0.0)										// not passing a vertex
		{
			Vector2D projRay = keyPoint2D - src2D, opRay = opPos - src2D;
			projRay.normalize(); opRay.normalize();
			double crossProd = projRay ^ opRay;
			double a[2][2], b[2];
			Vector2D newSrc2D;
			double lambda = keyPoint.pos / mesh->edges[keyPoint.edgeIndex].edgeLen;
			// test if the geodesic path passes the opposite vertex
			if (crossProd > -LDOUBLE_EPS && crossProd < LDOUBLE_EPS)	// pass through the opposite vertex
			{
				nextKeyPoint.edgeIndex = e2;
				nextKeyPoint.pos = 0.0;
				newSrc2D.x = (l0*l0 + l4*l4 - l3*l3) / (2.0*l4);
				newSrc2D.y = sqrt(fabs(l0*l0 - newSrc2D.x*newSrc2D.x));
				newSrc2D *= (1 - lambda); newSrc2D.x = l4 - newSrc2D.x;
			}
			else                                                       // not pass through the opposite vertex
			{
				Vector2D edgeVec, edgeStart;

				if (crossProd > LDOUBLE_EPS)							// toLeft() true, project on e1
				{
					edgeVec = opPos; edgeStart = Vector2D(0.0, 0.0);
					nextKeyPoint.edgeIndex = e1;
					newSrc2D.x = (l0*l0 + l3*l3 - l4*l4) / (2.0*l3);
					newSrc2D.y = sqrt(fabs(l0*l0 - newSrc2D.x*newSrc2D.x));
					newSrc2D *= lambda;
				}
				else														// toLeft() false, project on e2
				{
					edgeVec = endPoint1 - opPos; edgeStart = opPos;
					nextKeyPoint.edgeIndex = e2;
					newSrc2D.x = (l0*l0 + l4*l4 - l3*l3) / (2.0*l4);
					newSrc2D.y = sqrt(fabs(l0*l0 - newSrc2D.x*newSrc2D.x));
					newSrc2D *= (1 - lambda); newSrc2D.x = l4 - newSrc2D.x;
				}
				// solve 2-order linear system
				a[0][0] = projRay.x; a[0][1] = -edgeVec.x;
				a[1][0] = projRay.y; a[1][1] = -edgeVec.y;
				b[0] = edgeStart.x - src2D.x; b[1] = edgeStart.y - src2D.y;

				double det = a[0][0] * a[1][1] - a[0][1] * a[1][0];			// this cannot be zero
				double t0 = (a[1][1] * b[0] - a[0][1] * b[1]) / det;
				double t1 = (-a[1][0] * b[0] + a[0][0] * b[1]) / det;
				nextKeyPoint.pos = t1 * edgeVec.length();
			}
			// update src2D
			src2D = newSrc2D;
		}
		else                                                            // passing a vertex
		{
			// total sweep angle
			double sweepAngle = mesh->verts[mesh->edges[keyPoint.edgeIndex].verts[0]].angle / 2.0;
			// angle that already formed
			double curAngle = src2D.x / src2D.length();
			if (curAngle > 1.0) curAngle = 1.0; else if (curAngle < -1.0) curAngle = -1.0;
			curAngle = acos(curAngle);
			unsigned startEdge = keyPoint.edgeIndex;
			double len0, len1, len2;
			// while not reach the half of the total angle 
			while (curAngle < sweepAngle)
			{
				startEdge = mesh->edges[startEdge].twinEdge;
				// reach the boundary
				if (startEdge == -1) break;
				len0 = mesh->edges[startEdge].edgeLen;
				startEdge = mesh->edges[startEdge].nextEdge;
				len1 = mesh->edges[startEdge].edgeLen;
				len2 = mesh->edges[mesh->edges[startEdge].nextEdge].edgeLen;

				double curTriangleAngle = (len0*len0 + len1*len1 - len2*len2) / (2.0*len0*len1);
				if (curTriangleAngle > 1.0) curTriangleAngle = 1.0; else if (curTriangleAngle < -1.0) curTriangleAngle = -1.0;
				curTriangleAngle = acos(curTriangleAngle);

				if (curAngle + curTriangleAngle >= sweepAngle) break;
				curAngle += curTriangleAngle;
			}
			if (startEdge == -1) break;

			double localAngle = sweepAngle - curAngle; // bust off angle
			double triangleAngle = (len2*len2 + len0*len0 - len1*len1) / (2.0*len2*len0);
			if (triangleAngle > 1.0) triangleAngle = 1.0; else if (triangleAngle < -1.0) triangleAngle = -1.0;
			triangleAngle = acos(triangleAngle);
			double thirdAngle = PI - localAngle - triangleAngle;

			// nextKeyPoint
			nextKeyPoint.pos = sin(localAngle) * len0 / sin(thirdAngle);
			nextKeyPoint.pos = len2 - nextKeyPoint.pos;
			if (fabs(nextKeyPoint.pos) < LDOUBLE_EPS) nextKeyPoint.pos = 0.0;
			nextKeyPoint.edgeIndex = mesh->edges[startEdge].nextEdge;

			// src2D
			src2D.x = (len1*len1 + len2*len2 - len0*len0) / (2.0 * len2);
			src2D.y = sqrt(fabs(len1*len1 - src2D.x*src2D.x));
		}

		prevKeyPoint = keyPoint;
		keyPoint = nextKeyPoint;
		keyPoint2D = Vector2D(keyPoint.pos, 0.0);
		curLen += (src2D - keyPoint2D).length();
		++keyPointsNum;
	}

	if (curLen > geodesicLength)
	{
		Vector3D p0, p1;
		
		Vector3D v0 = mesh->verts[mesh->edges[keyPoint.edgeIndex].verts[0]].pos;
		Vector3D v1 = mesh->verts[mesh->edges[keyPoint.edgeIndex].verts[1]].pos;
		Vector3D univ = v1 - v0; univ.normalize(); p1 = v0 + univ * keyPoint.pos;

		if (keyPointsNum > 1)
		{
			v0 = mesh->verts[mesh->edges[prevKeyPoint.edgeIndex].verts[0]].pos;
			v1 = mesh->verts[mesh->edges[prevKeyPoint.edgeIndex].verts[1]].pos;
			univ = v1 - v0; univ.normalize(); p0 = v0 + univ * prevKeyPoint.pos;
			
		}
		else
		{
			p0 = startPointIndex == -1 ? startPointPos : mesh->verts[startPointIndex].pos;
		}

		univ = p0 - p1; univ.normalize();
		keyPoint.isInterior = true;
		keyPoint.faceIndex = mesh->edges[keyPoint.edgeIndex].faceId;
		keyPoint.facePos3D = p1 + (curLen - geodesicLength) * univ;
	}
	return keyPoint;
}

// __device__ void InitialValueGeodesic::OutputGeodesicPathWithMesh(const char* fileName)
// {
// 	ofstream output(fileName);
// 
// 	for (unsigned i = 0; i < mesh->m_nVertex; ++i)
// 		output << "v " << mesh->verts[i].pos << endl;
// 
// 	if (startPointIndex != -1) output << "v " << mesh->verts[startPointIndex].pos << endl;
// 	else output << "v " << startPointPos << endl;
// 
// 	for (auto iter = geodesicPath.begin(); iter != geodesicPath.end(); ++iter)
// 	{
// 		Vector3D keyPoint;
// 		if (!iter->isInterior)
// 		{
// 			Vector3D p0 = mesh->verts[mesh->edges[iter->edgeIndex].verts[0]].pos;
// 			Vector3D p1 = mesh->verts[mesh->edges[iter->edgeIndex].verts[1]].pos;
// 			Vector3D univ = p1 - p0; univ.normalize();
// 
// 			keyPoint = p0 + iter->pos * univ;
// 		}
// 		else keyPoint = iter->facePos3D;
// 		output << "v " << keyPoint << endl;
// 	}
// 	for (unsigned i = 0; i < mesh->m_nFace; ++i)
// 		output << "f " << mesh->faces[i].verts[0] + 1 << " "
// 		<< mesh->faces[i].verts[1] + 1 << " "
// 		<< mesh->faces[i].verts[2] + 1 << endl;
// 
// 	output << "l " << mesh->m_nVertex + 1 << " " << mesh->m_nVertex + 2 << endl;
// 	for (unsigned i = 0; i < geodesicPath.size()-1; ++i)
// 		output << "l " << i + mesh->m_nVertex + 2 << " " << i + mesh->m_nVertex + 3 << endl;
// 
// 	output.close();
// }

// __device__ InitialValueGeodesic::GeodesicKeyPoint InitialValueGeodesic::ProjectDirectionOnto1RingNeighOfSource()
// {
// 	Vector3D surfaceNormal = mesh->verts[startPointIndex].m_vNormal;
// 	Vector3D planeNormal = surfaceNormal ^ direction; planeNormal.normalize();
// 	Vector3D planePoint = mesh->verts[startPointIndex].pos;
// 
// 	double minAngle = 1e30;
// 	GeodesicKeyPoint keyPoint;
// 	unsigned curEdge = mesh->verts[startPointIndex].firstEdge;
// 	do
// 	{
// 		unsigned curFace = mesh->edges[curEdge].faceId;
// 		unsigned v0 = -1, v1 = -1, opEdge = -1;
// 		for (unsigned j = 0; j < 3; ++j)
// 			if (mesh->faces[curFace].verts[j] == startPointIndex) { v0 = j; break; }
// 		v0 = (v0 + 1) % 3; v1 = (v0 + 1) % 3;
// 		opEdge = mesh->faces[curFace].edges[v0];
// 		v0 = mesh->faces[curFace].verts[v0]; 
// 		v1 = mesh->faces[curFace].verts[v1];
// 
// 		Vector3D p0 = mesh->verts[v0].pos, p1 = mesh->verts[v1].pos;
// 		double dist0 = (planePoint - p0) * planeNormal, dist1 = (planePoint - p1) * planeNormal;
// 		if (dist0 * dist1 >= 0.0) continue;
// 
// 		double lambda = fabs(dist0) / (fabs(dist0) + fabs(dist1));
// 		Vector3D intersectPoint = (1 - lambda)*p0 + lambda*p1;
// 		Vector3D dstDirection = intersectPoint - mesh->verts[startPointIndex].pos;
// 		dstDirection.normalize();
// 
// 		double angle = direction * dstDirection;
// 		if (angle > 1.0) angle = 1.0; else if (angle < -1.0) angle = -1.0;
// 		angle = acos(angle);
// 
// 		Vector3D rotateAxis = direction ^ dstDirection; rotateAxis.normalize(); // this is the current rotate axis
// 		if (rotateAxis * planeNormal < 0.0) angle = 2.0 * PI - angle;
// 		if (angle > minAngle) continue;
// 		minAngle = angle;
// 		keyPoint.edgeIndex = opEdge;
// 		keyPoint.pos = lambda * mesh->edges[opEdge].edgeLen;
// 
// 		curEdge = mesh->edges[curEdge].twinEdge;
// 		if (curEdge != -1) curEdge = mesh->edges[curEdge].nextEdge;
// 	} while (curEdge != -1 && curEdge != mesh->verts[startPointIndex].firstEdge);
// 	return keyPoint;
// }
// 
// __device__ InitialValueGeodesic::GeodesicKeyPoint InitialValueGeodesic::ProjectDirectionOntoFace()
// {
// 	Vector3D surfaceNormal = mesh->faces[startPointFaceIndex].m_vNormal;
// 	Vector3D planeNormal = surfaceNormal ^ direction; planeNormal.normalize();
// 	Vector3D planePoint = startPointPos;
// 
// 	double minAngle = 1e30;
// 	GeodesicKeyPoint keyPoint;
// 	for (unsigned i = 0; i < 3; ++i)
// 	{
// 		Vector3D p0 = mesh->verts[mesh->faces[startPointFaceIndex].verts[i]].pos;
// 		Vector3D p1 = mesh->verts[mesh->faces[startPointFaceIndex].verts[(i+1)%3]].pos;
// 		double dist0 = (planePoint - p0) * planeNormal, dist1 = (planePoint - p1) * planeNormal;
// 		if (dist0 * dist1 >= 0.0) continue;
// 		double lambda = fabs(dist0) / (fabs(dist0) + fabs(dist1));
// 		Vector3D intersectPoint = (1 - lambda)*p0 + lambda*p1;
// 		Vector3D dstDirection = intersectPoint - startPointPos; dstDirection.normalize();
// 
// 		double angle = direction * dstDirection;
// 		if (angle > 1.0) angle = 1.0; else if (angle < -1.0) angle = -1.0;
// 		angle = acos(angle);
// 
// 		if (angle > minAngle) continue;
// 		minAngle = angle;
// 		keyPoint.edgeIndex = mesh->faces[startPointFaceIndex].edges[i];
// 		keyPoint.pos = lambda * mesh->edges[keyPoint.edgeIndex].edgeLen;
// 	}
// 	return keyPoint;
// }