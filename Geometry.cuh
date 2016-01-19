#ifndef GEOMETRY_CUH
#define GEOMETRY_CUH

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "Geometry.cuh"

#include <math.h>
#include <cassert>
#include <iostream>

#define DOUBLE_EPS 1e-12
#define LDOUBLE_EPS 1e-6
#define LHDOUBLE_EPS 1e-5
#define	EQUALZERO(x)	(fabs((x)) < DOUBLE_EPS)
#define	LEQUALZERO(x) (fabs((x)) < LDOUBLE_EPS)
#define PI	3.141592653589793238

/////////////////////////////////////////////////////////////
// Vector2D : 2D vector
/////////////////////////////////////////////////////////////
class Vector2D
{
public:
	double x, y;

public:
	__host__ __device__ Vector2D(){ x = 0;	y = 0; }
	// constructions
	__host__ __device__ Vector2D(double xx, double yy)	{ x = xx;	y = yy; }
	__host__ __device__ Vector2D(const Vector2D& v)	{ x = v.x;	y = v.y; }


	// operator
	__host__ __device__ double	  length() const		{ return sqrt(x*x + y*y); }
	__host__ __device__ double	  length2() const		{ return x*x + y*y; }
	__host__ __device__ double	  normalize()	{ double len = length();	if (!EQUALZERO(len)) { x /= len; y /= len; }	return len; }
	__host__ __device__ double normalizeStrict() { double len = length(); if (len != 0.0) { x /= len; y /= len; } return len; }
	__host__ __device__ Vector2D& operator=(const Vector2D& v);
	__host__ __device__ Vector2D& operator+=(const Vector2D& v);
	__host__ __device__ Vector2D& operator-=(const Vector2D& v);
	__host__ __device__ Vector2D& operator*=(double u);
	__host__ __device__ Vector2D& operator/=(double u);
	__host__ __device__ double& operator[](int idx) {
		assert(idx < 2);
		switch (idx)
		{
		case 0: return x;
		case 1: return y;
		}
	}
	__host__ __device__ bool operator == (const Vector2D& right) const
	{
		return x == right.x && y == right.y;
	}
	//Vector2D& operator^=(const Vector2D& v);

	__host__ __device__ bool	Intersect(Vector2D v1, Vector2D v2, Vector2D v3, Vector2D v4);
	__host__ __device__ bool	Intersect(Vector2D v1, Vector2D v2);

	__host__ __device__ friend Vector2D operator+(const Vector2D& lv, const Vector2D& rv);
	__host__ __device__ friend Vector2D operator-(const Vector2D& lv, const Vector2D& rv);
	__host__ __device__ friend Vector2D operator*(const double u, const Vector2D& rv);
	__host__ __device__ friend Vector2D operator*(const Vector2D& lv, const double u);
	__host__ __device__ friend Vector2D operator/(const Vector2D& lv, const double u);
	__host__ __device__ friend double   operator*(const Vector2D& lv, const Vector2D& rv);
	__host__ __device__ friend double operator^(const Vector2D& lv, const Vector2D& rv);
	__host__ __device__ friend std::ostream& operator<< (std::ostream &output, Vector2D& v);

	__host__ __device__ short	AtWhere(Vector2D v0, Vector2D v1);
	__host__ __device__ bool	AtRight(Vector2D v0, Vector2D v1);
	__host__ __device__ bool	AtLeft(Vector2D v0, Vector2D v1);
	__host__ __device__ bool	OnLine(Vector2D v0, Vector2D v1);
	__host__ __device__ double	GetArea(Vector2D v);


};

/////////////////////////////////////////////////////////////
// Vector3D : 3D vector
/////////////////////////////////////////////////////////////

class Vector3D
{
public:
	double x, y, z;

	// constructions
	__host__ __device__ Vector3D()	{ x = 0;	y = 0;	z = 0; }
	__host__ __device__ Vector3D(double xx, double yy, double zz)	{ x = xx;	y = yy;	z = zz; }
	__host__ __device__ Vector3D(const Vector3D& v)	{ x = v.x;	y = v.y;	z = v.z; }

	// operator
	__host__ __device__ double	  length() const		{ return sqrt(x*x + y*y + z*z); }
	__host__ __device__ double	  length2() const		{ return x*x + y*y + z*z; }
	__host__ __device__ double	  normalize()	{ double len = length();	if (!EQUALZERO(len)) { x /= len; y /= len; z /= len; }	return len; }
	__host__ __device__ double normalizeStrict() { double len = length(); if (len != 0.0) { x /= len; y /= len; z /= len; } return len; }

	__host__ __device__ Vector3D& operator=(const Vector3D& v);
	__host__ __device__ Vector3D& operator+=(const Vector3D& v);
	__host__ __device__ Vector3D& operator-=(const Vector3D& v);
	__host__ __device__ Vector3D& operator*=(double u);
	__host__ __device__ Vector3D& operator/=(double u);
	__host__ __device__ Vector3D& operator^=(const Vector3D& v);
	__host__ __device__ double& operator[](int idx) {
		assert(idx < 3);
		switch (idx)
		{
		case 0: return x;
		case 1: return y;
		case 2: return z;
		}
	}

	__host__ __device__ bool operator < (const Vector3D& right) const
	{
		return x - right.x < -LHDOUBLE_EPS ||
			fabs(x - right.x) <= LHDOUBLE_EPS && y - right.y < -LHDOUBLE_EPS ||
			fabs(x - right.x) <= LHDOUBLE_EPS && fabs(y - right.y) <= LHDOUBLE_EPS && z - right.z < -LHDOUBLE_EPS;
	}

	__host__ __device__ bool operator == (const Vector3D& right) const
	{
		return x == right.x && y == right.y && z == right.z;
	}

	__host__ __device__ friend Vector3D operator+(const Vector3D& lv, const Vector3D& rv);
	__host__ __device__ friend Vector3D operator-(const Vector3D& lv, const Vector3D& rv);
	__host__ __device__ friend Vector3D operator*(const double u, const Vector3D& rv);
	__host__ __device__ friend Vector3D operator*(const Vector3D& lv, const double u);
	__host__ __device__ friend Vector3D operator/(const Vector3D& lv, const double u);
	__host__ __device__ friend double   operator*(const Vector3D& lv, const Vector3D& rv);
	__host__ __device__ friend Vector3D operator^(const Vector3D& lv, const Vector3D& rv);
	__host__ __device__ friend std::ostream& operator<< (std::ostream &output, const Vector3D& v);

};

__host__ __device__ double Area2(Vector2D &a, Vector2D &b, Vector2D &c);

__host__ __device__ double SpcDivision(const double &a, const double &b);

__host__ __device__ Vector3D Rotate(Vector3D point, double angle, Vector3D rAxis);

__host__ __device__ bool toLeft(const Vector2D &p, const Vector2D &p0, const Vector2D &p1);

#endif