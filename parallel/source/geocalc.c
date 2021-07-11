/***************************************************************************
 *   Copyright (C) 2004-2008 by OpenFVM team                               *
 *   http://sourceforge.net/projects/openfvm/                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <math.h>

#include "mesh.h"
#include "globals.h"
#include "geocalc.h"

double
GeoMagVector (msh_vector n1)
{

  double mag;

  mag = sqrt (n1.x * n1.x + n1.y * n1.y + n1.z * n1.z);

  return mag;

}

msh_vector
GeoAddVectorVector (msh_vector n2, msh_vector n1)
{

  msh_vector rv;

  rv.x = n2.x + n1.x;
  rv.y = n2.y + n1.y;
  rv.z = n2.z + n1.z;

  return rv;

}

msh_vector
GeoSubVectorVector (msh_vector n2, msh_vector n1)
{

  msh_vector rv;

  rv.x = n2.x - n1.x;
  rv.y = n2.y - n1.y;
  rv.z = n2.z - n1.z;

  return rv;

}

msh_vector
GeoMultVectorVector (msh_vector n2, msh_vector n1)
{

  msh_vector rv;

  rv.x = n2.x * n1.x;
  rv.y = n2.y * n1.y;
  rv.z = n2.z * n1.z;

  return rv;

}

msh_vector
GeoDivVectorVector (msh_vector n2, msh_vector n1)
{

  msh_vector rv;

  rv.x = n2.x / n1.x;
  rv.y = n2.y / n1.y;
  rv.z = n2.z / n1.z;

  return rv;

}

msh_vector
GeoMultScalarVector (double s, msh_vector n1)
{

  msh_vector rv;

  rv.x = s * n1.x;
  rv.y = s * n1.y;
  rv.z = s * n1.z;

  return rv;

}

double
GeoDotVectorVector (msh_vector n2, msh_vector n1)
{

  double dot;

  dot = n2.x * n1.x + n2.y * n1.y + n2.z * n1.z;

  return dot;

}

msh_vector
GeoCrossVector (msh_vector v1, msh_vector v2)
{

  msh_vector rv;

  rv.x = v1.y * v2.z - v1.z * v2.y;
  rv.y = v1.z * v2.x - v1.x * v2.z;
  rv.z = v1.x * v2.y - v1.y * v2.x;

  return rv;

}

msh_vector
GeoNormalizeVector (msh_vector v1)
{

  msh_vector rv;

  double length;
  double factor;
  double min_normal_length;

  rv = v1;

  length = GeoMagVector (v1);

  min_normal_length = 0.000000000001f;

  if (length < min_normal_length)
    {
      rv.x = 1.0;
      rv.y = 0.0;
      rv.z = 0.0;

      return rv;
    }

  factor = 1.0 / length;

  rv.x *= factor;
  rv.y *= factor;
  rv.z *= factor;

  return rv;

}

double
GeoCalcAngle (msh_vector v1, msh_vector v2)
{

  double angle;
  double dot;

  v1 = GeoNormalizeVector (v1);
  v2 = GeoNormalizeVector (v2);

  dot = GeoDotVectorVector (v1, v2) / (GeoMagVector (v1) * GeoMagVector (v2));

  if (dot >= 1.0)
    return 0.0;

  if (dot <= -1.0)
    return PI;

  angle = acos (dot);

  return angle;
}

double
GeoDotProduct (msh_vector v1, msh_vector v2)
{

  double l1, l2;
  double ang;

  l1 = GeoMagVector (v1);
  l2 = GeoMagVector (v2);

  ang = GeoCalcAngle (v1, v2);

  return (l1 * l2 * cos (ang));

}

msh_vector
GeoCalcNormal (msh_vector n1, msh_vector n2, msh_vector n3)
{

  msh_vector rv;

  msh_vector v1;
  msh_vector v2;

  v1 = GeoSubVectorVector (n2, n1);
  v2 = GeoSubVectorVector (n3, n1);

  rv = GeoCrossVector (v1, v2);

  rv = GeoNormalizeVector (rv);

  return rv;

}

msh_vector
GeoCalcCentroid2 (msh_vector n1, msh_vector n2)
{

  msh_vector rv;

  rv.x = (n1.x + n2.x) / 2.0f;
  rv.y = (n1.y + n2.y) / 2.0f;
  rv.z = (n1.z + n2.z) / 2.0f;

  return rv;

}

msh_vector
GeoCalcCentroid3 (msh_vector n1, msh_vector n2, msh_vector n3)
{

  msh_vector rv;

  rv.x = (n1.x + n2.x + n3.x) / 3.0f;
  rv.y = (n1.y + n2.y + n3.y) / 3.0f;
  rv.z = (n1.z + n2.z + n3.z) / 3.0f;

  return rv;

}

msh_vector
GeoCalcCentroid4 (msh_vector n1, msh_vector n2, msh_vector n3, msh_vector n4)
{

  msh_vector rv;

  rv.x = (n1.x + n2.x + n3.x + n4.x) / 4.0f;
  rv.y = (n1.y + n2.y + n3.y + n4.y) / 4.0f;
  rv.z = (n1.z + n2.z + n3.z + n4.z) / 4.0f;

  return rv;

}

msh_vector
GeoCalcCentroid6 (msh_vector n1, msh_vector n2, msh_vector n3, msh_vector n4,
		  msh_vector n5, msh_vector n6)
{

  msh_vector rv;

  rv.x = (n1.x + n2.x + n3.x + n4.x + n5.x + n6.x) / 6.0f;
  rv.y = (n1.y + n2.y + n3.y + n4.y + n5.y + n6.y) / 6.0f;
  rv.z = (n1.z + n2.z + n3.z + n4.z + n5.z + n6.z) / 6.0f;

  return rv;

}

msh_vector
GeoCalcCentroid8 (msh_vector n1, msh_vector n2, msh_vector n3, msh_vector n4,
		  msh_vector n5, msh_vector n6, msh_vector n7, msh_vector n8)
{

  msh_vector rv;

  rv.x = (n1.x + n2.x + n3.x + n4.x + n5.x + n6.x + n7.x + n8.x) / 8.0f;
  rv.y = (n1.y + n2.y + n3.y + n4.y + n5.y + n6.y + n7.y + n8.y) / 8.0f;
  rv.z = (n1.z + n2.z + n3.z + n4.z + n5.z + n6.z + n7.z + n8.z) / 8.0f;

  return rv;

}

double
GeoCalcTriArea (msh_vector n1, msh_vector n2, msh_vector n3)
{

  double area;

  msh_vector c[3];
  msh_vector normal;
  msh_vector sum;

  c[0] = GeoCrossVector (n1, n2);
  c[1] = GeoCrossVector (n2, n3);
  c[2] = GeoCrossVector (n3, n1);

  sum.x = c[0].x + c[1].x + c[2].x;
  sum.y = c[0].y + c[1].y + c[2].y;
  sum.z = c[0].z + c[1].z + c[2].z;

  normal = GeoCalcNormal (n1, n2, n3);

  area =
    (0.5 * LABS (normal.x * sum.x + normal.y * sum.y + normal.z * sum.z));

  return area;
}

double
GeoCalcQuadArea (msh_vector n1, msh_vector n2, msh_vector n3, msh_vector n4)
{

  double area;

  area = 0.0;

  area += GeoCalcTriArea (n1, n2, n3);
  area += GeoCalcTriArea (n1, n3, n4);

  return area;

}

double
GeoCalcTetraVolume (msh_vector n1, msh_vector n2, msh_vector n3,
		    msh_vector n4)
{

  double volume;

  volume =
    LABS ((n2.x - n1.x) * ((n3.y - n1.y) * (n4.z - n1.z) -
			   (n4.y - n1.y) * (n3.z - n1.z)) - (n3.x -
							     n1.x) * ((n2.y -
								       n1.y) *
								      (n4.z -
								       n1.z) -
								      (n4.y -
								       n1.y) *
								      (n2.z -
								       n1.
								       z)) +
	  (n4.x - n1.x) * ((n2.y - n1.y) * (n3.z - n1.z) -
			   (n3.y - n1.y) * (n2.z - n1.z))) / 6.0f;

  return volume;

}

double
GeoCalcHexaVolume (msh_vector n1, msh_vector n2, msh_vector n3, msh_vector n4,
		   msh_vector n5, msh_vector n6, msh_vector n7, msh_vector n8)
{

  double volume;

  msh_vector c;

  volume = 0.0;

  c = GeoCalcCentroid8 (n1, n2, n3, n4, n5, n6, n7, n8);

  volume += GeoCalcTetraVolume (c, n1, n2, n3);
  volume += GeoCalcTetraVolume (c, n3, n4, n1);
  volume += GeoCalcTetraVolume (c, n8, n7, n6);
  volume += GeoCalcTetraVolume (c, n5, n6, n8);
  volume += GeoCalcTetraVolume (c, n6, n7, n3);
  volume += GeoCalcTetraVolume (c, n6, n2, n3);
  volume += GeoCalcTetraVolume (c, n4, n1, n5);
  volume += GeoCalcTetraVolume (c, n4, n5, n8);
  volume += GeoCalcTetraVolume (c, n1, n5, n6);
  volume += GeoCalcTetraVolume (c, n1, n6, n2);
  volume += GeoCalcTetraVolume (c, n8, n7, n3);
  volume += GeoCalcTetraVolume (c, n8, n3, n4);

  return volume;

}

double
GeoCalcPrismVolume (msh_vector n1, msh_vector n2, msh_vector n3,
		    msh_vector n4, msh_vector n5, msh_vector n6)
{

  double volume;

  volume = 0.0;

  volume += GeoCalcTetraVolume (n1, n2, n3, n6);
  volume += GeoCalcTetraVolume (n2, n4, n5, n6);
  volume += GeoCalcTetraVolume (n1, n2, n4, n6);
  volume += GeoCalcTetraVolume (n1, n2, n3, n4);
  volume += GeoCalcTetraVolume (n2, n4, n5, n3);
  volume += GeoCalcTetraVolume (n4, n5, n6, n3);
  volume *= 0.5;

  return volume;

}

void
GeoCalcRotation (msh_vector * v, msh_vector axis, msh_vector v0, double angle)
{

  double ax, bx, cx;
  double ay, by, cy;
  double az, bz, cz;

  ax = 1 + (1 - cos (angle)) * (axis.x * axis.x - 1);
  bx = -axis.z * sin (angle) + (1 - cos (angle)) * axis.x * axis.y;
  cx = axis.y * sin (angle) + (1 - cos (angle)) * axis.x * axis.z;

  ay = axis.z * sin (angle) + (1 - cos (angle)) * axis.x * axis.y;
  by = 1 + (1 - cos (angle)) * (axis.y * axis.y - 1);
  cy = -axis.x * sin (angle) + (1 - cos (angle)) * axis.y * axis.z;

  az = -axis.y * sin (angle) + (1 - cos (angle)) * axis.x * axis.z;
  bz = axis.x * sin (angle) + (1 - cos (angle)) * axis.y * axis.z;
  cz = 1 + (1 - cos (angle)) * (axis.z * axis.z - 1);

  v->x = (ax * v0.x + bx * v0.y + cx * v0.z);
  v->y = (ay * v0.x + by * v0.y + cy * v0.z);
  v->z = (az * v0.x + bz * v0.y + cz * v0.z);


}

char
GeoCalcSegSegIntersection (msh_vector a, msh_vector b, msh_vector c,
			   msh_vector d, msh_vector * p)
{

  double s, t;			/* The two parameters of the parametric eqns. */
  double num, denom;		/* Numerator and denoninator of equations. */

  char code = '?';		/* Return char characterizing intersection. */

  double C_EPS = 1E-8;

  denom = a.x * (d.y - c.y) +
    b.x * (c.y - d.y) + d.x * (b.y - a.y) + c.x * (a.y - b.y);

  /* If denom is zero, then segments are parallel: handle separately. */
  if (LABS (denom) < C_EPS)
    {
      code = '0';
      return code;
    }

  num = a.x * (d.y - c.y) + c.x * (a.y - d.y) + d.x * (c.y - a.y);

  if (LABS (num) < C_EPS)
    {
      if ((num > denom) && (num < denom))
	{
	  code = 'v';
	  return code;
	}
    }

  s = num / denom;

  num = -(a.x * (c.y - b.y) + b.x * (a.y - c.y) + c.x * (b.y - a.y));

  if (LABS (num) < C_EPS)
    {
      if ((num > denom) && (num < denom))
	{
	  code = 'v';
	  return code;
	}
    }

  t = num / denom;

  if ((s > C_EPS) && (s < 1 - C_EPS) && (t > C_EPS) && (t < 1 - C_EPS))
    code = '1';
  else if ((s < 0) || (s > 1) || (t < 0) || (t > 1))
    code = '0';

  p->x = (a.x + s * (b.x - a.x));
  p->y = (a.y + s * (b.y - a.y));
  p->z = 0.0;

  return code;

}

/*
char GeoCalcSegFaceIntersection(msh_face F, msh_vector q0, msh_vector r0)
{
 
	int code = '?';

	stl_vertex q;
	stl_vertex r;

   	stl_facet F1, F2;

	q.x = q0.x;
	q.y = q0.y;
	q.z = q0.z;

	r.x = r0.x;
	r.y = r0.y;
	r.z = r0.z;

	if (F.type == TRIANGLE)
	{
		F1.vertex[0].x = nodes[F.node[0]].x;
		F1.vertex[0].y = nodes[F.node[0]].y;
		F1.vertex[0].z = nodes[F.node[0]].z;

		F1.vertex[1].x = nodes[F.node[1]].x;
		F1.vertex[1].y = nodes[F.node[1]].y;
		F1.vertex[1].z = nodes[F.node[1]].z;

		F1.vertex[2].x = nodes[F.node[2]].x;
		F1.vertex[2].y = nodes[F.node[2]].y;
		F1.vertex[2].z = nodes[F.node[2]].z;

		code = StlSegTriInt(F1, q, r);

	}

	if (F.type == QUADRANGLE)
	{
		F1.vertex[0].x = nodes[F.node[0]].x;
		F1.vertex[0].y = nodes[F.node[0]].y;
		F1.vertex[0].z = nodes[F.node[0]].z;

		F1.vertex[1].x = nodes[F.node[1]].x;
		F1.vertex[1].y = nodes[F.node[1]].y;
		F1.vertex[1].z = nodes[F.node[1]].z;

		F1.vertex[2].x = nodes[F.node[2]].x;
		F1.vertex[2].y = nodes[F.node[2]].y;
		F1.vertex[2].z = nodes[F.node[2]].z;

		code = StlSegTriInt(F1, q, r);

		if (code != '1')
		{

			F2.vertex[0].x = nodes[F.node[0]].x;
			F2.vertex[0].y = nodes[F.node[0]].y;
			F2.vertex[0].z = nodes[F.node[0]].z;

			F2.vertex[1].x = nodes[F.node[2]].x;
			F2.vertex[1].y = nodes[F.node[2]].y;
			F2.vertex[1].z = nodes[F.node[2]].z;

			F2.vertex[2].x = nodes[F.node[3]].x;
			F2.vertex[2].y = nodes[F.node[3]].y;
			F2.vertex[2].z = nodes[F.node[3]].z;

			code = StlSegTriInt(F2, q, r);

		}

	}
    
    return code;

}

void GeoTransformGlobalLocalCoordinates(msh_vector *vl, msh_vector vg, msh_vector il, msh_vector jl, msh_vector kl)
{

	msh_vector ig, jg, kg;

	double t[3][3];

	ig.x = +1.0; ig.y = +0.0; ig.z = +0.0; 
	jg.x = +0.0; jg.y = +1.0; jg.z = +0.0; 
	kg.x = +0.0; kg.y = +0.0; kg.z = +1.0; 

	t[0][0] = GeoDotVectorVector(ig, il);
	t[1][0] = GeoDotVectorVector(ig, jl);
	t[2][0] = GeoDotVectorVector(ig, kl);

	t[0][1] = GeoDotVectorVector(jg, il);
	t[1][1] = GeoDotVectorVector(jg, jl);
	t[2][1] = GeoDotVectorVector(jg, kl);

	t[0][2] = GeoDotVectorVector(kg, il);
	t[1][2] = GeoDotVectorVector(kg, jl);
	t[2][2] = GeoDotVectorVector(kg, kl);

	vl->x = (vg.x * t[0][0] + vg.y * t[0][1] + vg.z * t[0][2]);
	vl->y = (vg.x * t[1][0] + vg.y * t[1][1] + vg.z * t[1][2]);
	vl->z = (vg.x * t[2][0] + vg.y * t[2][1] + vg.z * t[2][2]);

}


void GeoTransformLocalGlobalCoordinates(msh_vector *vg, msh_vector vl, msh_vector il, msh_vector jl, msh_vector kl)
{

	msh_vector ig, jg, kg;

	double t[3][3];

	ig.x = +1.0; ig.y = +0.0; ig.z = +0.0; 
	jg.x = +0.0; jg.y = +1.0; jg.z = +0.0; 
	kg.x = +0.0; kg.y = +0.0; kg.z = +1.0; 

	t[0][0] = GeoDotVectorVector(ig, il);
	t[0][1] = GeoDotVectorVector(ig, jl);
	t[0][2] = GeoDotVectorVector(ig, kl);

	t[1][0] = GeoDotVectorVector(jg, il);
	t[1][1] = GeoDotVectorVector(jg, jl);
	t[1][2] = GeoDotVectorVector(jg, kl);

	t[2][0] = GeoDotVectorVector(kg, il);
	t[2][1] = GeoDotVectorVector(kg, jl);
	t[2][2] = GeoDotVectorVector(kg, kl);

	vg->x = (vl.x * t[0][0] + vl.y * t[0][1] + vl.z * t[0][2]);
	vg->y = (vl.x * t[1][0] + vl.y * t[1][1] + vl.z * t[1][2]);
	vg->z = (vl.x * t[2][0] + vl.y * t[2][1] + vl.z * t[2][2]);

}

*/
