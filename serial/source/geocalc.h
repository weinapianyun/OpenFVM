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

double GeoMagVector (msh_vector n1);
msh_vector GeoAddVectorVector (msh_vector n2, msh_vector n1);
msh_vector GeoSubVectorVector (msh_vector n2, msh_vector n1);
msh_vector GeoMultVectorVector (msh_vector n2, msh_vector n1);
msh_vector GeoDivVectorVector (msh_vector n2, msh_vector n1);
msh_vector GeoMultScalarVector (double s, msh_vector n1);
double GeoDotVectorVector (msh_vector n2, msh_vector n1);
msh_vector GeoMakeVector (msh_vector n2, msh_vector n1);
msh_vector GeoNormalizeVector (msh_vector n1);
msh_vector GeoCrossVector (msh_vector v1, msh_vector v2);
double GeoCalcAngle (msh_vector v1, msh_vector v2);
msh_vector GeoCalcNormal (msh_vector n1, msh_vector n2, msh_vector n3);

double GeoCalcTriArea (msh_vector n1, msh_vector n2, msh_vector n3);
double GeoCalcQuadArea (msh_vector n1, msh_vector n2, msh_vector n3,
			msh_vector n4);

double GeoCalcTetraVolume (msh_vector n1, msh_vector n2, msh_vector n3,
			   msh_vector n4);
double GeoCalcHexaVolume (msh_vector n1, msh_vector n2, msh_vector n3,
			  msh_vector n4, msh_vector n5, msh_vector n6,
			  msh_vector n7, msh_vector n8);
double GeoCalcPrismVolume (msh_vector n1, msh_vector n2, msh_vector n3,
			   msh_vector n4, msh_vector n5, msh_vector n6);

msh_vector GeoCalcCentroid2 (msh_vector n1, msh_vector n2);
msh_vector GeoCalcCentroid3 (msh_vector n1, msh_vector n2, msh_vector n3);
msh_vector GeoCalcCentroid4 (msh_vector n1, msh_vector n2, msh_vector n3,
			     msh_vector n4);
msh_vector GeoCalcCentroid6 (msh_vector n1, msh_vector n2, msh_vector n3,
			     msh_vector n4, msh_vector n5, msh_vector n6);
msh_vector GeoCalcCentroid8 (msh_vector n1, msh_vector n2, msh_vector n3,
			     msh_vector n4, msh_vector n5, msh_vector n6,
			     msh_vector n7, msh_vector n8);

void GeoCalcShapeFunctionsTriangle (msh_vector n1, msh_vector n2,
				    msh_vector n3, double *dx, double *dy,
				    double *c1x);
void GeoCalcShapeFunctionsTetrahedron (msh_vector n1, msh_vector n2,
				       msh_vector n3, msh_vector n4,
				       double *dx, double *dy, double *dz,
				       double *c1x, double *c2x, double *c2y);
char GeoCalcSegSegIntersection (msh_vector a, msh_vector b, msh_vector c,
				msh_vector d, msh_vector * p);
char GeoCalcSegPlaneIntersection (msh_plane P, msh_vector q, msh_vector r,
				  msh_vector * p);
void GeoTransformLocalGlobalCoordinates (msh_vector * vg, msh_vector vl,
					 msh_vector il, msh_vector jl,
					 msh_vector kl);
void GeoTransformGlobalLocalCoordinates (msh_vector * vl, msh_vector vg,
					 msh_vector il, msh_vector jl,
					 msh_vector kl);
