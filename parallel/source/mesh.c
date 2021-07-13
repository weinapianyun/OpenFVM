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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// PETSc
#include "petscksp.h"

#include "globals.h"
#include "param.h"
#include "mesh.h"
#include "bcond.h"
#include "geocalc.h"
#include "octree.h"

void
MshFreeMemory ()
{

  int i;

  int patch, face, element;

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      if (faces[face].type == TRIANGLE || faces[face].type == QUADRANGLE)
	{
	  free (faces[face].node);
	}

    }

  for (i = 0; i < nbpatches; i++)
    {
      patch = i;

      free (patches[patch].node);

    }

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type == TETRAHEDRON)
	{

	  free (elements[element].b);
	  free (elements[element].c);
	  free (elements[element].d);

	}

    }

  if (nbnodes > 0)
    {
      nbnodes = 0;
      free (nodes);
    }

  if (nbfaces > 0)
    {
      nbfaces = 0;
      free (faces);
    }

  if (nbelements > 0)
    {
      nbelements = 0;
      free (elements);
    }

  if (nbpatches > 0)
    {
      nbpatches = 0;
      free (patches);
    }

}

void
MshCorrectNonOrthogonality ()
{

  int i;

  int face, pair;

  int element;
  
  //int neighbor;

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

      faces[face].kj = 0.0;

      /*
      if (pair != -1)
	{

	  neighbor = faces[pair].element;

	  faces[face].kj +=
	    GeoMagVector (GeoSubVectorVector
			  (elements[element].celement, faces[face].rpl));

	  faces[face].kj +=
	    GeoMagVector (GeoSubVectorVector
			  (elements[neighbor].celement, faces[pair].rpl));


	}
      */
    }


}

void
MshGetShapeFunctions () //网格形状判断
{

  int i;

  int element;

  //double detJ;

  double x14, y14, z14;
  double x24, y24, z24;
  double x34, y34, z34;

  double bi, bj, bk, bl;
  double ci, cj, ck, cl;
  double di, dj, dk, dl;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type == TETRAHEDRON)
	{

	  x14 =
	    nodes[elements[element].node[0]].x -
	    nodes[elements[element].node[3]].x;
	  y14 =
	    nodes[elements[element].node[0]].y -
	    nodes[elements[element].node[3]].y;
	  z14 =
	    nodes[elements[element].node[0]].z -
	    nodes[elements[element].node[3]].z;

	  x24 =
	    nodes[elements[element].node[1]].x -
	    nodes[elements[element].node[3]].x;
	  y24 =
	    nodes[elements[element].node[1]].y -
	    nodes[elements[element].node[3]].y;
	  z24 =
	    nodes[elements[element].node[1]].z -
	    nodes[elements[element].node[3]].z;

	  x34 =
	    nodes[elements[element].node[2]].x -
	    nodes[elements[element].node[3]].x;
	  y34 =
	    nodes[elements[element].node[2]].y -
	    nodes[elements[element].node[3]].y;
	  z34 =
	    nodes[elements[element].node[2]].z -
	    nodes[elements[element].node[3]].z;

	  /*
	     detJ = (x14 * y24 * z34 + y14 * z24 * x34 + z14 * x24 * y34) - 
	     (x34 * y24 * z14 + y34 * z24 * x14 + z34 * x24 * y14);

	     elements[element].Vp = LABS(detJ) / 6.0;
	   */

	  bi = (y24 * z34 - y34 * z24);
	  bj = (y34 * z14 - y14 * z34);
	  bk = (y14 * z24 - y24 * z14);
	  bl = -(bi + bj + bk);

	  ci = (z24 * x34 - z34 * x24);
	  cj = (z34 * x14 - z14 * x34);
	  ck = (z14 * x24 - z24 * x14);
	  cl = -(ci + cj + ck);

	  di = (x24 * y34 - x34 * y24);
	  dj = (x34 * y14 - x14 * y34);
	  dk = (x14 * y24 - x24 * y14);
	  dl = -(di + dj + dk);

	  // Allocate memory

	  elements[element].b =
	    calloc (elements[element].nbnodes, sizeof (double));
	  elements[element].c =
	    calloc (elements[element].nbnodes, sizeof (double));
	  elements[element].d =
	    calloc (elements[element].nbnodes, sizeof (double));

	  elements[element].b[0] = bi;
	  elements[element].b[1] = bj;
	  elements[element].b[2] = bk;
	  elements[element].b[3] = bl;

	  elements[element].c[0] = ci;
	  elements[element].c[1] = cj;
	  elements[element].c[2] = ck;
	  elements[element].c[3] = cl;

	  elements[element].d[0] = di;
	  elements[element].d[1] = dj;
	  elements[element].d[2] = dk;
	  elements[element].d[3] = dl;

	}
    }

}

void
MshGetElementTypes ()
{


  int i;

  int element, face, pair;

  nbtris = 0;
  nbquads = 0;
  nbtetras = 0;
  nbhexas = 0;
  nbprisms = 0;

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      pair = faces[face].pair;

      if (faces[face].type == TRIANGLE)
	{
	  nbtris++;
	}

      if (faces[face].type == QUADRANGLE)
	{
	  nbquads++;
	}

    }

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type == TETRAHEDRON)
	{
	  nbtetras++;
	}

      if (elements[element].type == HEXAHEDRON)
	{
	  nbhexas++;
	}

      if (elements[element].type == PRISM)
	{
	  nbprisms++;
	}
    }

}

void
MshCreateFaces ()
{

  int i, j;

  int element;

  msh_vector d;

  // Create faces
  faces = realloc (faces, nbelements * 8 * sizeof (msh_face));

  nbfaces = 0;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type == TRIANGLE)
	{

	  elements[element].nbfaces = 1;

	  // Allocate memory
	  elements[element].face =
	    calloc (elements[element].nbfaces, sizeof (int));

	  elements[element].face[0] = nbfaces;

	  faces[nbfaces].type = TRIANGLE;
	  faces[nbfaces].nbnodes = 3;

	  // Allocate memory
	  faces[nbfaces].node = calloc (faces[nbfaces].nbnodes, sizeof (int));

	  faces[nbfaces].node[0] = elements[element].node[0];
	  faces[nbfaces].node[1] = elements[element].node[1];
	  faces[nbfaces].node[2] = elements[element].node[2];

	  // Calculate area
	  faces[nbfaces].Aj = GeoCalcTriArea (nodes[faces[nbfaces].node[0]],
					      nodes[faces[nbfaces].node[1]],
					      nodes[faces[nbfaces].node[2]]);

	  // Calculate centroid
	  faces[nbfaces].cface =
	    GeoCalcCentroid3 (nodes[faces[nbfaces].node[0]],
			      nodes[faces[nbfaces].node[1]],
			      nodes[faces[nbfaces].node[2]]);

	  // Calculate normal of the face                         
	  faces[nbfaces].n = GeoCalcNormal (nodes[faces[nbfaces].node[0]],
					    nodes[faces[nbfaces].node[1]],
					    nodes[faces[nbfaces].node[2]]);

	  faces[nbfaces].element = element;
	  faces[nbfaces].physreg = -1;
	  faces[nbfaces].elemreg = -1;
	  faces[nbfaces].partition = -1;
	  faces[nbfaces].pair = -1;
	  faces[nbfaces].bc = NONE;

	  nbfaces++;

	}

      if (elements[element].type == QUADRANGLE)
	{

	  elements[element].nbfaces = 1;

	  // Allocate memory
	  elements[element].face =
	    calloc (elements[element].nbfaces, sizeof (int));

	  elements[element].face[0] = nbfaces;

	  faces[nbfaces].type = QUADRANGLE;
	  faces[nbfaces].nbnodes = 4;

	  // Allocate memory
	  faces[nbfaces].node = calloc (faces[nbfaces].nbnodes, sizeof (int));

	  faces[nbfaces].node[0] = elements[element].node[0];
	  faces[nbfaces].node[1] = elements[element].node[1];
	  faces[nbfaces].node[2] = elements[element].node[2];
	  faces[nbfaces].node[3] = elements[element].node[3];

	  // Calculate area
	  faces[nbfaces].Aj = GeoCalcQuadArea (nodes[faces[nbfaces].node[0]],
					       nodes[faces[nbfaces].node[1]],
					       nodes[faces[nbfaces].node[2]],
					       nodes[faces[nbfaces].node[3]]);

	  // Calculate centroid                   
	  faces[nbfaces].cface =
	    GeoCalcCentroid4 (nodes[faces[nbfaces].node[0]],
			      nodes[faces[nbfaces].node[1]],
			      nodes[faces[nbfaces].node[2]],
			      nodes[faces[nbfaces].node[3]]);

	  // Calculate normal of the face                         
	  faces[nbfaces].n = GeoCalcNormal (nodes[faces[nbfaces].node[0]],
					    nodes[faces[nbfaces].node[1]],
					    nodes[faces[nbfaces].node[2]]);

	  faces[nbfaces].element = element;
	  faces[nbfaces].physreg = -1;
	  faces[nbfaces].elemreg = -1;
	  faces[nbfaces].partition = -1;
	  faces[nbfaces].pair = -1;
	  faces[nbfaces].bc = NONE;

	  nbfaces++;

	}

      if (elements[element].type == TETRAHEDRON)
	{

	  elements[element].nbfaces = 4;

	  // Allocate memory
	  elements[element].face =
	    calloc (elements[element].nbfaces, sizeof (int));

	  for (j = 0; j < 4; j++)
	    {

	      elements[element].face[j] = nbfaces;

	      faces[nbfaces].type = TRIANGLE;
	      faces[nbfaces].nbnodes = 3;

	      switch (j)
		{
		case 0:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[0];
		  faces[nbfaces].node[1] = elements[element].node[1];
		  faces[nbfaces].node[2] = elements[element].node[2];

		  break;

		case 1:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[1];
		  faces[nbfaces].node[1] = elements[element].node[3];
		  faces[nbfaces].node[2] = elements[element].node[2];

		  break;

		case 2:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[2];
		  faces[nbfaces].node[1] = elements[element].node[3];
		  faces[nbfaces].node[2] = elements[element].node[0];

		  break;

		case 3:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[3];
		  faces[nbfaces].node[1] = elements[element].node[1];
		  faces[nbfaces].node[2] = elements[element].node[0];

		  break;
		}

	      // Calculate area
	      faces[nbfaces].Aj =
		GeoCalcTriArea (nodes[faces[nbfaces].node[0]],
				nodes[faces[nbfaces].node[1]],
				nodes[faces[nbfaces].node[2]]);

	      // Calculate centroid
	      faces[nbfaces].cface =
		GeoCalcCentroid3 (nodes[faces[nbfaces].node[0]],
				  nodes[faces[nbfaces].node[1]],
				  nodes[faces[nbfaces].node[2]]);

	      // Calculate normal of the face                         
	      faces[nbfaces].n = GeoCalcNormal (nodes[faces[nbfaces].node[0]],
						nodes[faces[nbfaces].node[1]],
						nodes[faces[nbfaces].
						      node[2]]);

	      d =
		GeoSubVectorVector (faces[nbfaces].cface,
				    elements[element].celement);

	      // Normal should point to neighbour
	      // Flip normal if necessary
	      if (GeoDotVectorVector (faces[nbfaces].n, d) < 0.0)
		{
		  faces[nbfaces].n.x *= -1;
		  faces[nbfaces].n.y *= -1;
		  faces[nbfaces].n.z *= -1;
		}

	      faces[nbfaces].element = element;
	      faces[nbfaces].physreg = -1;
	      faces[nbfaces].elemreg = -1;
              faces[nbfaces].partition = -1;
	      faces[nbfaces].pair = -1;
	      faces[nbfaces].bc = NONE;

	      nbfaces++;

	    }

	}

      if (elements[element].type == HEXAHEDRON)
	{

	  elements[element].nbfaces = 6;

	  // Allocate memory
	  elements[element].face =
	    calloc (elements[element].nbfaces, sizeof (int));

	  for (j = 0; j < 6; j++)
	    {

	      elements[element].face[j] = nbfaces;

	      faces[nbfaces].type = QUADRANGLE;
	      faces[nbfaces].nbnodes = 4;

	      switch (j)
		{
		case 0:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[0];
		  faces[nbfaces].node[1] = elements[element].node[1];
		  faces[nbfaces].node[2] = elements[element].node[2];
		  faces[nbfaces].node[3] = elements[element].node[3];

		  break;

		case 1:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[7];
		  faces[nbfaces].node[1] = elements[element].node[6];
		  faces[nbfaces].node[2] = elements[element].node[5];
		  faces[nbfaces].node[3] = elements[element].node[4];

		  break;

		case 2:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[5];
		  faces[nbfaces].node[1] = elements[element].node[6];
		  faces[nbfaces].node[2] = elements[element].node[2];
		  faces[nbfaces].node[3] = elements[element].node[1];

		  break;

		case 3:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[3];
		  faces[nbfaces].node[1] = elements[element].node[7];
		  faces[nbfaces].node[2] = elements[element].node[4];
		  faces[nbfaces].node[3] = elements[element].node[0];

		  break;

		case 4:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[0];
		  faces[nbfaces].node[1] = elements[element].node[4];
		  faces[nbfaces].node[2] = elements[element].node[5];
		  faces[nbfaces].node[3] = elements[element].node[1];

		  break;

		case 5:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[7];
		  faces[nbfaces].node[1] = elements[element].node[3];
		  faces[nbfaces].node[2] = elements[element].node[2];
		  faces[nbfaces].node[3] = elements[element].node[6];

		  break;

		}

	      // Calculate area
	      faces[nbfaces].Aj =
		GeoCalcQuadArea (nodes[faces[nbfaces].node[0]],
				 nodes[faces[nbfaces].node[1]],
				 nodes[faces[nbfaces].node[2]],
				 nodes[faces[nbfaces].node[3]]);

	      // Calculate centroid   
	      faces[nbfaces].cface =
		GeoCalcCentroid4 (nodes[faces[nbfaces].node[0]],
				  nodes[faces[nbfaces].node[1]],
				  nodes[faces[nbfaces].node[2]],
				  nodes[faces[nbfaces].node[3]]);

	      // Calculate normal of the face                         
	      faces[nbfaces].n = GeoCalcNormal (nodes[faces[nbfaces].node[0]],
						nodes[faces[nbfaces].node[1]],
						nodes[faces[nbfaces].
						      node[2]]);

	      d =
		GeoSubVectorVector (faces[nbfaces].cface,
				    elements[element].celement);

	      // Normal should point to neighbour
	      // Flip normal if necessary
	      if (GeoDotVectorVector (faces[nbfaces].n, d) < 0.0)
		{
		  faces[nbfaces].n.x *= -1;
		  faces[nbfaces].n.y *= -1;
		  faces[nbfaces].n.z *= -1;
		}

	      faces[nbfaces].element = element;
	      faces[nbfaces].physreg = -1;
	      faces[nbfaces].elemreg = -1;
              faces[nbfaces].partition = -1;
	      faces[nbfaces].pair = -1;
	      faces[nbfaces].bc = NONE;

	      nbfaces++;

	    }
	}

      if (elements[element].type == PRISM)
	{

	  elements[element].nbfaces = 5;

	  // Allocate memory
	  elements[element].face =
	    calloc (elements[element].nbfaces, sizeof (int));

	  for (j = 0; j < 2; j++)
	    {

	      elements[element].face[j] = nbfaces;

	      faces[nbfaces].type = TRIANGLE;
	      faces[nbfaces].nbnodes = 3;

	      switch (j)
		{
		case 0:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[0];
		  faces[nbfaces].node[1] = elements[element].node[1];
		  faces[nbfaces].node[2] = elements[element].node[2];

		  break;

		case 1:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[3];
		  faces[nbfaces].node[1] = elements[element].node[5];
		  faces[nbfaces].node[2] = elements[element].node[4];

		  break;

		}

	      // Calculate area
	      faces[nbfaces].Aj =
		GeoCalcTriArea (nodes[faces[nbfaces].node[0]],
				nodes[faces[nbfaces].node[1]],
				nodes[faces[nbfaces].node[2]]);

	      // Calculate centroid   
	      faces[nbfaces].cface =
		GeoCalcCentroid3 (nodes[faces[nbfaces].node[0]],
				  nodes[faces[nbfaces].node[1]],
				  nodes[faces[nbfaces].node[2]]);

	      // Calculate normal of the face                         
	      faces[nbfaces].n = GeoCalcNormal (nodes[faces[nbfaces].node[0]],
						nodes[faces[nbfaces].node[1]],
						nodes[faces[nbfaces].
						      node[2]]);

	      d =
		GeoSubVectorVector (faces[nbfaces].cface,
				    elements[element].celement);

	      // Normal should point to neighbour
	      // Flip normal if necessary
	      if (GeoDotVectorVector (faces[nbfaces].n, d) < 0.0)
		{
		  faces[nbfaces].n.x *= -1;
		  faces[nbfaces].n.y *= -1;
		  faces[nbfaces].n.z *= -1;
		}

	      faces[nbfaces].element = element;
	      faces[nbfaces].physreg = -1;
	      faces[nbfaces].elemreg = -1;
              faces[nbfaces].partition = -1;
	      faces[nbfaces].pair = -1;

	      nbfaces++;

	    }

	  for (j = 2; j < 5; j++)
	    {

	      elements[element].face[j] = nbfaces;

	      faces[nbfaces].type = QUADRANGLE;
	      faces[nbfaces].nbnodes = 4;

	      switch (j)
		{
		case 2:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[0];
		  faces[nbfaces].node[1] = elements[element].node[2];
		  faces[nbfaces].node[2] = elements[element].node[5];
		  faces[nbfaces].node[3] = elements[element].node[3];

		  break;

		case 3:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[1];
		  faces[nbfaces].node[1] = elements[element].node[4];
		  faces[nbfaces].node[2] = elements[element].node[5];
		  faces[nbfaces].node[3] = elements[element].node[2];

		  break;

		case 4:

		  // Allocate memory
		  faces[nbfaces].node =
		    calloc (faces[nbfaces].nbnodes, sizeof (int));

		  faces[nbfaces].node[0] = elements[element].node[0];
		  faces[nbfaces].node[1] = elements[element].node[3];
		  faces[nbfaces].node[2] = elements[element].node[4];
		  faces[nbfaces].node[3] = elements[element].node[1];

		  break;

		}

	      // Calculate area 
	      faces[nbfaces].Aj =
		GeoCalcQuadArea (nodes[faces[nbfaces].node[0]],
				 nodes[faces[nbfaces].node[1]],
				 nodes[faces[nbfaces].node[2]],
				 nodes[faces[nbfaces].node[3]]);

	      // Calculate centroid   
	      faces[nbfaces].cface =
		GeoCalcCentroid4 (nodes[faces[nbfaces].node[0]],
				  nodes[faces[nbfaces].node[1]],
				  nodes[faces[nbfaces].node[2]],
				  nodes[faces[nbfaces].node[3]]);

	      // Calculate normal of the face                         
	      faces[nbfaces].n = GeoCalcNormal (nodes[faces[nbfaces].node[0]],
						nodes[faces[nbfaces].node[1]],
						nodes[faces[nbfaces].
						      node[2]]);

	      d =
		GeoSubVectorVector (faces[nbfaces].cface,
				    elements[element].celement);

	      // Normal should point to neighbour
	      // Flip normal if necessary
	      if (GeoDotVectorVector (faces[nbfaces].n, d) < 0.0)
		{
		  faces[nbfaces].n.x *= -1;
		  faces[nbfaces].n.y *= -1;
		  faces[nbfaces].n.z *= -1;
		}

	      faces[nbfaces].element = element;
	      faces[nbfaces].physreg = -1;
	      faces[nbfaces].elemreg = -1;
              faces[nbfaces].partition = -1;
	      faces[nbfaces].pair = -1;
	      faces[nbfaces].bc = NONE;

	      nbfaces++;

	    }

	}

    }

  faces = realloc (faces, nbfaces * sizeof (msh_face));

}

void
MshAddBoundaryFaces ()
{

  int i;

  nbfaces += nbpatches;

  faces = realloc (faces, nbfaces * sizeof (msh_face));

  for (i = nbfaces - nbpatches; i < nbfaces; i++)
    {

      faces[i] = patches[i - nbfaces + nbpatches];

    }

}

void
MshConnectFaces ()
{

  int i, j, k;

  double min[3], max[3];

  int nb;
  int nbv1, nbv2;
  int face, another_face;
  double dist;
  msh_vector cent[2];

  oct_data *tab = NULL;

  // Find pairs of all faces 
  // Boundary faces have no pair: -1
  // Create an octree with face cfaces

  if (nbfaces == 0)
    return;

  tab = malloc (nbfaces * sizeof (oct_data));

  min[0] = faces[0].cface.x;
  min[1] = faces[0].cface.y;
  min[2] = faces[0].cface.z;

  max[0] = faces[0].cface.x;
  max[1] = faces[0].cface.y;
  max[2] = faces[0].cface.z;

  for (i = 0; i < nbfaces; i++)
    {

      tab[i].x = faces[i].cface.x;
      tab[i].y = faces[i].cface.y;
      tab[i].z = faces[i].cface.z;

      min[0] = LMIN (min[0], faces[i].cface.x);
      min[1] = LMIN (min[1], faces[i].cface.y);
      min[2] = LMIN (min[2], faces[i].cface.z);

      max[0] = LMAX (max[0], faces[i].cface.x);
      max[1] = LMAX (max[1], faces[i].cface.y);
      max[2] = LMAX (max[2], faces[i].cface.z);

    }

  OctCreateOctree (min, max, tab, nbfaces);

  for (i = 0; i < nbleafs; i++)
    {

      nb = leafs[i].nbentities;

      // For each node, we find the two corresponding faces
      for (j = 0; j < nb; j++)
	{

	  face = leafs[i].entities[j];
	  nbv1 = faces[face].nbnodes;
	  cent[0] = faces[face].cface;

	  if (faces[face].pair != -1)
	    continue;

	  for (k = j + 1; k < nb; k++)
	    {

	      another_face = leafs[i].entities[k];
	      nbv2 = faces[another_face].nbnodes;
	      cent[1] = faces[another_face].cface;

	      dist = GeoMagVector (GeoSubVectorVector (cent[0], cent[1]));

	      if (dist < EPSILON && nbv1 == nbv2 && face != another_face)
		{

		  faces[face].pair = another_face;
		  faces[another_face].pair = face;

		  break;
		}


	    }

	  leafs[i].entities[k] = leafs[i].entities[nb - j - 1];

	}
    }

  OctDestroyOctree ();

  free (tab);

}

void
MshRemoveBoundaryFaces ()
{

  int i;

  int face, pair;

  nbfaces -= nbpatches;

  for (i = 0; i < nbfaces; i++)
    {
      face = i;

      pair = faces[face].pair;

      if (pair != -1)
	{
	  faces[face].physreg = faces[pair].physreg;
	  faces[face].elemreg = faces[pair].elemreg;
	  faces[face].partition = faces[pair].partition;
	  faces[face].bc = faces[pair].bc;
	}
    }

  for (i = 0; i < nbfaces; i++)
    {
      face = i;

      faces[face].pair = -1;
    }

  faces = realloc (faces, nbfaces * sizeof (msh_face));

}

void
MshCalcNewElementCenters ()
{


  int i, j;


  int element, face;


  double x[8], y[8], z[8];


  for (i = 0; i < nbelements; i++)

    {


      element = i;


      if (elements[element].type == TETRAHEDRON)

	{


	  // Element centre = centre of sphere connecting face centers

	  elements[element].celement.x = 0.0;

	  elements[element].celement.y = 0.0;

	  elements[element].celement.z = 0.0;


	  for (j = 0; j < elements[element].nbfaces; j++)

	    {


	      face = elements[element].face[j];


	      x[j] = faces[face].cface.x;

	      y[j] = faces[face].cface.y;

	      z[j] = faces[face].cface.z;


	    }


	  elements[element].celement.x =
	    -0.5 * (y[2] * z[3] - y[2] * z[0] - y[0] * z[3] - z[2] * y[3] +
		    z[2] * y[0] + z[0] * y[3]) / (-x[2] * y[1] * z[0] +
						  x[1] * y[0] * z[3] -
						  x[1] * z[0] * y[3] +
						  x[2] * z[1] * y[0] -
						  x[3] * y[1] * z[2] -
						  x[1] * y[2] * z[3] +
						  x[1] * y[2] * z[0] +
						  x[1] * z[2] * y[3] -
						  x[1] * z[2] * y[0] +
						  x[0] * y[2] * z[3] -
						  x[0] * z[2] * y[3] +
						  x[2] * y[1] * z[3] -
						  x[2] * y[0] * z[3] -
						  x[2] * z[1] * y[3] +
						  x[2] * z[0] * y[3] -
						  x[0] * y[1] * z[3] +
						  x[0] * z[1] * y[3] +
						  x[3] * y[1] * z[0] +
						  x[3] * z[2] * y[0] +
						  x[3] * z[1] * y[2] -
						  x[3] * z[1] * y[0] -
						  x[3] * y[2] * z[0] +
						  x[0] * y[1] * z[2] -
						  x[0] * z[1] * y[2]) *
	    (x[1] * x[1] - x[0] * x[0] + y[1] * y[1] - y[0] * y[0] +
	     z[1] * z[1] - z[0] * z[0]) + 0.5 * (-y[1] * z[0] + z[1] * y[0] +
						 y[1] * z[3] - y[0] * z[3] -
						 z[1] * y[3] +
						 z[0] * y[3]) / (-x[2] *
								 y[1] * z[0] +
								 x[1] * y[0] *
								 z[3] -
								 x[1] * z[0] *
								 y[3] +
								 x[2] * z[1] *
								 y[0] -
								 x[3] * y[1] *
								 z[2] -
								 x[1] * y[2] *
								 z[3] +
								 x[1] * y[2] *
								 z[0] +
								 x[1] * z[2] *
								 y[3] -
								 x[1] * z[2] *
								 y[0] +
								 x[0] * y[2] *
								 z[3] -
								 x[0] * z[2] *
								 y[3] +
								 x[2] * y[1] *
								 z[3] -
								 x[2] * y[0] *
								 z[3] -
								 x[2] * z[1] *
								 y[3] +
								 x[2] * z[0] *
								 y[3] -
								 x[0] * y[1] *
								 z[3] +
								 x[0] * z[1] *
								 y[3] +
								 x[3] * y[1] *
								 z[0] +
								 x[3] * z[2] *
								 y[0] +
								 x[3] * z[1] *
								 y[2] -
								 x[3] * z[1] *
								 y[0] -
								 x[3] * y[2] *
								 z[0] +
								 x[0] * y[1] *
								 z[2] -
								 x[0] * z[1] *
								 y[2]) *
	    (x[2] * x[2] - x[0] * x[0] + y[2] * y[2] - y[0] * y[0] +
	     z[2] * z[2] - z[0] * z[0]) - 0.5 * (+y[1] * z[2] - y[1] * z[0] -
						 z[2] * y[0] - z[1] * y[2] +
						 z[1] * y[0] +
						 y[2] * z[0]) / (-x[2] *
								 y[1] * z[0] +
								 x[1] * y[0] *
								 z[3] -
								 x[1] * z[0] *
								 y[3] +
								 x[2] * z[1] *
								 y[0] -
								 x[3] * y[1] *
								 z[2] -
								 x[1] * y[2] *
								 z[3] +
								 x[1] * y[2] *
								 z[0] +
								 x[1] * z[2] *
								 y[3] -
								 x[1] * z[2] *
								 y[0] +
								 x[0] * y[2] *
								 z[3] -
								 x[0] * z[2] *
								 y[3] +
								 x[2] * y[1] *
								 z[3] -
								 x[2] * y[0] *
								 z[3] -
								 x[2] * z[1] *
								 y[3] +
								 x[2] * z[0] *
								 y[3] -
								 x[0] * y[1] *
								 z[3] +
								 x[0] * z[1] *
								 y[3] +
								 x[3] * y[1] *
								 z[0] +
								 x[3] * z[2] *
								 y[0] +
								 x[3] * z[1] *
								 y[2] -
								 x[3] * z[1] *
								 y[0] -
								 x[3] * y[2] *
								 z[0] +
								 x[0] * y[1] *
								 z[2] -
								 x[0] * z[1] *
								 y[2]) *
	    (x[3] * x[3] - x[0] * x[0] + y[3] * y[3] - y[0] * y[0] +
	     z[3] * z[3] - z[0] * z[0]);

	  elements[element].celement.y =
	    +0.5 * (x[2] * z[3] - x[2] * z[0] - x[0] * z[3] - z[2] * x[3] +
		    z[2] * x[0] + z[0] * x[3]) / (-x[2] * y[1] * z[0] +
						  x[1] * y[0] * z[3] -
						  x[1] * z[0] * y[3] +
						  x[2] * z[1] * y[0] -
						  x[3] * y[1] * z[2] -
						  x[1] * y[2] * z[3] +
						  x[1] * y[2] * z[0] +
						  x[1] * z[2] * y[3] -
						  x[1] * z[2] * y[0] +
						  x[0] * y[2] * z[3] -
						  x[0] * z[2] * y[3] +
						  x[2] * y[1] * z[3] -
						  x[2] * y[0] * z[3] -
						  x[2] * z[1] * y[3] +
						  x[2] * z[0] * y[3] -
						  x[0] * y[1] * z[3] +
						  x[0] * z[1] * y[3] +
						  x[3] * y[1] * z[0] +
						  x[3] * z[2] * y[0] +
						  x[3] * z[1] * y[2] -
						  x[3] * z[1] * y[0] -
						  x[3] * y[2] * z[0] +
						  x[0] * y[1] * z[2] -
						  x[0] * z[1] * y[2]) *
	    (x[1] * x[1] - x[0] * x[0] + y[1] * y[1] - y[0] * y[0] +
	     z[1] * z[1] - z[0] * z[0]) - 0.5 * (+x[1] * z[3] - x[1] * z[0] -
						 x[0] * z[3] - z[1] * x[3] +
						 z[1] * x[0] +
						 z[0] * x[3]) / (-x[2] *
								 y[1] * z[0] +
								 x[1] * y[0] *
								 z[3] -
								 x[1] * z[0] *
								 y[3] +
								 x[2] * z[1] *
								 y[0] -
								 x[3] * y[1] *
								 z[2] -
								 x[1] * y[2] *
								 z[3] +
								 x[1] * y[2] *
								 z[0] +
								 x[1] * z[2] *
								 y[3] -
								 x[1] * z[2] *
								 y[0] +
								 x[0] * y[2] *
								 z[3] -
								 x[0] * z[2] *
								 y[3] +
								 x[2] * y[1] *
								 z[3] -
								 x[2] * y[0] *
								 z[3] -
								 x[2] * z[1] *
								 y[3] +
								 x[2] * z[0] *
								 y[3] -
								 x[0] * y[1] *
								 z[3] +
								 x[0] * z[1] *
								 y[3] +
								 x[3] * y[1] *
								 z[0] +
								 x[3] * z[2] *
								 y[0] +
								 x[3] * z[1] *
								 y[2] -
								 x[3] * z[1] *
								 y[0] -
								 x[3] * y[2] *
								 z[0] +
								 x[0] * y[1] *
								 z[2] -
								 x[0] * z[1] *
								 y[2]) *
	    (x[2] * x[2] - x[0] * x[0] + y[2] * y[2] - y[0] * y[0] +
	     z[2] * z[2] - z[0] * z[0]) + 0.5 * (-x[1] * z[0] + x[1] * z[2] -
						 z[2] * x[0] - z[1] * x[2] +
						 x[2] * z[0] +
						 z[1] * x[0]) / (-x[2] *
								 y[1] * z[0] +
								 x[1] * y[0] *
								 z[3] -
								 x[1] * z[0] *
								 y[3] +
								 x[2] * z[1] *
								 y[0] -
								 x[3] * y[1] *
								 z[2] -
								 x[1] * y[2] *
								 z[3] +
								 x[1] * y[2] *
								 z[0] +
								 x[1] * z[2] *
								 y[3] -
								 x[1] * z[2] *
								 y[0] +
								 x[0] * y[2] *
								 z[3] -
								 x[0] * z[2] *
								 y[3] +
								 x[2] * y[1] *
								 z[3] -
								 x[2] * y[0] *
								 z[3] -
								 x[2] * z[1] *
								 y[3] +
								 x[2] * z[0] *
								 y[3] -
								 x[0] * y[1] *
								 z[3] +
								 x[0] * z[1] *
								 y[3] +
								 x[3] * y[1] *
								 z[0] +
								 x[3] * z[2] *
								 y[0] +
								 x[3] * z[1] *
								 y[2] -
								 x[3] * z[1] *
								 y[0] -
								 x[3] * y[2] *
								 z[0] +
								 x[0] * y[1] *
								 z[2] -
								 x[0] * z[1] *
								 y[2]) *
	    (x[3] * x[3] - x[0] * x[0] + y[3] * y[3] - y[0] * y[0] +
	     z[3] * z[3] - z[0] * z[0]);

	  elements[element].celement.z =
	    -0.5 * (x[2] * y[3] - x[2] * y[0] - x[0] * y[3] - y[2] * x[3] +
		    y[2] * x[0] + y[0] * x[3]) / (-x[2] * y[1] * z[0] +
						  x[1] * y[0] * z[3] -
						  x[1] * z[0] * y[3] +
						  x[2] * z[1] * y[0] -
						  x[3] * y[1] * z[2] -
						  x[1] * y[2] * z[3] +
						  x[1] * y[2] * z[0] +
						  x[1] * z[2] * y[3] -
						  x[1] * z[2] * y[0] +
						  x[0] * y[2] * z[3] -
						  x[0] * z[2] * y[3] +
						  x[2] * y[1] * z[3] -
						  x[2] * y[0] * z[3] -
						  x[2] * z[1] * y[3] +
						  x[2] * z[0] * y[3] -
						  x[0] * y[1] * z[3] +
						  x[0] * z[1] * y[3] +
						  x[3] * y[1] * z[0] +
						  x[3] * z[2] * y[0] +
						  x[3] * z[1] * y[2] -
						  x[3] * z[1] * y[0] -
						  x[3] * y[2] * z[0] +
						  x[0] * y[1] * z[2] -
						  x[0] * z[1] * y[2]) *
	    (x[1] * x[1] - x[0] * x[0] + y[1] * y[1] - y[0] * y[0] +
	     z[1] * z[1] - z[0] * z[0]) + 0.5 * (-y[1] * x[3] + x[1] * y[3] -
						 x[1] * y[0] - x[0] * y[3] +
						 y[0] * x[3] +
						 y[1] * x[0]) / (-x[2] *
								 y[1] * z[0] +
								 x[1] * y[0] *
								 z[3] -
								 x[1] * z[0] *
								 y[3] +
								 x[2] * z[1] *
								 y[0] -
								 x[3] * y[1] *
								 z[2] -
								 x[1] * y[2] *
								 z[3] +
								 x[1] * y[2] *
								 z[0] +
								 x[1] * z[2] *
								 y[3] -
								 x[1] * z[2] *
								 y[0] +
								 x[0] * y[2] *
								 z[3] -
								 x[0] * z[2] *
								 y[3] +
								 x[2] * y[1] *
								 z[3] -
								 x[2] * y[0] *
								 z[3] -
								 x[2] * z[1] *
								 y[3] +
								 x[2] * z[0] *
								 y[3] -
								 x[0] * y[1] *
								 z[3] +
								 x[0] * z[1] *
								 y[3] +
								 x[3] * y[1] *
								 z[0] +
								 x[3] * z[2] *
								 y[0] +
								 x[3] * z[1] *
								 y[2] -
								 x[3] * z[1] *
								 y[0] -
								 x[3] * y[2] *
								 z[0] +
								 x[0] * y[1] *
								 z[2] -
								 x[0] * z[1] *
								 y[2]) *
	    (x[2] * x[2] - x[0] * x[0] + y[2] * y[2] - y[0] * y[0] +
	     z[2] * z[2] - z[0] * z[0]) - 0.5 * (+x[1] * y[2] - x[1] * y[0] -
						 y[2] * x[0] - y[1] * x[2] +
						 y[1] * x[0] +
						 x[2] * y[0]) / (-x[2] *
								 y[1] * z[0] +
								 x[1] * y[0] *
								 z[3] -
								 x[1] * z[0] *
								 y[3] +
								 x[2] * z[1] *
								 y[0] -
								 x[3] * y[1] *
								 z[2] -
								 x[1] * y[2] *
								 z[3] +
								 x[1] * y[2] *
								 z[0] +
								 x[1] * z[2] *
								 y[3] -
								 x[1] * z[2] *
								 y[0] +
								 x[0] * y[2] *
								 z[3] -
								 x[0] * z[2] *
								 y[3] +
								 x[2] * y[1] *
								 z[3] -
								 x[2] * y[0] *
								 z[3] -
								 x[2] * z[1] *
								 y[3] +
								 x[2] * z[0] *
								 y[3] -
								 x[0] * y[1] *
								 z[3] +
								 x[0] * z[1] *
								 y[3] +
								 x[3] * y[1] *
								 z[0] +
								 x[3] * z[2] *
								 y[0] +
								 x[3] * z[1] *
								 y[2] -
								 x[3] * z[1] *
								 y[0] -
								 x[3] * y[2] *
								 z[0] +
								 x[0] * y[1] *
								 z[2] -
								 x[0] * z[1] *
								 y[2]) *
	    (x[3] * x[3] - x[0] * x[0] + y[3] * y[3] - y[0] * y[0] +
	     z[3] * z[3] - z[0] * z[0]);


	}


      if (elements[element].type == PRISM)

	{


	  elements[element].celement.x = 0.0;

	  elements[element].celement.y = 0.0;

	  elements[element].celement.z = 0.0;


	  for (j = 0; j < elements[element].nbfaces; j++)

	    {


	      face = elements[element].face[j];


	      x[j] = faces[face].cface.x;

	      y[j] = faces[face].cface.y;

	      z[j] = faces[face].cface.z;


	      elements[element].celement.x += x[j];

	      elements[element].celement.y += y[j];

	      elements[element].celement.z += z[j];


	    }


	  elements[element].celement.x /= elements[element].nbfaces;

	  elements[element].celement.y /= elements[element].nbfaces;

	  elements[element].celement.z /= elements[element].nbfaces;


	}


      if (elements[element].type == HEXAHEDRON)

	{


	  elements[element].celement.x = 0.0;

	  elements[element].celement.y = 0.0;

	  elements[element].celement.z = 0.0;


	  for (j = 0; j < elements[element].nbfaces; j++)

	    {


	      face = elements[element].face[j];


	      x[j] = faces[face].cface.x;

	      y[j] = faces[face].cface.y;

	      z[j] = faces[face].cface.z;


	      elements[element].celement.x += x[j];

	      elements[element].celement.y += y[j];

	      elements[element].celement.z += z[j];


	    }


	  elements[element].celement.x /= elements[element].nbfaces;

	  elements[element].celement.y /= elements[element].nbfaces;

	  elements[element].celement.z /= elements[element].nbfaces;


	}


    }


}

void
MshCalcPropMesh ()
{

  int i, j;

  int element, face, pair;

  double sum_area;
  double sum_volume;

  double total_area;
  double total_volume;

  total_area = 0.0;
  total_volume = 0.0;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      sum_area = 0.0;
      sum_volume = 0.0;

      // Volume calculation using Gauss theorem

      for (j = 0; j < elements[element].nbfaces; j++)
	{

	  face = elements[element].face[j];

	  sum_volume += faces[face].cface.x * faces[face].A.x;

	  pair = faces[face].pair;

	  if (pair == -1)
	    sum_area += faces[face].Aj;

	}

      if (sum_volume <= 0.0)
	{
	  PetscPrintf (PETSC_COMM_WORLD,
		       "\nError: Element with zero or negative volume\n");
	  exit (LOGICAL_ERROR);
	}

      elements[element].Vp = sum_volume;

      total_volume += elements[i].Vp;
      total_area += sum_area;

    }

  PetscPrintf (PETSC_COMM_WORLD, "Sub total surface area: \t\t\t%+.3E %s^2\n",
	       total_area, parameter.ulength);
  PetscPrintf (PETSC_COMM_WORLD, "Sub total volume: \t\t\t\t%+.3E %s^3\n",
	       total_volume, parameter.ulength);

}

void
MshCalcPropFaces ()
{

  int i, j, k;

  int element, face, pair, neighbor;

  int nb_internal_faces;
  int nb_total_faces;

  double distance_mean;
  double northo_mean, shapedev_mean;

  int node1, node2;
  double d, dmin, dmax;
  double aspectratio, aspectratio_min, aspectratio_max;

  nb_internal_faces = 0;
  nb_total_faces = 0;

  northo_mean = 0.0;
  distance_mean = 0.0;
  shapedev_mean = 0.0;

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

      faces[face].rpl = GeoSubVectorVector (faces[face].cface, GeoMultScalarVector (GeoDotVectorVector(GeoSubVectorVector(faces[face].cface, elements[element].celement), faces[face].n), faces[face].n));

      if (pair != -1)
	{

	  neighbor = faces[pair].element;

	  faces[face].rnl =
	    GeoSubVectorVector (faces[face].cface, GeoMultScalarVector (GeoDotVectorVector(GeoSubVectorVector(faces[face].cface, elements[neighbor].celement), faces[face].n), faces[face].n));

	  faces[face].A = GeoMultScalarVector (faces[face].Aj, faces[face].n);

	  faces[face].d = GeoSubVectorVector (faces[face].rnl, faces[face].rpl);

	  faces[face].dj = GeoMagVector (faces[face].d);

          if (faces[face].dj == 0)
	  {
		  PetscPrintf (PETSC_COMM_WORLD, "\nError: Problem with mesh\n");
		  exit (LOGICAL_ERROR);
	  }

	  distance_mean += faces[face].dj;
	  nb_internal_faces++;

	}
      else
	{

	  faces[face].A = GeoMultScalarVector (faces[face].Aj, faces[face].n);

	  faces[face].d = GeoSubVectorVector (faces[face].cface, faces[face].rpl);

	  faces[face].dj = GeoMagVector (faces[face].d);

          if (faces[face].dj == 0)
	  {
		  PetscPrintf (PETSC_COMM_WORLD, "\nError: Problem with mesh\n");
		  exit (LOGICAL_ERROR);
	  }

	}

      // Measure non-orthogonality of the mesh
      northo_mean +=
	GeoMagVector (GeoSubVectorVector
		      (elements[element].celement,
		       faces[face].rpl)) / faces[face].dj;
      nb_total_faces++;

    }

  if (nb_internal_faces > 0)
    distance_mean /= nb_internal_faces;

  for (i = 0; i < nbfaces; i++)
    {
      face = i;

      pair = faces[face].pair;

      if (pair != -1)
	{
	  shapedev_mean +=
	    LABS (faces[face].dj - distance_mean) / distance_mean;

	}

    }

  if (nb_total_faces > 0)
    northo_mean /= nb_total_faces;

  if (nb_internal_faces > 0)
    shapedev_mean /= nb_internal_faces;

  aspectratio_min = GREAT;
  aspectratio_max = SMALL;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      dmin = 1E+30f;
      dmax = 1E-30f;

      for (j = 0; j < elements[element].nbfaces; j++)
	{

	  face = elements[element].face[j];

	  for (k = 0; k < faces[face].nbnodes; k++)
	    {

	      node1 = faces[face].node[(k + 0) % faces[face].nbnodes];
	      node2 = faces[face].node[(k + 1) % faces[face].nbnodes];

	      d =
		(float)
		GeoMagVector (GeoSubVectorVector
			      (nodes[node2], nodes[node1]));

	      dmin = LMIN (dmin, d);
	      dmax = LMAX (dmax, d);

	    }

	}

      aspectratio = dmin / dmax;

      aspectratio_min = LMIN (aspectratio_min, aspectratio);
      aspectratio_max = LMAX (aspectratio_max, aspectratio);

    }

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "Non-orthogonality of the mesh: %f\n",
	       northo_mean);
  PetscPrintf (PETSC_COMM_WORLD, "Shape deviation: %f\n", shapedev_mean);
  PetscPrintf (PETSC_COMM_WORLD, "Aspect ratio (worst - best): %f - %f\n",
	       aspectratio_min, aspectratio_max);

}

int
MshImportMSH (char *file)
{
  int i, j;

  int nnod;
  int nele;

  int eindex;
  int maxindex;
  int etype;

  int ntags;
  int physreg, elemreg, partition;

  int ival;

  float fval;

  FILE *fp;
  char descr[512];

  fp = fopen (file, "r");

  if (fp == NULL)
    {
      PetscPrintf (PETSC_COMM_WORLD, "\nError: Mesh file not found!\n");
      PetscPrintf (PETSC_COMM_WORLD, "%s\n\n", file);
      exit (LOGICAL_ERROR);
    }

  PetscPrintf (PETSC_COMM_WORLD, "\nReading mesh file: %s ...\n", file);

  maxindex = 0;

  do
    {

      do
	{

	  strcpy (descr, "");

	  fscanf (fp, "%s", descr);

	  if (strcmp (descr, "$Nodes") == 0)
	    break;

	  if (strcmp (descr, "$Elements") == 0)
	    break;

	}
      while (!feof (fp));

      if (strcmp (descr, "$Nodes") == 0)
	{
	  fscanf (fp, "%d", &nnod);

	  nodes = realloc (nodes, nnod * sizeof (msh_vector));

	  nod_correlation_malloced = nnod;
	  nod_correlation =
	    realloc (nod_correlation,
		     nod_correlation_malloced * sizeof (int));

	  nbnodes = 0;

	  for (i = 0; i < nnod; i++)
	    {

	      fscanf (fp, "%d", &ival);

	      if (ival >= nod_correlation_malloced)
		{
		  nod_correlation_malloced *= 2;
		  nod_correlation =
		    realloc (nod_correlation,
			     nod_correlation_malloced * sizeof (int));
		}

	      nod_correlation[ival] = i;

	      fscanf (fp, "%f", &fval);
	      nodes[nbnodes].x = fval;

	      fscanf (fp, "%f", &fval);
	      nodes[nbnodes].y = fval;

	      fscanf (fp, "%f", &fval);
	      nodes[nbnodes].z = fval;

	      nbnodes++;

	    }
	}

      if (strcmp (descr, "$Elements") == 0)
	{

	  fscanf (fp, "%d", &nele);

	  patches = realloc (patches, nele * sizeof (msh_face));
	  elements = realloc (elements, nele * sizeof (msh_element));

	  nbpatches = 0;
	  nbelements = 0;

	  for (i = 0; i < nele; i++)
	    {

	      fscanf (fp, "%d %d", &eindex, &etype);

	      if (etype < 1)
		{
		  // Invalid element

		  PetscPrintf (PETSC_COMM_WORLD,
			       "\nError: Invalid element\n");
		  exit (LOGICAL_ERROR);

		}

	      if (etype == 1)
		{

		  // Beam element

		  PetscPrintf (PETSC_COMM_WORLD,
			       "\nError: Invalid element\n");
		  exit (LOGICAL_ERROR);

		}

	      if (etype == 2 || etype == 3)
		{

		  // Shell element or surface element

		  if (etype == 2) ival = 3;
		  if (etype == 3) ival = 4;

		  fscanf (fp, "%d", &ntags);

		  if (ntags != 3)
		  {
			printf ("\nError: Number of tags != 3\n");
		 	exit (LOGICAL_ERROR);
		  }

		  fscanf (fp, "%d %d %d", &physreg, &elemreg, &partition);

		  patches[nbpatches].nbnodes = ival;

		  patches[nbpatches].cface.x = 0.0;
		  patches[nbpatches].cface.y = 0.0;
		  patches[nbpatches].cface.z = 0.0;

		  // Allocate memory
		  patches[nbpatches].node =
		    calloc (patches[nbpatches].nbnodes, sizeof (int));

		  for (j = 0; j < patches[nbpatches].nbnodes; j++)
		    {
		      fscanf (fp, "%d", &ival);

		      patches[nbpatches].node[j] = nod_correlation[ival];

		      patches[nbpatches].cface.x +=
			nodes[patches[nbpatches].node[j]].x;
		      patches[nbpatches].cface.y +=
			nodes[patches[nbpatches].node[j]].y;
		      patches[nbpatches].cface.z +=
			nodes[patches[nbpatches].node[j]].z;

		    }

		  patches[nbpatches].cface.x /= patches[nbpatches].nbnodes;
		  patches[nbpatches].cface.y /= patches[nbpatches].nbnodes;
		  patches[nbpatches].cface.z /= patches[nbpatches].nbnodes;

		  patches[nbpatches].physreg = physreg;
		  patches[nbpatches].elemreg = elemreg;
		  patches[nbpatches].partition = partition;

		  patches[nbpatches].element = -1;
		  patches[nbpatches].pair = -1;
		  patches[nbpatches].bc = NONE;

		  if (etype == 2)
		    patches[nbpatches].type = TRIANGLE;

		  if (etype == 3)
		    patches[nbpatches].type = QUADRANGLE;

		  patches[nbpatches].index = eindex;

		  if (nbelements == 0)
		    maxindex = eindex;

		  nbpatches++;

		}

	      if (etype >= 4 && etype <= 6)
		{

		  // Solid elements or volume element

		  if (etype == 4) ival = 4;
		  if (etype == 5) ival = 8;
		  if (etype == 6) ival = 6;

		  fscanf (fp, "%d", &ntags);

		  if (ntags != 3)
		  {
			printf ("\nError: Number of tags != 3\n");
		 	exit (LOGICAL_ERROR);
		  }

		  fscanf (fp, "%d %d %d", &physreg, &elemreg, &partition);

		  elements[nbelements].nbnodes = ival;

		  elements[nbelements].celement.x = 0.0;
		  elements[nbelements].celement.y = 0.0;
		  elements[nbelements].celement.z = 0.0;

		  // Allocate memory
		  elements[nbelements].node =
		    calloc (elements[nbelements].nbnodes, sizeof (int));

		  for (j = 0; j < elements[nbelements].nbnodes; j++)
		    {
		      fscanf (fp, "%d", &ival);

		      elements[nbelements].node[j] = nod_correlation[ival];

		      elements[nbelements].celement.x +=
			nodes[elements[nbelements].node[j]].x;
		      elements[nbelements].celement.y +=
			nodes[elements[nbelements].node[j]].y;
		      elements[nbelements].celement.z +=
			nodes[elements[nbelements].node[j]].z;

		    }

		  elements[nbelements].celement.x /= elements[nbelements].
		    nbnodes;
		  elements[nbelements].celement.y /= elements[nbelements].
		    nbnodes;
		  elements[nbelements].celement.z /= elements[nbelements].
		    nbnodes;

		  elements[nbelements].physreg = physreg;
		  elements[nbelements].elemreg = elemreg;
		  elements[nbelements].partition = partition;

		  if (etype == 4)
		    elements[nbelements].type = TETRAHEDRON;

		  if (etype == 5)
		    elements[nbelements].type = HEXAHEDRON;

		  if (etype == 6)
		    elements[nbelements].type = PRISM;

		  elements[nbelements].index = eindex;

		  nbelements++;

		}

	      if (etype >= 7)
		{

		  // Unkown element

		  PetscPrintf (PETSC_COMM_WORLD, "\nError: Unkown element\n");
		  exit (LOGICAL_ERROR);

		}

	    }

	}

    }
  while (!feof (fp));

  fclose (fp);

  nodes = realloc (nodes, nbnodes * sizeof (msh_vector));
  patches = realloc (patches, nbpatches * sizeof (msh_face));
  elements = realloc (elements, nbelements * sizeof (msh_element));

  free (nod_correlation);

  // Assign patches for parallel processing

  for (i = 0; i < nbpatches; i++)
    {
      if (patches[i].index > maxindex)
	{
	  patches[i].bc = PROCESSOR;
	}
    }

  PetscPrintf (PETSC_COMM_WORLD, "Done.\n");

  MshCreateFaces ();
  MshAddBoundaryFaces ();
  MshConnectFaces ();
  MshRemoveBoundaryFaces ();
  MshConnectFaces ();
  MshCalcPropFaces ();
  //MshCalcNewElementCenters();

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "Number of nodes: \t\t\t\t%d\n", nbnodes);
  PetscPrintf (PETSC_COMM_WORLD, "Number of faces: \t\t\t\t%d\n", nbfaces);
  PetscPrintf (PETSC_COMM_WORLD, "Number of boundary faces: \t\t\t%d\n",
	       nbpatches);
  PetscPrintf (PETSC_COMM_WORLD, "Number of elements: \t\t\t\t%d\n",
	       nbelements);

  if (nbelements == 0)
    {
      PetscPrintf (PETSC_COMM_WORLD, "\nError: No solid elements.\n");
      exit (LOGICAL_ERROR);
    }

  MshCalcPropMesh ();
  MshGetElementTypes ();
  MshCorrectNonOrthogonality ();
  //MshGetShapeFunctions();

  return LOGICAL_TRUE;

}

int
MshExportMSH (char *file)
{

  int i, j;

  int node, patch, element;

  FILE *fp;

  fp = fopen (file, "w");

  if (fp == NULL)
    return LOGICAL_FALSE;

  fprintf (fp, "$MeshFormat\n");
  fprintf (fp, "2 0 8\n");
  fprintf (fp, "$EndMeshFormat\n");

  fprintf (fp, "$Nodes\n");
  fprintf (fp, "%d\n", nbnodes);

  for (i = 0; i < nbnodes; i++)
    {
      node = i;

      fprintf (fp, "%d %f %f %f", node + 1, nodes[node].x, nodes[node].y, nodes[node].z);
      fprintf (fp, "\n");
    }

  fprintf (fp, "$EndNodes\n");

  fprintf (fp, "$Elements\n");

  fprintf (fp, "%d\n", nbpatches + nbelements);

  for (i = 0; i < nbpatches; i++)
    {

      patch = i;

      fprintf (fp, "%d", patches[patch].index);

      if (patches[patch].type == TRIANGLE)
	fprintf (fp, " %d", 2);

      if (patches[patch].type == QUADRANGLE)
	fprintf (fp, " %d", 3);

      fprintf (fp, " %d", 3);

      fprintf (fp, " %d", patches[patch].physreg);
      fprintf (fp, " %d", patches[patch].elemreg);
      fprintf (fp, " %d", patches[patch].partition);

      for (j = 0; j < patches[patch].nbnodes; j++)
	{
	  fprintf (fp, " %d", patches[patch].node[j] + 1);
	}

      fprintf (fp, "\n");
    }

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      fprintf (fp, "%d", elements[element].index);

      if (elements[element].type == TETRAHEDRON)
	fprintf (fp, " %d", 4);

      if (elements[element].type == HEXAHEDRON)
	fprintf (fp, " %d", 5);

      if (elements[element].type == PRISM)
	fprintf (fp, " %d", 6);

      fprintf (fp, " %d", 3);

      fprintf (fp, " %d", elements[element].physreg);
      fprintf (fp, " %d", elements[element].elemreg);
      fprintf (fp, " %d", elements[element].partition);

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %d", elements[element].node[j] + 1);
	}

      fprintf (fp, "\n");
    }

  fprintf (fp, "$EndElements\n");

  fclose (fp);

  return LOGICAL_TRUE;

}

int
MshExportDecomposedMSH (char *file, int partition, int nbpartitions)
{

  int i, j, n;

  int node, face, pair, patch, element, neighbor;

  float b;

  FILE *fp;

  int *nflag;

  int nbnewnodes;
  int nbnewpatches;
  int nbnewelements;

  msh_vector *newnodes;
  msh_face *newpatches;
  msh_element *newelements;

  int *newindex;

  // Identify elements in partition

  newelements = calloc (nbelements, sizeof (msh_element));

  nbnewelements = 0;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].partition != partition)
	continue;

      newelements[nbnewelements] = elements[element];
      nbnewelements++;

    }

  if (nbnewelements == 0)
    {
      PetscPrintf (PETSC_COMM_WORLD, "\nError: Number of elements is zero\n");
      return LOGICAL_ERROR;
    }

  b = (float) nbnewelements / (float) nbelements *nbpartitions;

  printf ("Partition: %d, Balance: %.2f\n", partition, b);

  newelements = realloc (newelements, nbnewelements * sizeof (msh_element));

  // Remove unconnected and merged nodes

  nflag = calloc (nbnodes, sizeof (int));

  for (i = 0; i < nbnodes; i++)
    {
      node = i;

      nflag[node] = LOGICAL_FALSE;
    }

  for (i = 0; i < nbnewelements; i++)
    {

      element = i;

      for (j = 0; j < newelements[element].nbnodes; j++)
	{
	  node = newelements[element].node[j];

	  nflag[node] = LOGICAL_TRUE;
	}

    }

  newnodes = calloc (nbnodes, sizeof (msh_vector));

  newindex = calloc (nbnodes, sizeof (int));

  nbnewnodes = 0;

  for (i = 0; i < nbnodes; i++)
    {

      node = i;

      if (nflag[node] == LOGICAL_TRUE)
	{

	  newnodes[nbnewnodes] = nodes[node];

	  newindex[node] = nbnewnodes;

	  nbnewnodes++;

	}
    }

  newnodes = realloc (newnodes, nbnewnodes * sizeof (msh_vector));

  // Identify patches in partition

  newpatches = calloc (nbpatches + nbfaces, sizeof (msh_face));

  nbnewpatches = 0;

  for (i = 0; i < nbpatches; i++)
    {

      patch = i;

      n = LOGICAL_TRUE;

      for (j = 0; j < patches[patch].nbnodes; j++)
	{
	  node = patches[patch].node[j];

	  if (nflag[node] == LOGICAL_FALSE)
	    n = LOGICAL_FALSE;
	}

      if (n == LOGICAL_FALSE)
	continue;

      newpatches[nbnewpatches] = patches[patch];
      nbnewpatches++;
    }

  // Add patches which belong to cut 

  n = 0;

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

      if (pair != -1)
	{

	  neighbor = faces[pair].element;

	  if (elements[element].partition != partition)
	    continue;

	  if (elements[element].partition == elements[neighbor].partition)
	    continue;

	  newpatches[nbnewpatches] = faces[face];

	  // Set physical region pointing to neighbor
	  newpatches[nbnewpatches].physreg = elements[neighbor].index;

	  // Set element partition pointing to neighbor region
	  newpatches[nbnewpatches].partition = elements[neighbor].partition;

	  // Set index to the end of the list of elements
	  newpatches[nbnewpatches].index = nbelements + nbpatches + n + 1;

	  n++;

	  nbnewpatches++;

	}

    }

  newpatches = realloc (newpatches, nbnewpatches * sizeof (msh_face));

  // Renumber patches and elements        

  for (i = 0; i < nbnewelements; i++)
    {

      element = i;

      for (j = 0; j < newelements[element].nbnodes; j++)
	{
	  node = newelements[element].node[j];

	  newelements[element].node[j] = newindex[node];

	}
    }

  for (i = 0; i < nbnewpatches; i++)
    {

      patch = i;

      for (j = 0; j < newpatches[patch].nbnodes; j++)
	{
	  node = newpatches[patch].node[j];

	  newpatches[patch].node[j] = newindex[node];

	}
    }

  fp = fopen (file, "w");

  if (fp == NULL)
    return LOGICAL_FALSE;

  fprintf (fp, "$MeshFormat\n");
  fprintf (fp, "2 0 8\n");
  fprintf (fp, "$EndMeshFormat\n");

  fprintf (fp, "$Nodes\n");
  fprintf (fp, "%d\n", nbnewnodes);

  for (i = 0; i < nbnewnodes; i++)
    {
      node = i;

      fprintf (fp, "%d %f %f %f", node + 1, newnodes[node].x,
	       newnodes[node].y, newnodes[node].z);

      fprintf (fp, "\n");
    }

  fprintf (fp, "$EndNodes\n");

  fprintf (fp, "$Elements\n");

  fprintf (fp, "%d\n", nbnewpatches + nbnewelements);

  for (i = 0; i < nbnewpatches; i++)
    {

      patch = i;

      if (newpatches[patch].index > nbelements + nbpatches)
	continue;

      fprintf (fp, "%d", newpatches[patch].index);

      if (newpatches[patch].type == TRIANGLE)
	fprintf (fp, " %d", 2);

      if (newpatches[patch].type == QUADRANGLE)
	fprintf (fp, " %d", 3);

      fprintf (fp, " %d", 3);

      fprintf (fp, " %d", newpatches[patch].physreg);
      fprintf (fp, " %d", newpatches[patch].elemreg);
      fprintf (fp, " %d", newpatches[patch].partition);

      for (j = 0; j < newpatches[patch].nbnodes; j++)
	{
	  fprintf (fp, " %d", newpatches[patch].node[j] + 1);
	}

      fprintf (fp, "\n");
    }

  for (i = 0; i < nbnewelements; i++)
    {

      element = i;

      fprintf (fp, "%d", newelements[element].index);

      if (newelements[element].type == TETRAHEDRON)
	fprintf (fp, " %d", 4);

      if (newelements[element].type == HEXAHEDRON)
	fprintf (fp, " %d", 5);

      if (newelements[element].type == PRISM)
	fprintf (fp, " %d", 6);

      fprintf (fp, " %d", 3);

      fprintf (fp, " %d", newelements[element].physreg);
      fprintf (fp, " %d", newelements[element].elemreg);
      fprintf (fp, " %d", newelements[element].partition);

      for (j = 0; j < newelements[element].nbnodes; j++)
	{
	  fprintf (fp, " %d", newelements[element].node[j] + 1);
	}

      fprintf (fp, "\n");
    }

  for (i = 0; i < nbnewpatches; i++)
    {

      patch = i;

      if (newpatches[patch].index <= nbelements + nbpatches)
	continue;

      fprintf (fp, "%d", newpatches[patch].index);

      if (newpatches[patch].type == TRIANGLE)
	fprintf (fp, " %d", 2);

      if (newpatches[patch].type == QUADRANGLE)
	fprintf (fp, " %d", 3);

      fprintf (fp, " %d", 3);

      fprintf (fp, " %d", newpatches[patch].physreg);
      fprintf (fp, " %d", newpatches[patch].elemreg);
      fprintf (fp, " %d", newpatches[patch].partition);

      for (j = 0; j < newpatches[patch].nbnodes; j++)
	{
	  fprintf (fp, " %d", newpatches[patch].node[j] + 1);
	}

      fprintf (fp, "\n");
    }

  fprintf (fp, "$EndElements\n");

  fclose (fp);

  free (nflag);
  free (newnodes);
  free (newpatches);
  free (newelements);
  free (newindex);

  return LOGICAL_TRUE;

}
