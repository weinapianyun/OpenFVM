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

#include "variables.h"
#include "vector.h"

#include "ioutils.h"
#include "globals.h"
#include "mesh.h"
#include "param.h"
#include "bcond.h"
#include "geocalc.h"
#include "gradient.h"
#include "post.h"

#include "parallel.h"

// Function based on streamFunction.C developed by OpenFOAM
// modified and translated to C by the OpenFVM team in 08/02/2006

void
CalculateStreamFunction (double *streamFunction)
{

  int i, j;

  int face, pair, element;

  int found, finished;

  int *visitedPoint;

  int nVisited, nVisitedOld;

  int bPointFound, pointFound;

  double currentBStream, currentStream;

  msh_vector currentBStreamPoint, currentStreamPoint;
  msh_vector edgeHat;

  double vmin, vmax;

  visitedPoint = calloc (nbnodes, sizeof (int));

  nVisited = 0;
  nVisitedOld = 0;

  finished = LOGICAL_TRUE;

  do
    {

      // Find the boundary face with zero flux. set the stream function
      // to zero on that face

      found = LOGICAL_FALSE;

      // Boundary faces

      for (i = 0; i < nbfaces; i++)
	{

	  face = i;

	  pair = faces[face].pair;

	  if (pair != -1)
	    continue;

	  if (faces[face].bc == EMPTY)
	    continue;

	  if (LABS (V_GetCmp (&uf, face)) < SMALL)
	    {

	      // Zero flux face found
	      found = LOGICAL_TRUE;

	      for (j = 0; j < faces[face].nbnodes; j++)
		{

		  if (visitedPoint[faces[face].node[j]] == 1)
		    {
		      found = LOGICAL_FALSE;
		      break;
		    }

		}

	      if (found == LOGICAL_TRUE)
		{

		  //printf("Zero face: %d\n", face);              

		  for (j = 0; j < faces[face].nbnodes; j++)
		    {
		      streamFunction[faces[face].node[j]] = 0.0;
		      visitedPoint[faces[face].node[j]] = 1;
		      nVisited++;
		    }

		  break;
		}


	    }

	  if (found == LOGICAL_TRUE)
	    break;

	}

      if (found == LOGICAL_FALSE)
	{

	  for (i = 0; i < nbelements; i++)
	    {

	      element = i;

	      found = LOGICAL_TRUE;

	      for (j = 0; j < elements[element].nbnodes; j++)
		{

		  if (visitedPoint[elements[element].node[j]] == 1)
		    {
		      found = LOGICAL_FALSE;
		      break;
		    }

		}

	      if (found == LOGICAL_TRUE)
		{

		  for (j = 0; j < elements[element].nbnodes; j++)
		    {
		      streamFunction[elements[element].node[j]] = 0.0;
		      visitedPoint[elements[element].node[j]] = 1;
		      nVisited++;
		    }

		  break;
		}


	    }

	}

      // Loop through all faces. If one of the points on
      // the face has the streamFunction value different
      // from -1, all points with -1 on that face have the
      // streamFunction value equal to the face flux in
      // that point plus the value in the visited point

      do
	{

	  finished = LOGICAL_TRUE;

	  // Boundary faces

	  for (i = 0; i < nbfaces; i++)
	    {

	      face = i;

	      pair = faces[face].pair;

	      if (pair != -1)
		continue;

	      bPointFound = LOGICAL_FALSE;

	      currentBStream = 0.0;

	      currentBStreamPoint.x = 0.0;
	      currentBStreamPoint.y = 0.0;
	      currentBStreamPoint.z = 0.0;

	      for (j = 0; j < faces[face].nbnodes; j++)
		{

		  // Check if the point has been visited

		  if (visitedPoint[faces[face].node[j]] == 1)
		    {
		      // The point has been visited
		      currentBStream = streamFunction[faces[face].node[j]];
		      currentBStreamPoint = nodes[faces[face].node[j]];

		      bPointFound = LOGICAL_TRUE;

		      break;
		    }


		}

	      if (bPointFound == LOGICAL_TRUE)
		{

		  // Sort out other points on the face                            
		  for (j = 0; j < faces[face].nbnodes; j++)
		    {

		      // Check if the point has been visited
		      if (visitedPoint[faces[face].node[j]] == 0)
			{

			  if (faces[face].bc != EMPTY)
			    {

			      edgeHat =
				GeoSubVectorVector (nodes
						    [faces[face].node[j]],
						    currentBStreamPoint);

			      edgeHat.z = 0.0;

			      edgeHat = GeoNormalizeVector (edgeHat);

			      if (edgeHat.y > VSMALL)
				{
				  visitedPoint[faces[face].node[j]] = 1;
				  nVisited++;

				  streamFunction[faces[face].node[j]] =
				    currentBStream + V_GetCmp (&uf,
							       face) *
				    faces[face].Aj * LSGN (faces[face].n.x);
				}
			      else if (edgeHat.y < -VSMALL)
				{
				  visitedPoint[faces[face].node[j]] = 1;
				  nVisited++;

				  streamFunction[faces[face].node[j]] =
				    currentBStream - V_GetCmp (&uf,
							       face) *
				    faces[face].Aj * LSGN (faces[face].n.x);
				}
			      else
				{
				  if (edgeHat.x > VSMALL)
				    {
				      visitedPoint[faces[face].node[j]] = 1;
				      nVisited++;

				      streamFunction[faces[face].node[j]] =
					currentBStream + V_GetCmp (&uf,
								   face) *
					faces[face].Aj *
					LSGN (faces[face].n.y);
				    }
				  else if (edgeHat.x < -VSMALL)
				    {
				      visitedPoint[faces[face].node[j]] = 1;
				      nVisited++;

				      streamFunction[faces[face].node[j]] =
					currentBStream - V_GetCmp (&uf,
								   face) *
					faces[face].Aj *
					LSGN (faces[face].n.y);
				    }
				}
			    }
			}

		    }


		}
	      else
		{
		  finished = LOGICAL_FALSE;
		}

	    }

	  // Internal faces

	  for (i = 0; i < nbfaces; i++)
	    {

	      face = i;

	      pair = faces[face].pair;

	      if (pair == -1)
		continue;

	      pointFound = LOGICAL_FALSE;

	      currentStream = 0.0;

	      currentStreamPoint.x = 0.0;
	      currentStreamPoint.y = 0.0;
	      currentStreamPoint.z = 0.0;

	      for (j = 0; j < faces[face].nbnodes; j++)
		{

		  // Check if the point has been visited
		  if (visitedPoint[faces[face].node[j]] == 1)
		    {
		      // The point has been visited
		      currentStream = streamFunction[faces[face].node[j]];
		      currentStreamPoint = nodes[faces[face].node[j]];

		      pointFound = LOGICAL_TRUE;

		      break;
		    }

		}

	      if (pointFound == LOGICAL_TRUE)
		{

		  // Sort out other points on the face                            
		  for (j = 0; j < faces[face].nbnodes; j++)
		    {

		      // Check if the point has been visited
		      if (visitedPoint[faces[face].node[j]] == 0)
			{

			  edgeHat =
			    GeoSubVectorVector (nodes[faces[face].node[j]],
						currentStreamPoint);

			  edgeHat.z = 0.0;

			  edgeHat = GeoNormalizeVector (edgeHat);

			  if (edgeHat.y > VSMALL)
			    {
			      visitedPoint[faces[face].node[j]] = 1;
			      nVisited++;

			      streamFunction[faces[face].node[j]] =
				currentStream + V_GetCmp (&uf,
							  face) *
				faces[face].Aj * LSGN (faces[face].n.x);
			    }
			  else if (edgeHat.y < -VSMALL)
			    {
			      visitedPoint[faces[face].node[j]] = 1;
			      nVisited++;

			      streamFunction[faces[face].node[j]] =
				currentStream - V_GetCmp (&uf,
							  face) *
				faces[face].Aj * LSGN (faces[face].n.x);
			    }
			}

		    }


		}
	      else
		{
		  finished = LOGICAL_FALSE;
		}


	    }

	  if (nVisited == nVisitedOld)
	    {

	      //printf("Exhausted a seed. Looking for new seed.\n");

	      break;
	    }
	  else
	    {
	      nVisitedOld = nVisited;
	    }

	}
      while (finished == LOGICAL_FALSE);


    }
  while (finished == LOGICAL_FALSE);

  free (visitedPoint);

  // Get maximum and minimum values       
  vmin = +VGREAT;
  vmax = -VGREAT;

  for (i = 0; i < nbnodes; i++)
    {

      vmin = LMIN (vmin, streamFunction[i]);
      vmax = LMAX (vmax, streamFunction[i]);

    }

  // Normalize stream function

  if (vmax != 0.0)
    {
      for (i = 0; i < nbnodes; i++)
	{

	  streamFunction[i] = (streamFunction[i] - vmin) / (vmax - vmin);
	  //streamFunction[i] /= vmax; 

	}
    }

}

void
WriteProbeViews (FILE * fp, char *var, double curtime)
{

  int i, j, k;

  int node, element;

  double vs;
  double sd;

  double *se, *sf;

  vs = 0.0;

  se = calloc (nbelements, sizeof (double));
  sf = calloc (nbfaces, sizeof (double));  
  
  // Create probe views
  for (k = 0; k < nphi; k++)
    {

      if (parameter.probe[k] == LOGICAL_FALSE)
	continue;

      for (i = 0; i < nbelements; i++)
	{

	  element = i;

	  if (k == 0)
	    vs = V_GetCmp (&xul, element);
	  if (k == 1)
	    vs = V_GetCmp (&xvl, element);
	  if (k == 2)
	    vs = V_GetCmp (&xwl, element);
	  if (k == 3)
	    vs = V_GetCmp (&xpl, element);
	  if (k == 4)
	    vs = V_GetCmp (&xTl, element);
	  if (k == 5)
	    vs = V_GetCmp (&xsl, element);

	  se[element] = vs;

	}

      /*
         for (i = 0; i < nbfaces; i++)
         {

         face = i;

         if (k == 0)
         vs = V_GetCmp (&xuf, face);
         if (k == 1)
         vs = V_GetCmp (&xvf, face);
         if (k == 2)
         vs = V_GetCmp (&xwf, face);
         if (k == 3)
         vs = V_GetCmp (&xpf, face);
         if (k == 4)
         vs = V_GetCmp (&xTf, face);
         if (k == 5)
         vs = V_GetCmp (&xsf, face);

         sf[face] = vs;

         }
       */

      fprintf (fp, "View \"Variable: %c\"{\n", var[k]);

      fprintf (fp, "TIME { %f };\n", curtime);

      /*
         for (i = 0; i < nbfaces; i++)
         {

         face = i;

         if (faces[face].type == TRIANGLE) fprintf (fp, "ST(");
         if (faces[face].type == QUADRANGLE) fprintf (fp, "SQ(");

         for (j = 0; j < faces[face].nbnodes; j++)
         {
         node = faces[face].node[j];

         fprintf (fp, "%f, %f, %f", nodes[node].x, nodes[node].y, nodes[node].z);

         if (j != faces[face].nbnodes - 1)
         fprintf (fp, ",");

         }

         fprintf (fp, ")");

         fprintf (fp, "{");

         for (j = 0; j < faces[face].nbnodes; j++)
         {

         sd = sf[face];

         fprintf (fp, " %f", sd);

         if (j != faces[face].nbnodes - 1)
         fprintf (fp, ",");

         }

         fprintf (fp, "};\n");

         }
       */

      for (i = 0; i < nbelements; i++)
	{

	  element = i;

	  if (elements[element].type == TETRAHEDRON)
	    fprintf (fp, "SS(");
	  if (elements[element].type == HEXAHEDRON)
	    fprintf (fp, "SH(");
	  if (elements[element].type == PRISM)
	    fprintf (fp, "SI(");

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      node = elements[element].node[j];

	      fprintf (fp, "%f, %f, %f", nodes[node].x, nodes[node].y,
		       nodes[node].z);

	      if (j != elements[element].nbnodes - 1)
		fprintf (fp, ",");

	    }

	  fprintf (fp, ")");

	  fprintf (fp, "{");

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {

	      sd = se[element];

	      fprintf (fp, " %f", sd);

	      if (j != elements[element].nbnodes - 1)
		fprintf (fp, ",");

	    }

	  fprintf (fp, "};\n");

	}

      fprintf (fp, "};\n");

    }

  free (se);
  free (sf);

}

void
PrintAsciiHeaderScalarFace (FILE * fp, char *label)
{

  fprintf (fp, "%s %d\n", label, 1);
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-points nb-vector-points nb-tensor-points
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-lines nb-vector-lines nb-tensor-lines
  fprintf (fp, "%d %d %d\n", nbtris, 0, 0);	//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
  fprintf (fp, "%d %d %d\n", nbquads, 0, 0);	//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
  fprintf (fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

}

void
PrintAsciiHeaderScalarElement (FILE * fp, char *label)
{


  fprintf (fp, "%s %d\n", label, 1);
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-points nb-vector-points nb-tensor-points
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-lines nb-vector-lines nb-tensor-lines
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
  fprintf (fp, "%d %d %d\n", nbtetras, 0, 0);	//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
  fprintf (fp, "%d %d %d\n", nbhexas, 0, 0);	//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
  fprintf (fp, "%d %d %d\n", nbprisms, 0, 0);	//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
  fprintf (fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

}

void
PrintAsciiHeaderVectorElement (FILE * fp, char *label)
{


  fprintf (fp, "%s %d\n", label, 1);
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-points nb-vector-points nb-tensor-points
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-lines nb-vector-lines nb-tensor-lines
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
  fprintf (fp, "%d %d %d\n", 0, nbtetras, 0);	//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
  fprintf (fp, "%d %d %d\n", 0, nbhexas, 0);	//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
  fprintf (fp, "%d %d %d\n", 0, nbprisms, 0);	//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
  fprintf (fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

}

void
PrintAsciiHeaderVectorFace (FILE * fp, char *label)
{

  fprintf (fp, "%s %d\n", label, 1);
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-points nb-vector-points nb-tensor-points
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-lines nb-vector-lines nb-tensor-lines
  fprintf (fp, "%d %d %d\n", 0, nbtris, 0);	//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
  fprintf (fp, "%d %d %d\n", 0, nbquads, 0);	//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
  fprintf (fp, "%d %d %d\n", 0, 0, 0);	//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
  fprintf (fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

}

void
PrintBinaryHeaderScalarFace (FILE * fp, char *label)
{

  int one = 1;

  fprintf (fp, "%s %d ", label, 1);	//view-name nb-time-steps,
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-points nb-vector-points nb-tensor-points
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-lines nb-vector-lines nb-tensor-lines
  fprintf (fp, "%d %d %d ", nbtris, 0, 0);	//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
  fprintf (fp, "%d %d %d ", nbquads, 0, 0);	//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
  fprintf (fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

  fwrite (&one, sizeof (int), 1, fp);

}

void
PrintBinaryHeaderScalarElement (FILE * fp, char *label)
{

  int one = 1;

  fprintf (fp, "%s %d ", label, 1);	//view-name nb-time-steps,
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-points nb-vector-points nb-tensor-points
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-lines nb-vector-lines nb-tensor-lines
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
  fprintf (fp, "%d %d %d ", nbtetras, 0, 0);	//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
  fprintf (fp, "%d %d %d ", nbhexas, 0, 0);	//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
  fprintf (fp, "%d %d %d ", nbprisms, 0, 0);	//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
  fprintf (fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

  fwrite (&one, sizeof (int), 1, fp);

}

void
PrintBinaryHeaderVectorElement (FILE * fp, char *label)
{

  int one = 1;

  fprintf (fp, "%s %d ", label, 1);	//view-name nb-time-steps,
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-points nb-vector-points nb-tensor-points
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-lines nb-vector-lines nb-tensor-lines
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
  fprintf (fp, "%d %d %d ", 0, nbtetras, 0);	//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
  fprintf (fp, "%d %d %d ", 0, nbhexas, 0);	//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
  fprintf (fp, "%d %d %d ", 0, nbprisms, 0);	//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
  fprintf (fp, "%d %d %d ", 0, 0, 0);	//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
  fprintf (fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

  fwrite (&one, sizeof (int), 1, fp);

}

void
PrintAsciiTime (FILE * fp, double curtime)
{

  fprintf (fp, "%E\n", curtime);

}

void
PrintBinaryTime (FILE * fp, double curtime)
{

  double *times;

  times = calloc (1, sizeof (double));

  times[0] = curtime;

  fwrite (times, sizeof (double), 1, fp);

  free (times);

}

void
PrintAsciiScalarNode (FILE * fp, double *sn)
{

  int i, j;

  int element;

  double sd;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type != TETRAHEDRON)
	continue;

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].x);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].y);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].z);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{

	  sd = sn[elements[element].node[j]];

	  fprintf (fp, " %E", sd);
	}

      fprintf (fp, "\n");
    }

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type != HEXAHEDRON)
	continue;

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].x);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].y);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].z);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{

	  sd = sn[elements[element].node[j]];

	  fprintf (fp, " %E", sd);
	}

      fprintf (fp, "\n");
    }

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type != PRISM)
	continue;

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].x);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].y);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].z);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{

	  sd = sn[elements[element].node[j]];

	  fprintf (fp, " %E", sd);
	}

      fprintf (fp, "\n");
    }

}

void
PrintAsciiScalarElement (FILE * fp, double *se)
{

  int i, j;

  int element;

  double sd;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type != TETRAHEDRON)
	continue;

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].x);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].y);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].z);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{

	  sd = se[element];

	  fprintf (fp, " %E", sd);
	}

      fprintf (fp, "\n");
    }

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type != HEXAHEDRON)
	continue;

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].x);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].y);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].z);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{

	  sd = se[element];

	  fprintf (fp, " %E", sd);
	}

      fprintf (fp, "\n");
    }

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type != PRISM)
	continue;

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].x);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].y);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].z);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{

	  sd = se[element];

	  fprintf (fp, " %E", sd);
	}

      fprintf (fp, "\n");
    }

}

void
PrintAsciiVectorElement (FILE * fp, int iphix, int iphiy, int iphiz)
{

  int i, j;

  int element;

  double vx, vy, vz;

  vx = 0.0;
  vy = 0.0;
  vz = 0.0;
  
  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type != TETRAHEDRON)
	continue;

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].x);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].y);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].z);
	}

      if (iphix == 0)
	vx = V_GetCmp (&xul, element);
      if (iphix == 1)
	vx = V_GetCmp (&xvl, element);
      if (iphix == 2)
	vx = V_GetCmp (&xwl, element);
      if (iphix == 3)
	vx = V_GetCmp (&xpl, element);
      if (iphix == 4)
	vx = V_GetCmp (&xTl, element);
      if (iphix == 5)
	vx = V_GetCmp (&xsl, element);

      if (iphiy == 0)
	vy = V_GetCmp (&xul, element);
      if (iphiy == 1)
	vy = V_GetCmp (&xvl, element);
      if (iphiy == 2)
	vy = V_GetCmp (&xwl, element);
      if (iphiy == 3)
	vy = V_GetCmp (&xpl, element);
      if (iphiy == 4)
	vy = V_GetCmp (&xTl, element);
      if (iphiy == 5)
	vy = V_GetCmp (&xsl, element);

      if (iphiz == 0)
	vz = V_GetCmp (&xul, element);
      if (iphiz == 1)
	vz = V_GetCmp (&xvl, element);
      if (iphiz == 2)
	vz = V_GetCmp (&xwl, element);
      if (iphiz == 3)
	vz = V_GetCmp (&xpl, element);
      if (iphiz == 4)
	vz = V_GetCmp (&xTl, element);
      if (iphiz == 5)
	vz = V_GetCmp (&xsl, element);

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E %E %E", vx, vy, vz);
	}

      fprintf (fp, "\n");
    }

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type != HEXAHEDRON)
	continue;

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].x);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].y);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].z);
	}

      if (iphix == 0)
	vx = V_GetCmp (&xul, element);
      if (iphix == 1)
	vx = V_GetCmp (&xvl, element);
      if (iphix == 2)
	vx = V_GetCmp (&xwl, element);
      if (iphix == 3)
	vx = V_GetCmp (&xpl, element);
      if (iphix == 4)
	vx = V_GetCmp (&xTl, element);
      if (iphix == 5)
	vx = V_GetCmp (&xsl, element);

      if (iphiy == 0)
	vy = V_GetCmp (&xul, element);
      if (iphiy == 1)
	vy = V_GetCmp (&xvl, element);
      if (iphiy == 2)
	vy = V_GetCmp (&xwl, element);
      if (iphiy == 3)
	vy = V_GetCmp (&xpl, element);
      if (iphiy == 4)
	vy = V_GetCmp (&xTl, element);
      if (iphiy == 5)
	vy = V_GetCmp (&xsl, element);

      if (iphiz == 0)
	vz = V_GetCmp (&xul, element);
      if (iphiz == 1)
	vz = V_GetCmp (&xvl, element);
      if (iphiz == 2)
	vz = V_GetCmp (&xwl, element);
      if (iphiz == 3)
	vz = V_GetCmp (&xpl, element);
      if (iphiz == 4)
	vz = V_GetCmp (&xTl, element);
      if (iphiz == 5)
	vz = V_GetCmp (&xsl, element);

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E %E %E", vx, vy, vz);
	}

      fprintf (fp, "\n");
    }

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type != PRISM)
	continue;

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].x);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].y);
	}

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[elements[element].node[j]].z);
	}

      if (iphix == 0)
	vx = V_GetCmp (&xul, element);
      if (iphix == 1)
	vx = V_GetCmp (&xvl, element);
      if (iphix == 2)
	vx = V_GetCmp (&xwl, element);
      if (iphix == 3)
	vx = V_GetCmp (&xpl, element);
      if (iphix == 4)
	vx = V_GetCmp (&xTl, element);
      if (iphix == 5)
	vx = V_GetCmp (&xsl, element);

      if (iphiy == 0)
	vy = V_GetCmp (&xul, element);
      if (iphiy == 1)
	vy = V_GetCmp (&xvl, element);
      if (iphiy == 2)
	vy = V_GetCmp (&xwl, element);
      if (iphiy == 3)
	vy = V_GetCmp (&xpl, element);
      if (iphiy == 4)
	vy = V_GetCmp (&xTl, element);
      if (iphiy == 5)
	vy = V_GetCmp (&xsl, element);

      if (iphiz == 0)
	vz = V_GetCmp (&xul, element);
      if (iphiz == 1)
	vz = V_GetCmp (&xvl, element);
      if (iphiz == 2)
	vz = V_GetCmp (&xwl, element);
      if (iphiz == 3)
	vz = V_GetCmp (&xpl, element);
      if (iphiz == 4)
	vz = V_GetCmp (&xTl, element);
      if (iphiz == 5)
	vz = V_GetCmp (&xsl, element);

      for (j = 0; j < elements[element].nbnodes; j++)
	{
	  fprintf (fp, " %E %E %E", vx, vy, vz);
	}

      fprintf (fp, "\n");
    }

}

void
PrintAsciiScalarFace (FILE * fp, double *sf)
{

  int i, j;

  int face;

  double sd;

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      if (faces[face].type != TRIANGLE)
	continue;

      for (j = 0; j < faces[face].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[faces[face].node[j]].x);
	}

      for (j = 0; j < faces[face].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[faces[face].node[j]].y);
	}

      for (j = 0; j < faces[face].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[faces[face].node[j]].z);
	}

      for (j = 0; j < faces[face].nbnodes; j++)
	{

	  sd = sf[face];

	  fprintf (fp, " %E", sd);
	}

      fprintf (fp, "\n");
    }

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      if (faces[face].type != QUADRANGLE)
	continue;

      for (j = 0; j < faces[face].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[faces[face].node[j]].x);
	}

      for (j = 0; j < faces[face].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[faces[face].node[j]].y);
	}

      for (j = 0; j < faces[face].nbnodes; j++)
	{
	  fprintf (fp, " %E", nodes[faces[face].node[j]].z);
	}

      for (j = 0; j < faces[face].nbnodes; j++)
	{

	  sd = sf[face];

	  fprintf (fp, " %E", sd);
	}

      fprintf (fp, "\n");
    }

}

void
PrintBinaryScalarNode (FILE * fp, double *sn)
{

  int i, j;

  int element;

  double *tetras;
  double *hexas;
  double *prisms;

  int nt, nh, np;

  tetras = calloc (4 * nbtetras * 4, sizeof (double));
  hexas = calloc (4 * nbhexas * 8, sizeof (double));
  prisms = calloc (4 * nbprisms * 6, sizeof (double));

  nt = 0;
  nh = 0;
  np = 0;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type == TETRAHEDRON)
	{

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      tetras[nt] = nodes[elements[element].node[j]].x;
	      nt++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      tetras[nt] = nodes[elements[element].node[j]].y;
	      nt++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      tetras[nt] = nodes[elements[element].node[j]].z;
	      nt++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {

	      tetras[nt] = sn[elements[element].node[j]];
	      nt++;

	    }
	}

      if (elements[element].type == HEXAHEDRON)
	{

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      hexas[nh] = nodes[elements[element].node[j]].x;
	      nh++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      hexas[nh] = nodes[elements[element].node[j]].y;
	      nh++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      hexas[nh] = nodes[elements[element].node[j]].z;
	      nh++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {

	      hexas[nh] = sn[elements[element].node[j]];
	      nh++;

	    }
	}

      if (elements[element].type == PRISM)
	{

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      prisms[np] = nodes[elements[element].node[j]].x;
	      np++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      prisms[np] = nodes[elements[element].node[j]].y;
	      np++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      prisms[np] = nodes[elements[element].node[j]].z;
	      np++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {

	      prisms[np] = sn[elements[element].node[j]];
	      np++;

	    }
	}
    }

  fwrite (tetras, sizeof (double), nt, fp);
  fwrite (hexas, sizeof (double), nh, fp);
  fwrite (prisms, sizeof (double), np, fp);

  free (tetras);
  free (hexas);
  free (prisms);

}

void
PrintBinaryScalarFace (FILE * fp, double *sf)
{

  int i, j;

  int face;

  double *tris;
  double *quads;

  int nt, nq;

  tris = calloc (4 * nbtris * 3, sizeof (double));
  quads = calloc (4 * nbquads * 4, sizeof (double));

  nt = 0;
  nq = 0;

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      if (faces[face].type == TRIANGLE)
	{

	  for (j = 0; j < faces[face].nbnodes; j++)
	    {
	      tris[nt] = nodes[faces[face].node[j]].x;
	      nt++;
	    }

	  for (j = 0; j < faces[face].nbnodes; j++)
	    {
	      tris[nt] = nodes[faces[face].node[j]].y;
	      nt++;
	    }

	  for (j = 0; j < faces[face].nbnodes; j++)
	    {
	      tris[nt] = nodes[faces[face].node[j]].z;
	      nt++;
	    }

	  for (j = 0; j < faces[face].nbnodes; j++)
	    {

	      tris[nt] = sf[face];

	      nt++;

	    }
	}

      if (faces[face].type == QUADRANGLE)
	{

	  for (j = 0; j < faces[face].nbnodes; j++)
	    {
	      quads[nq] = nodes[faces[face].node[j]].x;
	      nq++;
	    }

	  for (j = 0; j < faces[face].nbnodes; j++)
	    {
	      quads[nq] = nodes[faces[face].node[j]].y;
	      nq++;
	    }

	  for (j = 0; j < faces[face].nbnodes; j++)
	    {
	      quads[nq] = nodes[faces[face].node[j]].z;
	      nq++;
	    }

	  for (j = 0; j < faces[face].nbnodes; j++)
	    {

	      quads[nq] = sf[face];

	      nq++;

	    }
	}

    }

  fwrite (tris, sizeof (double), nt, fp);
  fwrite (quads, sizeof (double), nq, fp);

  free (tris);
  free (quads);

}

void
PrintBinaryScalarElement (FILE * fp, double *se)
{

  int i, j;

  int element;

  double *tetras;
  double *hexas;
  double *prisms;

  int nt, nh, np;

  tetras = calloc (4 * nbtetras * 4, sizeof (double));
  hexas = calloc (4 * nbhexas * 8, sizeof (double));
  prisms = calloc (4 * nbprisms * 6, sizeof (double));

  nt = 0;
  nh = 0;
  np = 0;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type == TETRAHEDRON)
	{

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      tetras[nt] = nodes[elements[element].node[j]].x;
	      nt++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      tetras[nt] = nodes[elements[element].node[j]].y;
	      nt++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      tetras[nt] = nodes[elements[element].node[j]].z;
	      nt++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {

	      tetras[nt] = se[element];

	      nt++;

	    }
	}

      if (elements[element].type == HEXAHEDRON)
	{

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      hexas[nh] = nodes[elements[element].node[j]].x;
	      nh++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      hexas[nh] = nodes[elements[element].node[j]].y;
	      nh++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      hexas[nh] = nodes[elements[element].node[j]].z;
	      nh++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {

	      hexas[nh] = se[element];

	      nh++;

	    }
	}

      if (elements[element].type == PRISM)
	{

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      prisms[np] = nodes[elements[element].node[j]].x;
	      np++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      prisms[np] = nodes[elements[element].node[j]].y;
	      np++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      prisms[np] = nodes[elements[element].node[j]].z;
	      np++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {

	      prisms[np] = se[element];

	      np++;

	    }
	}
    }

  fwrite (tetras, sizeof (double), nt, fp);
  fwrite (hexas, sizeof (double), nh, fp);
  fwrite (prisms, sizeof (double), np, fp);

  free (tetras);
  free (hexas);
  free (prisms);

}

void
PrintBinaryVectorElement (FILE * fp, int iphix, int iphiy, int iphiz)
{

  int i, j;

  int element;

  double vx, vy, vz;
  
  double *tetras;
  double *hexas;
  double *prisms;

  int nt, nh, np;

  vx = 0.0;
  vy = 0.0;
  vz = 0.0;
  
  tetras = calloc (6 * nbtetras * 4, sizeof (double));
  hexas = calloc (6 * nbhexas * 8, sizeof (double));
  prisms = calloc (6 * nbprisms * 6, sizeof (double));

  nt = 0;
  nh = 0;
  np = 0;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (elements[element].type == TETRAHEDRON)
	{

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      tetras[nt] = nodes[elements[element].node[j]].x;
	      nt++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      tetras[nt] = nodes[elements[element].node[j]].y;
	      nt++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      tetras[nt] = nodes[elements[element].node[j]].z;
	      nt++;
	    }

	  if (iphix == 0)
	    vx = V_GetCmp (&xul, element);
	  if (iphix == 1)
	    vx = V_GetCmp (&xvl, element);
	  if (iphix == 2)
	    vx = V_GetCmp (&xwl, element);
	  if (iphix == 3)
	    vx = V_GetCmp (&xpl, element);
	  if (iphix == 4)
	    vx = V_GetCmp (&xTl, element);
	  if (iphix == 5)
	    vx = V_GetCmp (&xsl, element);

	  if (iphiy == 0)
	    vy = V_GetCmp (&xul, element);
	  if (iphiy == 1)
	    vy = V_GetCmp (&xvl, element);
	  if (iphiy == 2)
	    vy = V_GetCmp (&xwl, element);
	  if (iphiy == 3)
	    vy = V_GetCmp (&xpl, element);
	  if (iphiy == 4)
	    vy = V_GetCmp (&xTl, element);
	  if (iphiy == 5)
	    vy = V_GetCmp (&xsl, element);

	  if (iphiz == 0)
	    vz = V_GetCmp (&xul, element);
	  if (iphiz == 1)
	    vz = V_GetCmp (&xvl, element);
	  if (iphiz == 2)
	    vz = V_GetCmp (&xwl, element);
	  if (iphiz == 3)
	    vz = V_GetCmp (&xpl, element);
	  if (iphiz == 4)
	    vz = V_GetCmp (&xTl, element);
	  if (iphiz == 5)
	    vz = V_GetCmp (&xsl, element);

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {

	      tetras[nt] = vx;
	      nt++;

	      tetras[nt] = vy;
	      nt++;

	      tetras[nt] = vz;
	      nt++;

	    }
	}

      if (elements[element].type == HEXAHEDRON)
	{

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      hexas[nh] = nodes[elements[element].node[j]].x;
	      nh++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      hexas[nh] = nodes[elements[element].node[j]].y;
	      nh++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      hexas[nh] = nodes[elements[element].node[j]].z;
	      nh++;
	    }

	  if (iphix == 0)
	    vx = V_GetCmp (&xul, element);
	  if (iphix == 1)
	    vx = V_GetCmp (&xvl, element);
	  if (iphix == 2)
	    vx = V_GetCmp (&xwl, element);
	  if (iphix == 3)
	    vx = V_GetCmp (&xpl, element);
	  if (iphix == 4)
	    vx = V_GetCmp (&xTl, element);
	  if (iphix == 5)
	    vx = V_GetCmp (&xsl, element);

	  if (iphiy == 0)
	    vy = V_GetCmp (&xul, element);
	  if (iphiy == 1)
	    vy = V_GetCmp (&xvl, element);
	  if (iphiy == 2)
	    vy = V_GetCmp (&xwl, element);
	  if (iphiy == 3)
	    vy = V_GetCmp (&xpl, element);
	  if (iphiy == 4)
	    vy = V_GetCmp (&xTl, element);
	  if (iphiy == 5)
	    vy = V_GetCmp (&xsl, element);

	  if (iphiz == 0)
	    vz = V_GetCmp (&xul, element);
	  if (iphiz == 1)
	    vz = V_GetCmp (&xvl, element);
	  if (iphiz == 2)
	    vz = V_GetCmp (&xwl, element);
	  if (iphiz == 3)
	    vz = V_GetCmp (&xpl, element);
	  if (iphiz == 4)
	    vz = V_GetCmp (&xTl, element);
	  if (iphiz == 5)
	    vz = V_GetCmp (&xsl, element);

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {

	      hexas[nh] = vx;
	      nh++;

	      hexas[nh] = vy;
	      nh++;

	      hexas[nh] = vz;
	      nh++;

	    }
	}

      if (elements[element].type == PRISM)
	{

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      prisms[np] = nodes[elements[element].node[j]].x;
	      np++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      prisms[np] = nodes[elements[element].node[j]].y;
	      np++;
	    }

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {
	      prisms[np] = nodes[elements[element].node[j]].z;
	      np++;
	    }

	  if (iphix == 0)
	    vx = V_GetCmp (&xul, element);
	  if (iphix == 1)
	    vx = V_GetCmp (&xvl, element);
	  if (iphix == 2)
	    vx = V_GetCmp (&xwl, element);
	  if (iphix == 3)
	    vx = V_GetCmp (&xpl, element);
	  if (iphix == 4)
	    vx = V_GetCmp (&xTl, element);
	  if (iphix == 5)
	    vx = V_GetCmp (&xsl, element);

	  if (iphiy == 0)
	    vy = V_GetCmp (&xul, element);
	  if (iphiy == 1)
	    vy = V_GetCmp (&xvl, element);
	  if (iphiy == 2)
	    vy = V_GetCmp (&xwl, element);
	  if (iphiy == 3)
	    vy = V_GetCmp (&xpl, element);
	  if (iphiy == 4)
	    vy = V_GetCmp (&xTl, element);
	  if (iphiy == 5)
	    vy = V_GetCmp (&xsl, element);

	  if (iphiz == 0)
	    vz = V_GetCmp (&xul, element);
	  if (iphiz == 1)
	    vz = V_GetCmp (&xvl, element);
	  if (iphiz == 2)
	    vz = V_GetCmp (&xwl, element);
	  if (iphiz == 3)
	    vz = V_GetCmp (&xpl, element);
	  if (iphiz == 4)
	    vz = V_GetCmp (&xTl, element);
	  if (iphiz == 5)
	    vz = V_GetCmp (&xsl, element);

	  for (j = 0; j < elements[element].nbnodes; j++)
	    {

	      prisms[np] = vx;
	      np++;

	      prisms[np] = vy;
	      np++;

	      prisms[np] = vz;
	      np++;

	    }
	}
    }

  fwrite (tetras, sizeof (double), nt, fp);
  fwrite (hexas, sizeof (double), nh, fp);
  fwrite (prisms, sizeof (double), np, fp);

  free (tetras);
  free (hexas);
  free (prisms);

}

void
WriteAsciiScalarElement (FILE * fp, char *label, int iphi, double curtime)
{

  int i;

  int element;

  double vs;

  double *se;

  vs = 0.0;
  
  fprintf (fp, "$View\n");

  se = calloc (nbelements, sizeof (double));

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (iphi == 0)
	vs = V_GetCmp (&xul, element);
      if (iphi == 1)
	vs = V_GetCmp (&xvl, element);
      if (iphi == 2)
	vs = V_GetCmp (&xwl, element);
      if (iphi == 3)
	vs = V_GetCmp (&xpl, element);
      if (iphi == 4)
	vs = V_GetCmp (&xTl, element);
      if (iphi == 5)
	vs = V_GetCmp (&xsl, element);

      se[element] = vs;

    }

  PrintAsciiHeaderScalarElement (fp, label);
  PrintAsciiTime (fp, curtime);
  PrintAsciiScalarElement (fp, se);

  fprintf (fp, "$EndView\n");

  free (se);

}

void
WriteAsciiVectorElement (FILE * fp, char *label, int iphix, int iphiy,
			 int iphiz, double curtime)
{

  fprintf (fp, "$View\n");

  PrintAsciiHeaderVectorElement (fp, label);
  PrintAsciiTime (fp, curtime);
  PrintAsciiVectorElement (fp, iphix, iphiy, iphiz);

  fprintf (fp, "$EndView\n");

}

void
WriteAsciiStreamFunction (FILE * fp, char *label, double curtime)
{

  double *sn;

  fprintf (fp, "$View\n");

  sn = calloc (nbnodes, sizeof (double));

  CalculateStreamFunction (sn);

  PrintAsciiHeaderScalarElement (fp, label);
  PrintAsciiTime (fp, curtime);
  PrintAsciiScalarNode (fp, sn);

  fprintf (fp, "$EndView\n");

  free (sn);

}


void
WriteAsciiVorticity (FILE * fp, char *label, int axis, double curtime)
{

  int i;

  int element;

  double vs;

  double *se;
  
  msh_vector gradu, gradv, gradw;

  vs = 0.0;
  
  fprintf (fp, "$View\n");

  se = calloc (nbelements, sizeof (double));

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      gradu = Gradient (&xul, &xuf, LOGICAL_FALSE, element);
      gradv = Gradient (&xvl, &xvf, LOGICAL_FALSE, element);
      gradw = Gradient (&xwl, &xwf, LOGICAL_FALSE, element);

      if (axis == 0)
	{
	  vs = gradw.y - gradv.z;
	}

      if (axis == 1)
	{
	  vs = gradu.z - gradw.x;
	}

      if (axis == 2)
	{
	  vs = gradv.x - gradu.y;
	}

      se[element] = vs;

    }

  PrintAsciiHeaderScalarElement (fp, label);
  PrintAsciiTime (fp, curtime);
  PrintAsciiScalarElement (fp, se);

  fprintf (fp, "\n$EndView\n");

  free (se);

}

void
WriteAsciiScalarFace (FILE * fp, char *label, int iphi, double curtime)
{

  int i;

  int face;

  double vs;

  double *sf;

  vs = 0.0;
  
  fprintf (fp, "$View\n");

  sf = calloc (nbfaces, sizeof (double));

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      if (iphi == 0)
	vs = V_GetCmp (&xuf, face);
      if (iphi == 1)
	vs = V_GetCmp (&xvf, face);
      if (iphi == 2)
	vs = V_GetCmp (&xwf, face);
      if (iphi == 3)
	vs = V_GetCmp (&xpf, face);
      if (iphi == 4)
	vs = V_GetCmp (&xTf, face);
      if (iphi == 5)
	vs = V_GetCmp (&xsf, face);

      sf[face] = vs;

    }

  PrintAsciiHeaderScalarFace (fp, label);
  PrintAsciiTime (fp, curtime);
  PrintAsciiScalarFace (fp, sf);

  fprintf (fp, "\n$EndView\n");

  free (sf);

}

void
WriteAsciiVectorMagnitudeFace (FILE * fp, char *label,
			       int iphix, int iphiy, int iphiz,
			       double curtime)
{

  int i;

  int face;

  double vx, vy, vz;
  
  double *sf;

  vx = 0.0;
  vy = 0.0;
  vz = 0.0;
  
  sf = calloc (nbfaces, sizeof (double));

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      if (iphix == 0)
	vx = V_GetCmp (&xuf, face);
      if (iphix == 1)
	vx = V_GetCmp (&xvf, face);
      if (iphix == 2)
	vx = V_GetCmp (&xwf, face);
      if (iphix == 3)
	vx = V_GetCmp (&xpf, face);
      if (iphix == 4)
	vx = V_GetCmp (&xTf, face);
      if (iphix == 5)
	vx = V_GetCmp (&xsf, face);

      if (iphiy == 0)
	vy = V_GetCmp (&xuf, face);
      if (iphiy == 1)
	vy = V_GetCmp (&xvf, face);
      if (iphiy == 2)
	vy = V_GetCmp (&xwf, face);
      if (iphiy == 3)
	vy = V_GetCmp (&xpf, face);
      if (iphiy == 4)
	vy = V_GetCmp (&xTf, face);
      if (iphiy == 5)
	vy = V_GetCmp (&xsf, face);

      if (iphiz == 0)
	vz = V_GetCmp (&xuf, face);
      if (iphiz == 1)
	vz = V_GetCmp (&xvf, face);
      if (iphiz == 2)
	vz = V_GetCmp (&xwf, face);
      if (iphiz == 3)
	vz = V_GetCmp (&xpf, face);
      if (iphiz == 4)
	vz = V_GetCmp (&xTf, face);
      if (iphiz == 5)
	vz = V_GetCmp (&xsf, face);

      sf[face] = sqrt (vx * vx + vy * vy + vz * vz);

    }

  fprintf (fp, "$View\n");

  PrintAsciiHeaderScalarFace (fp, label);
  PrintAsciiTime (fp, curtime);
  PrintAsciiScalarFace (fp, sf);

  fprintf (fp, "\n$EndView\n");

  free (sf);

}

void
WriteAsciiFlux (FILE * fp, char *label, double curtime)
{

  int i;

  int face;

  double vs;

  double *sf;

  vs = 0.0;
  
  fprintf (fp, "$View\n");

  sf = calloc (nbfaces, sizeof (double));

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      vs = V_GetCmp (&uf, face);

      sf[face] = vs;

    }

  PrintAsciiHeaderScalarFace (fp, label);
  PrintAsciiTime (fp, curtime);
  PrintAsciiScalarFace (fp, sf);

  fprintf (fp, "\n$EndView\n");

  free (sf);

}

void
WriteBinaryFlux (FILE * fp, char *label, double curtime)
{

  int i;

  int face;

  double vs;

  double *sf;

  vs = 0.0;
  
  fprintf (fp, "$View\n");

  sf = calloc (nbfaces, sizeof (double));

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      vs = V_GetCmp (&uf, face);

      sf[face] = vs;

    }

  PrintBinaryHeaderScalarFace (fp, label);
  PrintBinaryTime (fp, curtime);
  PrintBinaryScalarFace (fp, sf);

  fprintf (fp, "\n$EndView\n");

  free (sf);

}

void
WriteBinaryScalarFace (FILE * fp, char *label, int iphi, double curtime)
{

  int i;

  int face;

  double vs;

  double *sf;

  vs = 0.0;
  
  fprintf (fp, "$View\n");

  sf = calloc (nbfaces, sizeof (double));

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      if (iphi == 0)
	vs = V_GetCmp (&xuf, face);
      if (iphi == 1)
	vs = V_GetCmp (&xvf, face);
      if (iphi == 2)
	vs = V_GetCmp (&xwf, face);
      if (iphi == 3)
	vs = V_GetCmp (&xpf, face);
      if (iphi == 4)
	vs = V_GetCmp (&xTf, face);
      if (iphi == 5)
	vs = V_GetCmp (&xsf, face);

      sf[face] = vs;

    }

  PrintBinaryHeaderScalarFace (fp, label);
  PrintBinaryTime (fp, curtime);
  PrintBinaryScalarFace (fp, sf);

  fprintf (fp, "\n$EndView\n");

  free (sf);

}

void
WriteBinaryVectorMagnitudeFace (FILE * fp, char *label,
				int iphix, int iphiy, int iphiz,
				double curtime)
{

  int i;

  int face;

  double vx, vy, vz;

  double *sf;
  
  vx = 0.0;
  vy = 0.0;
  vz = 0.0;

  sf = calloc (nbfaces, sizeof (double));
  
  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      if (iphix == 0)
	vx = V_GetCmp (&xuf, face);
      if (iphix == 1)
	vx = V_GetCmp (&xvf, face);
      if (iphix == 2)
	vx = V_GetCmp (&xwf, face);
      if (iphix == 3)
	vx = V_GetCmp (&xpf, face);
      if (iphix == 4)
	vx = V_GetCmp (&xTf, face);
      if (iphix == 5)
	vx = V_GetCmp (&xsf, face);

      if (iphiy == 0)
	vy = V_GetCmp (&xuf, face);
      if (iphiy == 1)
	vy = V_GetCmp (&xvf, face);
      if (iphiy == 2)
	vy = V_GetCmp (&xwf, face);
      if (iphiy == 3)
	vy = V_GetCmp (&xpf, face);
      if (iphiy == 4)
	vy = V_GetCmp (&xTf, face);
      if (iphiy == 5)
	vy = V_GetCmp (&xsf, face);

      if (iphiz == 0)
	vz = V_GetCmp (&xuf, face);
      if (iphiz == 1)
	vz = V_GetCmp (&xvf, face);
      if (iphiz == 2)
	vz = V_GetCmp (&xwf, face);
      if (iphiz == 3)
	vz = V_GetCmp (&xpf, face);
      if (iphiz == 4)
	vz = V_GetCmp (&xTf, face);
      if (iphiz == 5)
	vz = V_GetCmp (&xsf, face);

      sf[face] = sqrt (vx * vx + vy * vy + vz * vz);

    }

  fprintf (fp, "$View\n");

  PrintBinaryHeaderScalarFace (fp, label);
  PrintBinaryTime (fp, curtime);
  PrintBinaryScalarFace (fp, sf);

  fprintf (fp, "\n$EndView\n");

  free (sf);

}

void
WriteBinaryScalarElement (FILE * fp, char *label, int iphi, double curtime)
{

  int i;

  int element;

  double vs;

  double *se;

  vs = 0.0;
  
  fprintf (fp, "$View\n");

  se = calloc (nbelements, sizeof (double));

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (iphi == 0)
	vs = V_GetCmp (&xul, element);
      if (iphi == 1)
	vs = V_GetCmp (&xvl, element);
      if (iphi == 2)
	vs = V_GetCmp (&xwl, element);
      if (iphi == 3)
	vs = V_GetCmp (&xpl, element);
      if (iphi == 4)
	vs = V_GetCmp (&xTl, element);
      if (iphi == 5)
	vs = V_GetCmp (&xsl, element);

      se[element] = vs;

    }

  PrintBinaryHeaderScalarElement (fp, label);
  PrintBinaryTime (fp, curtime);
  PrintBinaryScalarElement (fp, se);

  fprintf (fp, "\n$EndView\n");

  free (se);

}

void
WriteBinaryVectorElement (FILE * fp, char *label, int iphix, int iphiy,
			  int iphiz, double curtime)
{

  fprintf (fp, "$View\n");

  PrintBinaryHeaderVectorElement (fp, label);
  PrintBinaryTime (fp, curtime);
  PrintBinaryVectorElement (fp, iphix, iphiy, iphiz);

  fprintf (fp, "\n$EndView\n");

}

void
WriteBinaryStreamFunction (FILE * fp, char *label, double curtime)
{

  double *sn;

  fprintf (fp, "$View\n");

  sn = calloc (nbnodes, sizeof (double));

  CalculateStreamFunction (sn);

  PrintBinaryHeaderScalarElement (fp, label);
  PrintBinaryTime (fp, curtime);
  PrintBinaryScalarNode (fp, sn);

  fprintf (fp, "\n$EndView\n");

  free (sn);
}

void
WriteBinaryVorticity (FILE * fp, char *label, int axis, double curtime)
{

  int i;

  int element;

  double omega = 0.0; 

  double *se;

  msh_vector gradu, gradv, gradw;

  fprintf (fp, "$View\n");

  se = calloc (nbelements, sizeof (double));

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      gradu = Gradient (&xul, &xuf, LOGICAL_FALSE, element);
      gradv = Gradient (&xvl, &xvf, LOGICAL_FALSE, element);
      gradw = Gradient (&xwl, &xwf, LOGICAL_FALSE, element);

      if (axis == 0)
	{
	  omega = gradw.y - gradv.z;
	}

      if (axis == 1)
	{
	  omega = gradu.z - gradw.x;
	}

      if (axis == 2)
	{
	  omega = gradv.x - gradu.y;
	}

      se[element] = omega;

    }

  PrintBinaryHeaderScalarElement (fp, label);
  PrintBinaryTime (fp, curtime);
  PrintBinaryScalarElement (fp, se);

  fprintf (fp, "\n$EndView\n");

  free (se);

}

void
WriteResults (FILE * fp, int onFace, int inElement, double curtime)
{

  VecGhostGetLocalForm (xu, &xul);
  VecGhostGetLocalForm (xv, &xvl);
  VecGhostGetLocalForm (xw, &xwl);
  VecGhostGetLocalForm (xp, &xpl);
  VecGhostGetLocalForm (xT, &xTl);
  VecGhostGetLocalForm (xs, &xsl);

  if (parameter.savflux == LOGICAL_TRUE)
    {
      if (parameter.wbinary == LOGICAL_TRUE)
	{
	  sprintf (viewname, "Flux_[%s/%s]-%03d", parameter.ulength,
		   parameter.utime, processor);
	  WriteBinaryFlux (fp, viewname, curtime);
	}
      else
	{
	  sprintf (viewname, "Flux_[%s/%s]-%03d", parameter.ulength,
		   parameter.utime, processor);
	  WriteAsciiFlux (fp, viewname, curtime);
	}
    }

  if (onFace == LOGICAL_TRUE)
    {

      if (parameter.wbinary == LOGICAL_TRUE)
	{

	  // Write velocity-U to file 
	  sprintf (viewname, "VelocityU(Face)_[%s/%s]-%03d",
		   parameter.ulength, parameter.utime, processor);
	  if (parameter.fsav[iu] == LOGICAL_TRUE)
	    WriteBinaryScalarFace (fp, viewname, iu, curtime);

	  // Write velocity-V to file 
	  sprintf (viewname, "VelocityV(Face)_[%s/%s]-%03d",
		   parameter.ulength, parameter.utime, processor);
	  if (parameter.fsav[iv] == LOGICAL_TRUE)
	    WriteBinaryScalarFace (fp, viewname, iv, curtime);

	  // Write velocity-W to file 
	  sprintf (viewname, "VelocityW(Face)_[%s/%s]-%03d",
		   parameter.ulength, parameter.utime, processor);
	  if (parameter.fsav[iw] == LOGICAL_TRUE)
	    WriteBinaryScalarFace (fp, viewname, iw, curtime);

	  // Write pressure to file
	  sprintf (viewname, "Pressure(Face)_[%s/(%s.%s.%s)]-%03d",
		   parameter.umass, parameter.ulength, parameter.utime,
		   parameter.utime, processor);
	  if (parameter.fsav[ip] == LOGICAL_TRUE)
	    WriteBinaryScalarFace (fp, viewname, ip, curtime);

	  // Write temperature to file
	  sprintf (viewname, "Temperature(Face)_[%s]-%03d",
		   parameter.utemperature, processor);
	  if (parameter.fsav[iT] == LOGICAL_TRUE)
	    WriteBinaryScalarFace (fp, viewname, iT, curtime);

	  // Write indicator function to file
	  sprintf (viewname, "Gamma(Face)_[-]-%03d", processor);
	  if (parameter.fsav[is] == LOGICAL_TRUE)
	    WriteBinaryScalarFace (fp, viewname, is, curtime);

	  // Write velocity to file 
	  sprintf (viewname, "Velocity(Face)_[%s/%s]-%03d", parameter.ulength,
		   parameter.utime, processor);
	  if (parameter.fvec == LOGICAL_TRUE)
	    WriteBinaryVectorMagnitudeFace (fp, viewname, iu, iv, iw,
					    curtime);

	}
      else
	{

	  // Write velocity-U to file 
	  sprintf (viewname, "VelocityU(Face)_[%s/%s]-%03d",
		   parameter.ulength, parameter.utime, processor);
	  if (parameter.fsav[iu] == LOGICAL_TRUE)
	    WriteAsciiScalarFace (fp, viewname, iu, curtime);

	  // Write velocity-V to file 
	  sprintf (viewname, "VelocityV(Face)_[%s/%s]-%03d",
		   parameter.ulength, parameter.utime, processor);
	  if (parameter.fsav[iv] == LOGICAL_TRUE)
	    WriteAsciiScalarFace (fp, viewname, iv, curtime);

	  // Write velocity-W to file 
	  sprintf (viewname, "VelocityW(Face)_[%s/%s]-%03d",
		   parameter.ulength, parameter.utime, processor);
	  if (parameter.fsav[iw] == LOGICAL_TRUE)
	    WriteAsciiScalarFace (fp, viewname, iw, curtime);

	  // Write pressure to file
	  sprintf (viewname, "Pressure(Face)_[%s/(%s.%s.%s)]-%03d",
		   parameter.umass, parameter.ulength, parameter.utime,
		   parameter.utime, processor);
	  if (parameter.fsav[ip] == LOGICAL_TRUE)
	    WriteAsciiScalarFace (fp, viewname, ip, curtime);

	  // Write temperature to file
	  sprintf (viewname, "Temperature(Face)_[%s]-%03d",
		   parameter.utemperature, processor);
	  if (parameter.fsav[iT] == LOGICAL_TRUE)
	    WriteAsciiScalarFace (fp, viewname, iT, curtime);

	  // Write indicator function to file
	  sprintf (viewname, "Gamma(Face)_[-]-%03d", processor);
	  if (parameter.fsav[is] == LOGICAL_TRUE)
	    WriteAsciiScalarFace (fp, viewname, is, curtime);

	  // Write velocity to file 
	  sprintf (viewname, "Velocity(Face)_[%s/%s]-%03d", parameter.ulength,
		   parameter.utime, processor);
	  if (parameter.fvec == LOGICAL_TRUE)
	    WriteBinaryVectorMagnitudeFace (fp, viewname, iu, iv, iw,
					    curtime);

	}
    }

  if (inElement == LOGICAL_TRUE)
    {

      if (parameter.wbinary == LOGICAL_TRUE)
	{

	  // Write velocity-U to file 
	  sprintf (viewname, "VelocityU(Cell)_[%s/%s]-%03d",
		   parameter.ulength, parameter.utime, processor);
	  if (parameter.csav[iu] == LOGICAL_TRUE)
	    WriteBinaryScalarElement (fp, viewname, iu, curtime);

	  // Write velocity-V to file 
	  sprintf (viewname, "VelocityV(Cell)_[%s/%s]-%03d",
		   parameter.ulength, parameter.utime, processor);
	  if (parameter.csav[iv] == LOGICAL_TRUE)
	    WriteBinaryScalarElement (fp, viewname, iv, curtime);

	  // Write velocity-W to file 
	  sprintf (viewname, "VelocityW(Cell)_[%s/%s]-%03d",
		   parameter.ulength, parameter.utime, processor);
	  if (parameter.csav[iw] == LOGICAL_TRUE)
	    WriteBinaryScalarElement (fp, viewname, iw, curtime);

	  // Write pressure to file
	  sprintf (viewname, "Pressure(Cell)_[%s/(%s.%s.%s)]-%03d",
		   parameter.umass, parameter.ulength, parameter.utime,
		   parameter.utime, processor);
	  if (parameter.csav[ip] == LOGICAL_TRUE)
	    WriteBinaryScalarElement (fp, viewname, ip, curtime);

	  // Write indicator function to file
	  sprintf (viewname, "Temperature(Cell)_[%s]-%03d",
		   parameter.utemperature, processor);
	  if (parameter.csav[iT] == LOGICAL_TRUE)
	    WriteBinaryScalarElement (fp, viewname, iT, curtime);

	  // Write indicator function to file
	  sprintf (viewname, "Gamma(Cell)_[-]-%03d", processor);
	  if (parameter.csav[is] == LOGICAL_TRUE)
	    WriteBinaryScalarElement (fp, viewname, is, curtime);

	  // Write velocity to file 
	  sprintf (viewname, "Velocity(Cell)_[%s/%s]-%03d", parameter.ulength,
		   parameter.utime, processor);
	  if (parameter.cvec == LOGICAL_TRUE)
	    WriteBinaryVectorElement (fp, viewname, iu, iv, iw, curtime);

	  // Write vorticity to file                      
	  sprintf (viewname, "VorticityU(Cell)_[1/%s]-%03d", parameter.utime,
		   processor);
	  if (parameter.vortex[0] == LOGICAL_TRUE)
	    WriteBinaryVorticity (fp, viewname, 0, curtime);

	  sprintf (viewname, "VorticityV(Cell)_[1/%s]-%03d", parameter.utime,
		   processor);
	  if (parameter.vortex[1] == LOGICAL_TRUE)
	    WriteBinaryVorticity (fp, viewname, 1, curtime);

	  sprintf (viewname, "VorticityW(Cell)_[1/%s]-%03d", parameter.utime,
		   processor);
	  if (parameter.vortex[2] == LOGICAL_TRUE)
	    WriteBinaryVorticity (fp, viewname, 2, curtime);

	  sprintf (viewname, "StreamFunction(Cell)_[-]-%03d", processor);
	  if (parameter.streamf == LOGICAL_TRUE)
	    WriteBinaryStreamFunction (fp, viewname, curtime);

	}
      else
	{

	  // Write velocity-U to file 
	  sprintf (viewname, "VelocityU(Cell)_[%s/%s]-%03d",
		   parameter.ulength, parameter.utime, processor);
	  if (parameter.csav[iu] == LOGICAL_TRUE)
	    WriteAsciiScalarElement (fp, viewname, iu, curtime);

	  // Write velocity-V to file 
	  sprintf (viewname, "VelocityV(Cell)_[%s/%s]-%03d",
		   parameter.ulength, parameter.utime, processor);
	  if (parameter.csav[iv] == LOGICAL_TRUE)
	    WriteAsciiScalarElement (fp, viewname, iv, curtime);

	  // Write velocity-W to file                     
	  sprintf (viewname, "VelocityW(Cell)_[%s/%s]-%03d",
		   parameter.ulength, parameter.utime, processor);
	  if (parameter.csav[iw] == LOGICAL_TRUE)
	    WriteAsciiScalarElement (fp, viewname, iw, curtime);

	  // Write pressure to file
	  sprintf (viewname, "Pressure(Cell)_[%s/(%s.%s.%s)]-%03d",
		   parameter.umass, parameter.ulength, parameter.utime,
		   parameter.utime, processor);
	  if (parameter.csav[ip] == LOGICAL_TRUE)
	    WriteAsciiScalarElement (fp, viewname, ip, curtime);

	  // Write indicator function to file
	  sprintf (viewname, "Temperature(Cell)_[%s]-%03d",
		   parameter.utemperature, processor);
	  if (parameter.csav[iT] == LOGICAL_TRUE)
	    WriteAsciiScalarElement (fp, viewname, iT, curtime);

	  // Write indicator function to file
	  sprintf (viewname, "Gamma(Cell)_[-]-%03d", processor);
	  if (parameter.csav[is] == LOGICAL_TRUE)
	    WriteAsciiScalarElement (fp, viewname, is, curtime);

	  // Write velocity to file 
	  sprintf (viewname, "Velocity(Cell)_[%s/%s]-%03d", parameter.ulength,
		   parameter.utime, processor);
	  if (parameter.cvec == LOGICAL_TRUE)
	    WriteAsciiVectorElement (fp, viewname, iu, iv, iw, curtime);

	  // Write vorticity to file                      
	  sprintf (viewname, "VorticityU(Cell)_[1/s]-%03d", processor);
	  if (parameter.vortex[0] == LOGICAL_TRUE)
	    WriteAsciiVorticity (fp, viewname, 0, curtime);

	  sprintf (viewname, "VorticityV(Cell)_[1/s]-%03d", processor);
	  if (parameter.vortex[1] == LOGICAL_TRUE)
	    WriteAsciiVorticity (fp, viewname, 1, curtime);

	  sprintf (viewname, "VorticityW(Cell)_[1/s]-%03d", processor);
	  if (parameter.vortex[2] == LOGICAL_TRUE)
	    WriteAsciiVorticity (fp, viewname, 2, curtime);

	  sprintf (viewname, "StreamFunction(Cell)_[-]-%03d", processor);
	  if (parameter.streamf == LOGICAL_TRUE)
	    WriteAsciiStreamFunction (fp, viewname, curtime);

	}

    }

  VecGhostRestoreLocalForm (xu, &xul);
  VecGhostRestoreLocalForm (xv, &xvl);
  VecGhostRestoreLocalForm (xw, &xwl);
  VecGhostRestoreLocalForm (xp, &xpl);
  VecGhostRestoreLocalForm (xT, &xTl);
  VecGhostRestoreLocalForm (xs, &xsl);

}

void
WriteResidual (FILE * fp, int iter, double *res)
{

  int i;

  PetscFPrintf (PETSC_COMM_WORLD, fp, "%d", iter);

  for (i = 0; i < nphi; i++)
    {
      PetscFPrintf (PETSC_COMM_WORLD, fp, " \t %E", res[i]);
    }

  PetscFPrintf (PETSC_COMM_WORLD, fp, "\n");

}
