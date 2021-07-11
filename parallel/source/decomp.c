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
#include <malloc.h>
#include <string.h>

#include "petsc.h"

#include "mesh.h"
#include "globals.h"
#include "decomp.h"

void
DecomposeMesh (char *path, int nbpartitions)
{

  int i, j, l, m, n;

  int nx, ny;
  double dxmin, dxmax, dymin, dymax;

  int element, neighbor, face, pair;

  int freeAdj;

  int newindex;

  int *xadj;
  int *adjncy;

  char *file;

  // No weight information provided
  int wgtFlag = 0;

  // C style numbering
  int numFlag = 0;

  // Decomposition options. options[0] = 0 - use defaults
  int options[5];

  // Processor weights
  float *weights;

  // Final decomposition
  int *decomp;

  // Number of edges cut
  int edgeCut = 0;

  printf ("\nDecomposing mesh into %d partitions...\n", nbpartitions);

  xadj = calloc (nbelements + 1, sizeof (int));

  n = 0;

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      pair = faces[face].pair;

      if (pair != -1)
	{
	  n++;
	}
    }

  adjncy = calloc (n, sizeof (int));

  freeAdj = 0;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      xadj[element] = freeAdj;

      for (j = 0; j < elements[element].nbfaces; j++)
	{

	  face = elements[element].face[j];

	  pair = faces[face].pair;

	  if (pair != -1)
	    {

	      neighbor = faces[pair].element;

	      adjncy[freeAdj++] = neighbor;

	    }
	}
    }

  xadj[nbelements] = freeAdj;

  // Use even weighting
  weights = calloc (nbpartitions, sizeof (float));

  for (j = 0; j < nbpartitions; j++)
    {
      weights[j] = 1.0 / (float) nbpartitions;
    }

  decomp = calloc (nbelements, sizeof (int));

  if (1 == 1)
    {
      // Decompose mesh
      if (nbpartitions > 8)
	{
	  options[0] = 1;
	  options[1] = 3;
	  options[2] = 1;
	  options[3] = 3;
	  options[4] = 0;
	  METIS_WPartGraphKway (&nbelements, xadj, adjncy, 3, 0, &wgtFlag,
				&numFlag, &nbpartitions, weights, options,
				&edgeCut, decomp);
	}
      else
	{
	  if (nbpartitions > 1)
	    {
	      options[0] = 1;
	      options[1] = 3;
	      options[2] = 1;
	      options[3] = 1;
	      options[4] = 0;
	      METIS_WPartGraphRecursive (&nbelements, xadj, adjncy, 3, 0,
					 &wgtFlag, &numFlag, &nbpartitions,
					 weights, options, &edgeCut, decomp);
	    }
	  else
	    {
	      for (i = 0; i < nbelements; i++)
		{
		  element = i;

		  decomp[element] = 0;
		}
	    }
	}
    }
  else
    {

      // Decompose grid
      nx = sqrt (nbpartitions);
      ny = sqrt (nbpartitions);

      for (i = 0; i < nbelements; i++)
	{

	  element = i;

	  for (l = 0; l < nx; l++)
	    {

	      dxmin = (float) (l + 0) / (float) nx;
	      dxmax = (float) (l + 1) / (float) nx;

	      if (elements[element].celement.x < dxmin
		  || elements[element].celement.x > dxmax)
		continue;

	      for (m = 0; m < ny; m++)
		{

		  dymin = (float) (m + 0) / (float) ny;
		  dymax = (float) (m + 1) / (float) ny;

		  if (elements[element].celement.y < dymin
		      || elements[element].celement.y > dymax)
		    continue;

		  decomp[element] = l * nx + m;

		  //printf("l * nx + m: %d \n", l * nx + m);

		}

	    }

	}

      edgeCut = nx * ny;

    }

  for (i = 0; i < nbelements; i++)
    {
      element = i;

      elements[element].partition = decomp[element];
    }

  // Create new indexes PETSc style
  newindex = 0;

  for (j = 0; j < nbpartitions; j++)
    {

      for (i = 0; i < nbelements; i++)
	{
	  element = i;

	  if (elements[element].partition == j)
	    {
	      elements[element].index = newindex;
	      newindex++;
	    }
	}

    }

  // Export mesh files
  file = calloc (strlen (path) + 9, sizeof (char));

  sprintf (file, "%s.ppp.msh", path);

  MshExportMSH (file);

  if (edgeCut == 0)
    PetscPrintf (PETSC_COMM_WORLD, "\nWarning: No cuts were made\n");

  for (j = 0; j < nbpartitions; j++)
    {
      sprintf (file, "%s.%03d.msh", path, j);

      MshExportDecomposedMSH (file, j, nbpartitions);
    }

  free (file);

  free (xadj);
  free (adjncy);
  free (weights);
  free (decomp);

}
