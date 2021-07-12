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
#include "reorder.h"
#include "rcm.h"

void
ReorderMesh (char *path)
{

  int i, j, n;

  int element, neighbor, face, pair;

  int freeAdj;

  signed char *mask;
  int *deg;

  int *xadj;
  int *adjncy;

  int *perm;

  msh_element *newelements;

  char *file;

  printf ("\nReordering mesh ...\n");

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

  perm = calloc (nbelements, sizeof (int));
  mask = calloc (nbelements, sizeof (signed char));
  deg = calloc (nbelements, sizeof (int));

  genrcmi (nbelements, 0, xadj, adjncy, perm, mask, deg);

  file = calloc (strlen (path) + 9, sizeof (char));

  // Create indexes

  newelements = calloc (nbelements, sizeof (msh_element));

  for (i = 0; i < nbelements; i++)
    {
      element = i;

      newelements[element] = elements[perm[element]];

      newelements[element].index = i;

    }

  for (i = 0; i < nbelements; i++)
    {
      element = i;

      elements[element] = newelements[element];
    }

  // Export mesh file

  sprintf (file, "%s.msh", path);

  MshExportMSH (file);

  free (file);

  free (newelements);

  free (xadj);
  free (adjncy);
  free (perm);
  free (mask);
  free (deg);

}
