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

#include <string.h>
#include <malloc.h>
#include <math.h>

#include "variables.h"
#include "vector.h"

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"

#include "globals.h"
#include "geocalc.h"
#include "parser.h"

void
SetInitialConditions ()
{

  int i, j;

  int element, volume;

  int rv;
  double value;

  gs = calloc (MAXL, sizeof (char));

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      elements[element].bc = NONE;

      V_SetCmp (&xu, elements[element].index, 0.0);
      V_SetCmp (&xv, elements[element].index, 0.0);
      V_SetCmp (&xw, elements[element].index, 0.0);
      V_SetCmp (&xp, elements[element].index, 0.0);
      V_SetCmp (&xT, elements[element].index, 0.0);
      V_SetCmp (&xs, elements[element].index, 0.0);

    }

  for (j = 0; j < nbbcvolumes; j++)
    {

      volume = j;

      for (i = 0; i < nbelements; i++)
	{

	  element = i;

	  if (elements[element].physreg == bcvolumes[volume].physreg)
	    {

	      elements[element].bc = bcvolumes[volume].bc;

	      cx = elements[element].celement.x;
	      cy = elements[element].celement.y;
	      cz = elements[element].celement.z;

	      strcpy (gs, bcvolumes[volume].fu);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xu, elements[element].index, value);

	      strcpy (gs, bcvolumes[volume].fv);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xv, elements[element].index, value);

	      strcpy (gs, bcvolumes[volume].fw);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xw, elements[element].index, value);

	      strcpy (gs, bcvolumes[volume].fp);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xp, elements[element].index, value);

	      strcpy (gs, bcvolumes[volume].fT);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xT, elements[element].index, value);

	      strcpy (gs, bcvolumes[volume].fs);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xs, elements[element].index, value);

	    }

	}

    }

  free (gs);

  for (j = 0; j < nbbcvolumes; j++)
    {

      volume = j;

      free (bcvolumes[volume].fu);
      free (bcvolumes[volume].fv);
      free (bcvolumes[volume].fw);
      free (bcvolumes[volume].fp);
      free (bcvolumes[volume].fT);
      free (bcvolumes[volume].fs);

    }

  free (bcvolumes);

  VecAssemblyBegin (xu);
  VecAssemblyEnd (xu);
  VecAssemblyBegin (xv);
  VecAssemblyEnd (xv);
  VecAssemblyBegin (xw);
  VecAssemblyEnd (xw);
  VecAssemblyBegin (xp);
  VecAssemblyEnd (xp);
  VecAssemblyBegin (xT);
  VecAssemblyEnd (xT);
  VecAssemblyBegin (xs);
  VecAssemblyEnd (xs);

  VecGhostUpdateBegin (xu, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (xu, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateBegin (xv, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (xv, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateBegin (xw, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (xw, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateBegin (xp, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (xp, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateBegin (xT, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (xT, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateBegin (xs, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (xs, INSERT_VALUES, SCATTER_FORWARD);

  VecCopy (xu, xu0);
  VecCopy (xv, xv0);
  VecCopy (xw, xw0);
  VecCopy (xp, xp0);
  VecCopy (xT, xT0);
  VecCopy (xs, xs0);

}

void
SetInitialFlux ()
{

  int i;

  int face, pair;
  int element, neighbor;

  msh_element ghost;

  //double dNf, dPf;
  double lambda;

  VecGhostGetLocalForm (xu, &xul);
  VecGhostGetLocalForm (xv, &xvl);
  VecGhostGetLocalForm (xw, &xwl);

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

      // Internal face - neighbor in this processor
      if (pair != -1)
	{

	  neighbor = faces[pair].element;

          /*
	  dNf =
	    GeoMagVector (GeoSubVectorVector
			  (elements[neighbor].celement, faces[face].cface));
	  dPf =
	    GeoMagVector (GeoSubVectorVector
			  (elements[element].celement, faces[face].cface));

	  lambda = dPf / (dPf + dNf);
	  */

	  lambda = 0.5;

	  V_SetCmp (&uf, face,
		    (V_GetCmp (&xul, neighbor) * lambda +
		     V_GetCmp (&xul, element) * (1 - lambda)) *
		    faces[face].n.x +
		    (V_GetCmp (&xvl, neighbor) * lambda +
		     V_GetCmp (&xvl, element) * (1 - lambda)) *
		    faces[face].n.y +
		    (V_GetCmp (&xwl, neighbor) * lambda +
		     V_GetCmp (&xwl, element) * (1 - lambda)) *
		    faces[face].n.z);

	}
      else
	{

	  if (faces[face].bc == PROCESSOR)
	    {

	      // Face between processors
	      ghost.index = faces[face].physreg;

	      ghost.celement.x = V_GetCmp (&cexl, faces[face].ghost);
	      ghost.celement.y = V_GetCmp (&ceyl, faces[face].ghost);
	      ghost.celement.z = V_GetCmp (&cezl, faces[face].ghost);

	      /*
	      dNf =
		GeoMagVector (GeoSubVectorVector
			      (ghost.celement, faces[face].cface));
	      dPf =
		GeoMagVector (GeoSubVectorVector
			      (elements[element].celement,
			       faces[face].cface));

	      lambda = dPf / (dPf + dNf);
	      */

	      lambda = 0.5;

	      V_SetCmp (&uf, face,
			(V_GetCmp (&xul, faces[face].ghost) * lambda +
			 V_GetCmp (&xul, element) * (1 - lambda)) *
			faces[face].n.x +
			(V_GetCmp (&xvl, faces[face].ghost) * lambda +
			 V_GetCmp (&xvl, element) * (1 - lambda)) *
			faces[face].n.y +
			(V_GetCmp (&xwl, faces[face].ghost) * lambda +
			 V_GetCmp (&xwl, element) * (1 - lambda)) *
			faces[face].n.z);


	    }
	  else
	    {

	      // Boundary face

	      V_SetCmp (&uf, face,
			V_GetCmp (&xul, element) *
			faces[face].n.x +
			V_GetCmp (&xvl, element) *
			faces[face].n.y +
			V_GetCmp (&xwl, element) * faces[face].n.z);

	    }

	}

    }

  VecGhostRestoreLocalForm (xu, &xul);
  VecGhostRestoreLocalForm (xv, &xvl);
  VecGhostRestoreLocalForm (xw, &xwl);

  VecAssemblyBegin (uf);
  VecAssemblyEnd (uf);

}

void
SetBoundary ()
{

  int i, j, n;

  int face, pair, surface;

  int rv;
  double value;

  int cyclic[2];

  gs = calloc (MAXL, sizeof (char));

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      pair = faces[face].pair;

      if (pair != -1 && faces[face].bc != PROCESSOR)
	{
	  faces[face].bc = NONE;
	}

      V_SetCmp (&xuf, face, 0.0);
      V_SetCmp (&xvf, face, 0.0);
      V_SetCmp (&xwf, face, 0.0);
      V_SetCmp (&xpf, face, 0.0);
      V_SetCmp (&xTf, face, 0.0);
      V_SetCmp (&xsf, face, 0.0);

    }

  for (j = 0; j < nbbcsurfaces; j++)
    {

      surface = j;

      for (i = 0; i < nbfaces; i++)
	{

	  face = i;

	  if (faces[face].bc == PROCESSOR)
	    continue;

	  pair = faces[face].pair;

	  if (pair != -1)
	    continue;

	  if (faces[face].physreg == bcsurfaces[surface].physreg)
	    {

	      faces[face].bc = bcsurfaces[surface].bc;

	      cx = faces[face].cface.x;
	      cy = faces[face].cface.y;
	      cz = faces[face].cface.z;

	      strcpy (gs, bcsurfaces[surface].fu);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xuf, face, value);

	      strcpy (gs, bcsurfaces[surface].fv);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xvf, face, value);

	      strcpy (gs, bcsurfaces[surface].fw);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xwf, face, value);

	      strcpy (gs, bcsurfaces[surface].fp);
	      rv = evaluate (gs, &value);
	      strcat (gs, "\n");
	      V_SetCmp (&xpf, face, value);

	      strcpy (gs, bcsurfaces[surface].fT);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xTf, face, value);

	      strcpy (gs, bcsurfaces[surface].fs);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xsf, face, value);

	    }

	}

    }

  free (gs);

  for (j = 0; j < nbbcsurfaces; j++)
    {

      surface = j;

      free (bcsurfaces[surface].fu);
      free (bcsurfaces[surface].fv);
      free (bcsurfaces[surface].fw);
      free (bcsurfaces[surface].fp);
      free (bcsurfaces[surface].fT);
      free (bcsurfaces[surface].fs);

    }

  free (bcsurfaces);

  n = 0;

  for (i = 0; i < nbfaces; i++)
    {
      face = i;

      if (faces[face].bc == CYCLIC)
	{
	  if (n < 2)
	    {
	      cyclic[n] = face;
	      n++;

	    }
	}
    }

  if (n == 2)
    {
      faces[cyclic[0]].pair = cyclic[1];
      faces[cyclic[1]].pair = cyclic[0];

      faces[cyclic[0]].dj *= 2.0;
      faces[cyclic[1]].dj *= 2.0;

      //PetscPrintf(PETSC_COMM_WORLD, "element: %d, neighbor: %d\n", faces[cyclic[0]].element, faces[cyclic[1]].element);

    }

  VecAssemblyBegin (xuf);
  VecAssemblyEnd (xuf);
  VecAssemblyBegin (xvf);
  VecAssemblyEnd (xvf);
  VecAssemblyBegin (xwf);
  VecAssemblyEnd (xwf);
  VecAssemblyBegin (xpf);
  VecAssemblyEnd (xpf);
  VecAssemblyBegin (xTf);
  VecAssemblyEnd (xTf);
  VecAssemblyBegin (xsf);
  VecAssemblyEnd (xsf);

}

void
SetMaterialProperties ()
{

  int i;

  int element;

  double fr[2];

  VecGhostGetLocalForm (xs0, &xs0l);
  VecGhostGetLocalForm (xs, &xsl);

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      fr[1] = LMIN(LMAX(V_GetCmp (&xsl, element) * 0.5 + V_GetCmp (&xs0l, element) * 0.5, 0.0), 1.0);
      fr[0] = 1.0 - fr[1];

      V_SetCmp (&dens, elements[element].index,
		material.dens[0].cdens * fr[0] +
		material.dens[1].cdens * fr[1]);
      V_SetCmp (&visc, elements[element].index,
		material.visc[0].cvisc * fr[0] +
		material.visc[1].cvisc * fr[1]);
      V_SetCmp (&spheat, elements[element].index,
		material.therm[0].cspheat * fr[0] +
		material.therm[1].cspheat * fr[1]);
      V_SetCmp (&thcond, elements[element].index,
		material.therm[0].cthcond * fr[0] +
		material.therm[1].cthcond * fr[1]);

    }

  VecGhostRestoreLocalForm (xs, &xsl);
  VecGhostRestoreLocalForm (xs0, &xs0l);

  VecAssemblyBegin (dens);
  VecAssemblyEnd (dens);
  VecAssemblyBegin (visc);
  VecAssemblyEnd (visc);
  VecAssemblyBegin (thcond);
  VecAssemblyEnd (thcond);
  VecAssemblyBegin (spheat);
  VecAssemblyEnd (spheat);

  VecGhostUpdateBegin (dens, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (dens, INSERT_VALUES, SCATTER_FORWARD);

  VecGhostUpdateBegin (visc, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (visc, INSERT_VALUES, SCATTER_FORWARD);

  VecGhostUpdateBegin (spheat, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (spheat, INSERT_VALUES, SCATTER_FORWARD);

  VecGhostUpdateBegin (thcond, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (thcond, INSERT_VALUES, SCATTER_FORWARD);

}

void
SetCenters ()
{

  int i;

  int element;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      V_SetCmp (&cex, elements[element].index, elements[element].celement.x);
      V_SetCmp (&cey, elements[element].index, elements[element].celement.y);
      V_SetCmp (&cez, elements[element].index, elements[element].celement.z);

    }

  VecAssemblyBegin (cex);
  VecAssemblyEnd (cex);
  VecAssemblyBegin (cey);
  VecAssemblyEnd (cey);
  VecAssemblyBegin (cez);
  VecAssemblyEnd (cez);

  VecGhostUpdateBegin (cex, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (cex, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateBegin (cey, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (cey, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateBegin (cez, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (cez, INSERT_VALUES, SCATTER_FORWARD);

}

void
SetGhosts ()
{

  int i;

  int face;

  nbghosts = 0;

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      if (faces[face].pair != -1)
	continue;

      if (faces[face].bc == PROCESSOR)
	{
	  faces[face].ghost = nbelements + nbghosts;
	  nbghosts++;
	}

    }

  ghosts = calloc (nbghosts, sizeof (int));

  nbghosts = 0;

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      if (faces[face].pair != -1)
	continue;

      if (faces[face].bc == PROCESSOR)
	{

	  ghosts[nbghosts] = faces[face].physreg;

	  nbghosts++;

	}

    }

}
