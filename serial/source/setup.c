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

#include "globals.h"

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"

#include "variables.h"
#include "gradient.h"
#include "geocalc.h"
#include "parser.h"

#include "setup.h"

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

      V_SetCmp (&xu, element + 1, 0.0);
      V_SetCmp (&xv, element + 1, 0.0);
      V_SetCmp (&xw, element + 1, 0.0);
      V_SetCmp (&xp, element + 1, 0.0);
      V_SetCmp (&xT, element + 1, 0.0);
      V_SetCmp (&xs, element + 1, 0.0);

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
	      V_SetCmp (&xu, element + 1, value);

	      strcpy (gs, bcvolumes[volume].fv);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xv, element + 1, value);

	      strcpy (gs, bcvolumes[volume].fw);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xw, element + 1, value);

	      strcpy (gs, bcvolumes[volume].fp);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xp, element + 1, value);

	      strcpy (gs, bcvolumes[volume].fT);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xT, element + 1, value);

	      strcpy (gs, bcvolumes[volume].fs);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xs, element + 1, value);

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

  Asgn_VV (&xu0, &xu);
  Asgn_VV (&xv0, &xv);
  Asgn_VV (&xw0, &xw);
  Asgn_VV (&xp0, &xp);
  Asgn_VV (&xT0, &xT);
  Asgn_VV (&xs0, &xs);

}

void
SetInitialFlux ()
{

  int i;

  int face, pair;
  int element, neighbor;

  //double dNf, dPf;
  double lambda;

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

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

	  // Element face velocity 
	  V_SetCmp (&uf, face + 1,
		    (V_GetCmp (&xu, neighbor + 1) * lambda + V_GetCmp (&xu, element + 1) * (1.0 - lambda)) * faces[face].n.x +
		    (V_GetCmp (&xv, neighbor + 1) * lambda + V_GetCmp (&xv, element + 1) * (1.0 - lambda)) * faces[face].n.y +
		    (V_GetCmp (&xw, neighbor + 1) * lambda + V_GetCmp (&xw, element + 1) * (1.0 - lambda)) * faces[face].n.z);

	}
      else
	{

	  // Element face velocity 
	  V_SetCmp (&uf, face + 1,
		    V_GetCmp (&xu, element + 1) * faces[face].n.x +
		    V_GetCmp (&xv, element + 1) * faces[face].n.y +
		    V_GetCmp (&xw, element + 1) * faces[face].n.z);


	}
    }

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

      if (pair != -1)
	{
	  faces[face].bc = NONE;
	}

      V_SetCmp (&xuf, face + 1, 0.0);
      V_SetCmp (&xvf, face + 1, 0.0);
      V_SetCmp (&xwf, face + 1, 0.0);
      V_SetCmp (&xpf, face + 1, 0.0);
      V_SetCmp (&xsf, face + 1, 0.0);
      V_SetCmp (&xTf, face + 1, 0.0);

    }

  for (j = 0; j < nbbcsurfaces; j++)
    {

      surface = j;

      for (i = 0; i < nbfaces; i++)
	{

	  face = i;

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
	      V_SetCmp (&xuf, face + 1, value);

	      strcpy (gs, bcsurfaces[surface].fv);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xvf, face + 1, value);

	      strcpy (gs, bcsurfaces[surface].fw);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xwf, face + 1, value);

	      strcpy (gs, bcsurfaces[surface].fp);
	      rv = evaluate (gs, &value);
	      strcat (gs, "\n");
	      V_SetCmp (&xpf, face + 1, value);

	      strcpy (gs, bcsurfaces[surface].fT);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xTf, face + 1, value);

	      strcpy (gs, bcsurfaces[surface].fs);
	      strcat (gs, "\n");
	      rv = evaluate (gs, &value);
	      V_SetCmp (&xsf, face + 1, value);

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

      //printf("element: %d, neighbor: %d\n", faces[cyclic[0]].element, faces[cyclic[1]].element);
    }

}

void
SetMaterialProperties ()
{

  int i;

  int element;

  double fr[2];

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      fr[1] = LMIN(LMAX(V_GetCmp (&xs, element + 1) * 0.5 + V_GetCmp (&xs0, element + 1) * 0.5, 0.0), 1.0);
      fr[0] = 1.0 - fr[1];

      V_SetCmp (&dens, element + 1,
		material.dens[0].cdens * fr[0] +
		material.dens[1].cdens * fr[1]);

      V_SetCmp (&visc, element + 1,
		material.visc[0].cvisc * fr[0] +
		material.visc[1].cvisc * fr[1]);

      V_SetCmp (&spheat, element + 1,
		material.therm[0].cspheat * fr[0] +
		material.therm[1].cspheat * fr[1]);

      V_SetCmp (&thcond, element + 1,
		material.therm[0].cthcond * fr[0] +
		material.therm[1].cthcond * fr[1]);

    }

}
