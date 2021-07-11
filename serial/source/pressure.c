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

#include "globals.h"

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"

#include "variables.h"
#include "gradient.h"
#include "geocalc.h"
#include "setup.h"
#include "msolver.h"

#include "pressure.h"

void
CorrectFaceP ()
{

  unsigned int i, k;

  unsigned int node;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  //double dNf, dPf;
  double lambda;

  double ppl;

  double phij;

  msh_vector gradpp;

  double apj;

  double ghf;

  msh_vector g;

  g.x = parameter.g[0];
  g.y = parameter.g[1];
  g.z = parameter.g[2];

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

      if (parameter.orthof != 0.0)
	gradpp = Gradient (&xp, &xpf, LOGICAL_TRUE, element);

      if (pair != -1)
	{

	  neighbor = faces[pair].element;

	  // Cell-based linear interpolation

	  /*
	  dNf = GeoMagVector (GeoSubVectorVector(elements[neighbor].celement, faces[face].cface));
	  dPf = GeoMagVector (GeoSubVectorVector(elements[element].celement, faces[face].cface));

	  lambda = dPf / (dPf + dNf);
	  */

	  lambda = 0.5;

	  phij = V_GetCmp (&xp, neighbor + 1) * lambda +
	    V_GetCmp (&xp, element + 1) * (1.0 - lambda);

	  V_SetCmp (&xpf, face + 1, phij);

	}
      else
	{

	  ppl = V_GetCmp (&xp, element + 1);

	  apj = V_GetCmp (&ap, element + 1);

	  if (parameter.orthof != 0.0)
	    {
	      ppl += parameter.orthof * GeoDotVectorVector (gradpp,
							    GeoSubVectorVector
							    (faces[face].rpl,
							     elements
							     [element].
							     celement));
	    }

	  ghf =
	    V_GetCmp (&dens, element + 1) * GeoDotVectorVector (g, faces[face].d);

	  ppl += ghf;

	  if (faces[face].bc == PERMEABLE)
	    {

	      V_SetCmp (&xpf, face + 1, V_GetCmp (&xpf, face + 1) * (1.0 - V_GetCmp (&xs, element + 1)) + ppl * V_GetCmp (&xs, element + 1));

	    }

	  if (faces[face].bc == INLET)
	    {

	      V_SetCmp (&xpf, face + 1, ppl - V_GetCmp (&uf, face + 1) * apj * (faces[face].dj + faces[face].kj));

	    }

	  if (faces[face].bc == WALL ||
	      faces[face].bc == MOVINGWALL ||
	      faces[face].bc == ADIABATICWALL ||
	      faces[face].bc == SLIP || faces[face].bc == SURFACE)
	    {

	      V_SetCmp (&xpf, face + 1, ppl);

	    }

	}

    }

}

void
BuildContinuityMatrix (double dt)
{

  unsigned int i, j, n, nj;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  double acp;
  double acn[MAXFACES];
  unsigned int ani[MAXFACES];
  double bcp;

  double apj;

  double Huj, Hvj, Hwj;
  double Hf;

  //double dNf, dPf;
  double lambda;

  msh_vector gradpp;
  msh_vector gradpn;

  // Equation: div(U) = 0

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      acp = 0.0;

      bcp = 0.0;

      n = 0;

      if (parameter.orthof != 0.0)
	gradpp = Gradient (&xp, &xpf, LOGICAL_TRUE, element);

      for (j = 0; j < elements[element].nbfaces; j++)
	{

	  face = elements[element].face[j];

	  pair = faces[face].pair;

	  if (pair != -1)
	    {

	      neighbor = faces[pair].element;

              /*
	      dNf = GeoMagVector (GeoSubVectorVector (elements[neighbor].celement, faces[face].cface));
	      dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement,  faces[face].cface));

	      lambda = dPf / (dPf + dNf);
	      */

	      lambda = 0.5;

	      apj = V_GetCmp (&ap, neighbor + 1) * lambda + V_GetCmp (&ap, element + 1) * (1.0 - lambda);

	      Huj = V_GetCmp (&hu, neighbor + 1) * lambda + V_GetCmp (&hu, element + 1) * (1.0 - lambda);
	      Hvj = V_GetCmp (&hv, neighbor + 1) * lambda + V_GetCmp (&hv, element + 1) * (1.0 - lambda);
	      Hwj = V_GetCmp (&hw, neighbor + 1) * lambda + V_GetCmp (&hw, element + 1) * (1.0 - lambda);

	      Hf = Huj * faces[face].n.x + Hvj * faces[face].n.y + Hwj * faces[face].n.z;

	      acp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj;

	      acn[n] = 1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj;

	      ani[n] = neighbor;
	      n++;

	      V_SetCmp (&uf, face + 1, Hf / apj);

	      bcp += V_GetCmp (&uf, face + 1) * faces[face].Aj;

	      // Non-orthogonal correction term  
	      if (parameter.orthof != 0.0)
		gradpn = Gradient (&xp, &xpf, LOGICAL_TRUE, neighbor);

	      if (parameter.orthof != 0.0)
		{
		  bcp += -1.0 * parameter.orthof / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj *
		    (GeoDotVectorVector (gradpn, GeoSubVectorVector (faces[face].rnl, elements[neighbor].celement)) -
		     GeoDotVectorVector (gradpp, GeoSubVectorVector (faces[face].rpl, elements[element].celement)));
		}

	    }
	  else
	    {

	      apj = V_GetCmp (&ap, element + 1);

	      Huj = V_GetCmp (&hu, element + 1);
	      Hvj = V_GetCmp (&hv, element + 1);
	      Hwj = V_GetCmp (&hw, element + 1);

	      Hf = Huj * faces[face].n.x + Hvj * faces[face].n.y + Hwj * faces[face].n.z;

	      if (faces[face].bc == PERMEABLE)
		{

		  acp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj * (1.0 - V_GetCmp (&xs, element + 1));
		  bcp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj * V_GetCmp (&xpf, face + 1) * (1.0 - V_GetCmp (&xs, element + 1));

		  V_SetCmp (&uf, face + 1, Hf / apj * (1.0 - V_GetCmp (&xs, element + 1)));

		  bcp += V_GetCmp (&uf, face + 1) * faces[face].Aj;

		  // Non-orthogonal correction term
		  if (parameter.orthof != 0.0)
		    {
		      bcp += 1.0 * parameter.orthof / (apj * (faces[face].dj + faces[face].kj)) * 
			(GeoDotVectorVector (gradpp,  GeoSubVectorVector (faces[face].rpl, elements[element].celement))) * 
			faces[face].Aj * (1.0 - V_GetCmp (&xs, element + 1));
		    }

		}

	      if (faces[face].bc == OUTLET)
		{

		  // velocity gradient = 0
		  // specified pressure         

		  acp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj;
		  bcp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj * V_GetCmp (&xpf, face + 1);

		  V_SetCmp (&uf, face + 1, Hf / apj);

		  bcp += V_GetCmp (&uf, face + 1) * faces[face].Aj;

		  // Non-orthogonal correction term     
		  if (parameter.orthof != 0.0)
		    {
		      bcp += 1.0 * parameter.orthof / (apj * (faces[face].dj + faces[face].kj)) *
			(GeoDotVectorVector (gradpp, GeoSubVectorVector (faces[face].rpl, elements[element].celement))) * faces[face].Aj;
		    }

		}

	      if (faces[face].bc == PRESSUREINLET)
		{

		  // specified pressure 
		  // velocity gradient = 0

		  acp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj;
		  bcp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj * V_GetCmp (&xpf, face + 1);

		  V_SetCmp (&uf, face + 1, Hf / apj);

		  bcp += V_GetCmp (&uf, face + 1) * faces[face].Aj;

		  // Non-orthogonal correction term
		  if (parameter.orthof != 0.0)
		    {
		      bcp += 1.0 * parameter.orthof / (apj * (faces[face].dj + faces[face].kj)) *
			(GeoDotVectorVector (gradpp, GeoSubVectorVector (faces[face].rpl, elements[element].celement))) * faces[face].Aj;
		    }
		}

	      if (faces[face].bc == INLET ||
		  faces[face].bc == MOVINGWALL ||
		  faces[face].bc == WALL ||
		  faces[face].bc == ADIABATICWALL ||
		  faces[face].bc == SURFACE)
		{

		  // pressure gradient = 0
		  // specified velocity

		  V_SetCmp (&uf, face + 1,
			    V_GetCmp (&xuf, face + 1) * faces[face].n.x +
			    V_GetCmp (&xvf, face + 1) * faces[face].n.y +
			    V_GetCmp (&xwf, face + 1) * faces[face].n.z);

		  bcp += V_GetCmp (&uf, face + 1) * faces[face].Aj;

		}

	    }

	}

      if (acp == 0.0 || acp != acp)
	{
	  printf ("\nError: Problem setting up continuity matrix\n");
	  exit (LOGICAL_ERROR);
	}

      nj = 0;

      for (j = 0; j < n; j++)
	{
	  if (ani[j] > element)
	    nj++;
	}

      Q_SetLen (&Ac, element + 1, nj + 1);

      Q_SetEntry (&Ac, element + 1, 0, element + 1, acp);

      nj = 0;

      for (j = 0; j < n; j++)
	{
	  if (ani[j] > element)
	    {
	      Q_SetEntry (&Ac, element + 1, nj + 1, ani[j] + 1, acn[j]);
	      nj++;
	    }
	}

      V_SetCmp (&bp, element + 1, bcp);

    }

}

void
CalculatePressure (char *var, int *fiter, double dt, double maxCp,
		   int verbose, int pchecks)
{

  unsigned int i, j;

  double mres;
  int miter;
  double mtime;

  double presc;

  if (parameter.calc[ip] == LOGICAL_FALSE)
    return;

  V_Constr (&xpp, "Pressure at cell center - previous iteration", nbelements, Normal, True);

  // Store previous time step values
  Asgn_VV (&xp0, &xp);

  fiter[ip]++;

  for (i = 0; i <= parameter.northocor; i++)
    {

      Q_Constr (&Ac, "Continuity matrix", nbelements, True, Rowws, Normal, True);
      V_Constr (&bp, "Continuity source", nbelements, Normal, True);

      // Store previous iteration values
      if (parameter.northocor > 0)
      {
        Asgn_VV (&xpp, &xp);
      }

      // Build the continuity matrix (mass conservation)
      BuildContinuityMatrix (dt);

      if (pchecks == LOGICAL_TRUE)
	{
	  if (!CheckIfDiagonalMatrix (&Ac))
	    {
	      printf
		("\nWarning: Continuity matrix is not diagonal dominant\n");
	      WriteMatrix (&Ac, LOGICAL_TRUE);
	      WriteVector (&bp);
	      //exit (LOGICAL_ERROR);
	    }
	}

      // Set matrix solution accuracy
      SetRTCAccuracy (parameter.mtol[ip]);

      // Solve matrix to get pressure p
      SolveMatrix (&Ac, &xp, &bp, &miter, &mres, &mtime,
		   parameter.msolver[ip], parameter.mprecond[ip],
		   parameter.miter[ip]);

      if (verbose == LOGICAL_TRUE)
	printf
	  ("\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
	   var[ip], miter, mres, mtime);

      if ((mres > parameter.mtol[ip] && miter == parameter.miter[ip])
	  || LASResult () != LASOK)
	{
	  printf ("\nError: Problem solving matrix %c\n", var[ip]);
	  exit (LOGICAL_ERROR);
	}

      presc = 0.0;

      // Calculate pressure convergence
      if (parameter.northocor > 0)
      {
        presc = l2Norm_V (Sub_VV (&xpp, &xp));
      }

      if (verbose == LOGICAL_TRUE)
	printf ("\nNon-orthogonality error %d (continuity): %+E\n", i, presc);

      CorrectFaceP ();

      Q_Destr (&Ac);
      V_Destr (&bp);

      if (presc < parameter.mtol[ip])
	break;

    }

  V_Destr (&xpp);

}
