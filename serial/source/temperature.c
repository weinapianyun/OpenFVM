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

#include "temperature.h"

void
CorrectFaceT ()
{

  unsigned int i, j;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  //double dNf, dPf;
  double lambda;

  double Tpl;

  msh_vector gradTp;

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

      if (pair != -1)
	{

	  neighbor = faces[pair].element;

	  /*
	  dNf = GeoMagVector (GeoSubVectorVector (elements[neighbor].celement, faces[face].cface));
	  dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement, faces[face].cface));

	  lambda = dPf / (dPf + dNf);
	  */
	
          lambda = 0.5;

	  V_SetCmp (&xTf, face + 1, V_GetCmp (&xT, neighbor + 1) * lambda + V_GetCmp (&xT, element + 1) * (1.0 - lambda));

	}
      else
	{

	  if (faces[face].bc == ADIABATICWALL || faces[face].bc == OUTLET || faces[face].bc == EMPTY)
	    {

	      // zero gradient

	      Tpl = V_GetCmp (&xT, element + 1);

	      if (parameter.orthof != 0.0)
		{

		  gradTp = Gradient (&xT, &xTf, LOGICAL_TRUE, element);

		  Tpl += parameter.orthof * GeoDotVectorVector (gradTp, GeoSubVectorVector (faces[face].rpl, elements[element].celement));
		}

	      V_SetCmp (&xTf, face + 1, Tpl);

	    }

	}

    }

}

void
BuildEnergyMatrix (double dt, double schemefactor)
{

  int i, j, n;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  double aep;
  double aen[MAXFACES];
  unsigned int ani[MAXFACES];

  double bep;

  //double dNf, dPf;
  double lambda;
  double xsi;

  msh_vector gradTp;
  msh_vector gradTn;

  double densp;
  double spheatp;
  double thcondj;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      aep = 0.0;
      bep = 0.0;

      n = 0;

      if (parameter.orthof != 0.0)
	gradTp = Gradient (&xT, &xTf, LOGICAL_TRUE, element);

      densp = V_GetCmp (&dens, element + 1);
      spheatp = V_GetCmp (&spheat, element + 1);

      for (j = 0; j < elements[element].nbfaces; j++)
	{

	  face = elements[element].face[j];

	  pair = faces[face].pair;

	  if (pair != -1)
	    {

	      neighbor = faces[pair].element;

              /*
	      dNf = GeoMagVector (GeoSubVectorVector (elements[neighbor].celement, faces[face].cface));
	      dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement, faces[face].cface));

	      lambda = dPf / (dPf + dNf);
	      */
	
	      lambda = 0.5;

	      thcondj = V_GetCmp (&thcond, element + 1) * (1.0 - lambda) + V_GetCmp (&thcond, neighbor + 1) * lambda;

	      // Conduction 
	      aep += schemefactor * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;
	      aen[n] = -schemefactor * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;

	      bep += -(1.0 - schemefactor) * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp * V_GetCmp (&xT0, element + 1);
	      bep += +(1.0 - schemefactor) * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp * V_GetCmp (&xT0, neighbor + 1);

	      // Convection 
	      if (parameter.scheme[iT] == UDS)
		{
		  // UDS
		  if (V_GetCmp (&uf, face + 1) > 0.0)
		    xsi = 0.0;
		  else
		    xsi = 1.0;
		}
	      else
		{
		  //CDS
		  xsi = lambda;
		}

	      // Convection 
	      aep += schemefactor * (1.0 - xsi) * densp * spheatp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp;
	      aen[n] += schemefactor * xsi * densp * spheatp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp;

	      ani[n] = neighbor;
	      n++;

	      bep += -(1.0 - schemefactor) * (1.0 - xsi) * densp * spheatp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp * V_GetCmp (&xT0, element + 1);
	      bep += -(1.0 - schemefactor) * xsi * densp * spheatp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp * V_GetCmp (&xT0, neighbor + 1);

	      if (parameter.orthof != 0.0)
		{

		  gradTn = Gradient (&xT, &xTf, LOGICAL_TRUE, neighbor);

		  // Non-orthogonal correction term         
		  bep += parameter.orthof * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
		    elements[element].Vp * (GeoDotVectorVector (gradTn, GeoSubVectorVector (faces[face].rnl, elements[neighbor].celement)) - GeoDotVectorVector (gradTp, GeoSubVectorVector (faces[face].rpl, elements[element].celement)));
		}

	    }
	  else
	    {

	      thcondj = material.bthcond;

	      if (faces[face].bc != EMPTY && faces[face].bc != ADIABATICWALL)
		{

		  // Conduction
		  aep += schemefactor * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;
		  bep += schemefactor * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp * V_GetCmp (&xTf, face + 1);

		  // Convection
		  if (parameter.scheme[iT] == UDS)
		    {
		      // UDS
		      if (V_GetCmp (&uf, face + 1) > 0.0)
			{
			  aep += schemefactor * densp * spheatp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp;
			}
		      else
			{
			  bep += -schemefactor * densp * spheatp * V_GetCmp (&uf, face +  1) * faces[face].Aj / elements[element].Vp * V_GetCmp (&xTf, face + 1);
			}

		    }
		  else
		    {
		      // CDS
		      bep += -schemefactor * densp * spheatp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp * V_GetCmp (&xTf, face + 1);
		    }

		  if (parameter.orthof != 0.0)
		    {
		      // Non-orthogonal correction term
		      bep += parameter.orthof * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp * GeoDotVectorVector (gradTp, GeoSubVectorVector (faces[face].rpl, elements[element].celement));
		    }

		}

	    }

	}

      if (dt > 0)
	{

	  // Unsteady term
	  aep += densp * spheatp / dt;
	  bep += densp * spheatp / dt * V_GetCmp (&xT0, element + 1);

	}

      if (aep == 0.0 || aep != aep)
	{
	  printf ("\nError: Problem setting up energy matrix\n");
	  exit (LOGICAL_ERROR);
	}

      Q_SetLen (&Ae, element + 1, n + 1);

      Q_SetEntry (&Ae, element + 1, 0, element + 1, aep);

      for (j = 0; j < n; j++)
	{
	  Q_SetEntry (&Ae, element + 1, j + 1, ani[j] + 1, aen[j]);
	}

      V_SetCmp (&bT, element + 1, bep);

    }

}

void
SolveEnergyExplicit (double dt)
{

  unsigned int i, j, n;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  double aep;
  double aen[MAXFACES];
  unsigned int ani[MAXFACES];

  double bep;

  //double dNf, dPf;
  double lambda;
  double xsi;

  msh_vector gradTp;
  msh_vector gradTn;

  double densp;
  double spheatp;
  double thcondj;

  double sumT;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      aep = 0.0;
      bep = 0.0;

      n = 0;

      if (parameter.orthof != 0.0)
	gradTp = Gradient (&xT, &xTf, LOGICAL_TRUE, element);

      densp = V_GetCmp (&dens, element + 1);
      spheatp = V_GetCmp (&spheat, element + 1);

      for (j = 0; j < elements[element].nbfaces; j++)
	{

	  face = elements[element].face[j];

	  pair = faces[face].pair;

	  if (pair != -1)
	    {

	      neighbor = faces[pair].element;

              /*
	      dNf = GeoMagVector (GeoSubVectorVector (elements[neighbor].celement, faces[face].cface));
	      dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement, faces[face].cface));

	      lambda = dPf / (dPf + dNf);
              */

 	      lambda = 0.5;

	      thcondj = V_GetCmp (&thcond, element + 1) * (1.0 - lambda) + V_GetCmp (&thcond, neighbor + 1) * lambda;

	      // Conduction 
	      aep += thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;
	      aen[n] = -thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;

	      // Convection 
	      if (parameter.scheme[iT] == UDS)
		{
		  // UDS
		  if (V_GetCmp (&uf, face + 1) > 0.0)
		    xsi = 0.0;
		  else
		    xsi = 1.0;
		}
	      else
		{
		  //CDS
		  xsi = lambda;
		}

	      // Convection 
	      aep += (1.0 - xsi) * densp * spheatp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp;
	      aen[n] += xsi * densp * spheatp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp;

	      ani[n] = neighbor;
	      n++;

	      if (parameter.orthof != 0.0)
		{

		  gradTn = Gradient (&xT, &xTf, LOGICAL_TRUE, neighbor);

		  // Non-orthogonal correction term         
		  bep += parameter.orthof * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp *
		    (GeoDotVectorVector (gradTn, GeoSubVectorVector (faces[face].rnl, elements[neighbor].celement)) - 
		     GeoDotVectorVector (gradTp, GeoSubVectorVector (faces[face].rpl, elements[element].celement)));
		}

	    }
	  else
	    {

	      thcondj = material.bthcond;

	      if (faces[face].bc != EMPTY && faces[face].bc != ADIABATICWALL)
		{

		  // Conduction
		  aep += thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;
		  bep += thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp * V_GetCmp (&xTf, face + 1);
		  
		  // Convection
		  if (parameter.scheme[iT] == UDS)
		    {
		      // UDS
		      if (V_GetCmp (&uf, face + 1) > 0.0)
			{
			  aep += densp * spheatp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp;
			}
		      else
			{
			  bep += -densp * spheatp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp * V_GetCmp (&xTf, face + 1);
			}

		    }
		  else
		    {
		      // CDS
		      bep += -densp * spheatp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp * V_GetCmp (&xTf, face + 1);
		    }

		  if (parameter.orthof != 0.0)
		    {
		      // Non-orthogonal correction term           
		      bep += parameter.orthof * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp * GeoDotVectorVector (gradTp, GeoSubVectorVector (faces[face].rpl, elements[element].celement));
		    }

		}

	    }

	}

      if (dt > 0)
	{

	  // Unsteady term - Euler
	  aep += densp * spheatp / dt;

	  bep += densp * spheatp / dt * V_GetCmp (&xT0, element + 1);

	}

      sumT = 0.0;

      for (j = 0; j < n; j++)
	{
	  sumT += aen[j] * V_GetCmp (&xT0, ani[j] + 1);
	}

      V_SetCmp (&xT, element + 1, (bep - sumT) / aep);

    }

}

void
CalculateTemperature (char *var, int *fiter, double dt, double maxCp,
		      int verbose, int pchecks)
{

  unsigned int i;

  double mres;
  int miter;
  double mtime;

  double tempc;

  if (parameter.calc[iT] == LOGICAL_FALSE)
    return;

  V_Constr (&xTp, "Temperature at cell center - previous iteration", nbelements, Normal, True);

  // Store previous time step values 
  Asgn_VV (&xT0, &xT);

  fiter[iT]++;

  for (i = 0; i <= parameter.northocor; i++)
    {

      if (parameter.timemethod[iT] == IMPLICITEULER
	  || parameter.timemethod[iT] == CRANKNICOLSON)
	{
	  Q_Constr (&Ae, "Energy matrix", nbelements, False, Rowws, Normal, True);
	  V_Constr (&bT, "Energy source", nbelements, Normal, True);
	}

      // Store previous iteration values
      if (parameter.northocor > 0)
      {
        Asgn_VV (&xTp, &xT);
      }

      if (parameter.timemethod[iT] == EXPLICITEULER)
	{
	  SolveEnergyExplicit (dt);
	}

      if (parameter.timemethod[iT] == IMPLICITEULER)
	{
	  BuildEnergyMatrix (dt, 1.0);
	}

      if (parameter.timemethod[iT] == CRANKNICOLSON)
	{
	  BuildEnergyMatrix (dt, 0.5);
	}

      if (parameter.timemethod[iT] == IMPLICITEULER
	  || parameter.timemethod[iT] == CRANKNICOLSON)
	{

	  if (pchecks == LOGICAL_TRUE)
	    {
	      if (!CheckIfDiagonalMatrix (&Ae))
		{
		  printf
		    ("\nWarning: Energy matrix is not diagonal dominant\n");
		  WriteMatrix (&Ae, LOGICAL_FALSE);
		  WriteVector (&bT);
		  //exit (LOGICAL_ERROR);
		}
	    }

	  // Set matrix solution accuracy
	  SetRTCAccuracy (parameter.mtol[iT]);

	  // Solve matrix to get temperature T
	  SolveMatrix (&Ae, &xT, &bT, &miter, &mres, &mtime,
		       parameter.msolver[iT], parameter.mprecond[iT],
		       parameter.miter[iT]);

	  if (verbose == LOGICAL_TRUE)
	    printf
	      ("\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
	       var[iT], miter, mres);

	  if ((mres > parameter.mtol[iT] && miter == parameter.miter[iT])
	      || LASResult () != LASOK)
	    {
	      printf ("\nError: Problem solving matrix %c\n", var[iT]);
	      exit (LOGICAL_ERROR);
	    }

	}

      tempc = 0.0;

      // Calculate temperature convergence
      if (parameter.northocor > 0)
      {
        tempc = l2Norm_V (Sub_VV (&xTp, &xT));
      }

      if (verbose == LOGICAL_TRUE)
	printf ("\nNon-orthogonality error %d (energy): %+E\n", i, tempc);

      CorrectFaceT (&xT, &xTf);

      if (parameter.timemethod[iT] == IMPLICITEULER
	  || parameter.timemethod[iT] == CRANKNICOLSON)
	{
	  Q_Destr (&Ae);
	  V_Destr (&bT);
	}

      if (tempc < parameter.mtol[iT])
	break;

    }

  V_Destr (&xTp);

}
