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

#include "variables.h"
#include "vector.h"
#include "matrix.h"
#include "itersolv.h"

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"

#include "parallel.h"
#include "gradient.h"
#include "geocalc.h"
#include "globals.h"
#include "setup.h"
#include "msolver.h"

#include "gamma.h"

double
CalculateMaxCourantNumber (double dt, int interface) //计算最大库朗数
{

  unsigned int i, j;

  unsigned int element, face;

  double Cpp;
  double Cj;

  double s, cs;
  
  double maxCp;

  maxCp = 0.0;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      if (interface == 1)
      {     
      	s = LMIN(LMAX(V_GetCmp (&xsl, element), 0.0), 1.0);
      	cs = (1.0 - s) * (1.0 - s) * s * s * 16.0;    
      }
      else
      {
        cs = 1.0;
      }
      
      Cpp = 0.0;
       
      for (j = 0; j < elements[element].nbfaces; j++)
	{

	  face = elements[element].face[j];
   
          Cj = LMAX (-V_GetCmp (&uf, face) * faces[face].Aj * dt / elements[element].Vp, 0.0);
	  Cpp += cs * Cj;
	      
	}

      maxCp = LMAX (maxCp, Cpp);

      V_SetCmp (&Co, elements[element].index, Cpp);

    }
    
  VecAssemblyBegin (Co);
  VecAssemblyEnd (Co);

  VecGhostUpdateBegin (Co, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (Co, INSERT_VALUES, SCATTER_FORWARD);

  return maxCp;

}

void
CorrectFaceS () 
{

  int i;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  msh_element ghost;

  double betaj;

  VecGhostGetLocalForm (xs, &xsl);
  
  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

      if (pair != -1)
	{

	  neighbor = faces[pair].element;

	  if (V_GetCmp (&uf, face) > 0.0)
	    betaj = V_GetCmp (&betaf, face);
	  else
	    betaj = 1.0 - V_GetCmp (&betaf, face);

	  V_SetCmp (&xsf, face, LMAX (LMIN((1.0 - betaj) * V_GetCmp (&xsl, element) + betaj * V_GetCmp (&xsl, neighbor), 1.0), 0.0));

	}
      else
	{

	  if (faces[face].bc == PROCESSOR)
	    {

	      ghost.index = faces[face].physreg;

	      ghost.celement.x = V_GetCmp (&cexl, faces[face].ghost);
	      ghost.celement.y = V_GetCmp (&ceyl, faces[face].ghost);
	      ghost.celement.z = V_GetCmp (&cezl, faces[face].ghost);

	      if (V_GetCmp (&uf, face) > 0.0)
		betaj = V_GetCmp (&betaf, face);
	      else
		betaj = 1.0 - V_GetCmp (&betaf, face);

	      V_SetCmp (&xsf, face, LMAX (LMIN((1.0 - betaj) * V_GetCmp (&xsl, element) + betaj * V_GetCmp (&xsl, faces[face].ghost), 1.0), 0.0));

	    }
	  else
	    {

	      if (faces[face].bc == OUTLET)
		{

		  // zero gradient

		  V_SetCmp (&xsf, face, V_GetCmp (&xsl, element));

		}

	    }

	}

    }

  VecGhostRestoreLocalForm (xs, &xsl);
    
  VecAssemblyBegin (xsf);
  VecAssemblyEnd (xsf);

}

void
PredictBeta ()
{

  int i;

  register unsigned int face, pair;
  register unsigned int element, neighbor;
  register unsigned int donor, acceptor;

  double dot, l1, l2;

  double su, sdn, sjnCBC, sjnUQ, sjn;

  double qj, tetaj;

  double betaj;

  double Cod;

  double ang;

  msh_vector grads;

  VecGhostGetLocalForm (xs, &xsl);
  VecGhostGetLocalForm (Co, &Col);
  
  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

      betaj = 0.0;

      if (parameter.ncicsamcor != 0)
	{

	  if (pair != -1)
	    {

	      neighbor = faces[pair].element;

	      if (V_GetCmp (&uf, face) != 0.0)
		{

		  if (V_GetCmp (&uf, face) > 0.0)
		    {

		      acceptor = neighbor;
		      donor = element;

		      grads = Gradient (&xsl, &xsf, LOGICAL_FALSE, donor);

		      dot = GeoDotVectorVector (grads, faces[face].d);

		      l1 = GeoMagVector (grads);

		      l2 = GeoMagVector (faces[face].d);

		    }
		  else
		    {

		      acceptor = element;
		      donor = neighbor;

		      grads = Gradient (&xsl, &xsf, LOGICAL_FALSE, donor);

		      dot = GeoDotVectorVector (grads, faces[pair].d);

		      l1 = GeoMagVector (grads);

		      l2 = GeoMagVector (faces[pair].d);

		    }

		  su =
		    LMIN (LMAX (V_GetCmp (&xsl, acceptor) - 2 * dot, 0.0),
			  1.0);

		  Cod = LMIN (V_GetCmp (&Col, donor), 1.0);

		  if (LABS (V_GetCmp (&xsl, acceptor) - su) > SMALL)
		    {

		      sdn =
			(V_GetCmp (&xsl, donor) -
			 su) / (V_GetCmp (&xsl, acceptor) - su);

		      if (sdn >= 0.0 && sdn <= 1.0 && LABS (Cod) > SMALL)
			sjnCBC = LMIN (1, sdn / Cod);
		      else
			sjnCBC = sdn;

		      if (sdn >= 0.0 && sdn <= 1.0)
			sjnUQ =
			  LMIN ((8.0 * Cod * sdn +
				 (1.0 - Cod) * (6.0 * sdn +
						3.0)) / 8.0, sjnCBC);
		      else
			sjnUQ = sdn;

		      if (LABS (l1 * l2) > SMALL)
			ang = LABS (dot / (l1 * l2));
		      else
			ang = LABS (dot / SMALL);

		      if (ang > 1.0)
			ang = 1.0;

		      tetaj = acos (ang);

		      qj =
			LMIN (parameter.kq * 0.5 * (cos (2 * tetaj) + 1.0),
			      1.0);

		      sjn = qj * sjnCBC + (1 - qj) * sjnUQ;

		      if (LABS (1.0 - sdn) > SMALL)
			{
			  betaj =
			    LMIN (LMAX ((sjn - sdn) / (1.0 - sdn), 0.0), 1.0);
			}

		    }

		}

	    }
	  else
	    {

	      if (faces[face].bc == PROCESSOR)
		{

		  if (V_GetCmp (&uf, face) != 0.0)
		    {

		      if (V_GetCmp (&uf, face) > 0.0)
			{

			  acceptor = faces[face].ghost;
			  donor = element;

			  grads.x = 0.0;
			  grads.y = 0.0;
			  grads.z = 0.0;

			  dot = GeoDotVectorVector (grads, faces[face].d);

			  l1 = GeoMagVector (grads);

			  l2 = GeoMagVector (faces[face].d);

			}
		      else
			{

			  acceptor = element;
			  donor = faces[face].ghost;

			  grads.x = 0.0;
			  grads.y = 0.0;
			  grads.z = 0.0;

			  dot = GeoDotVectorVector (grads, faces[face].d);

			  l1 = GeoMagVector (grads);

			  l2 = GeoMagVector (faces[face].d);

			}

		      su = LMIN (LMAX (V_GetCmp (&xsl, acceptor) - 2 * dot, 0.0), 1.0);

		      Cod = LMIN (V_GetCmp (&Col, donor), 1.0);

		      if (LABS (V_GetCmp (&xsl, acceptor) - su) > SMALL)
			{

			  sdn = (V_GetCmp (&xsl, donor) - su) / (V_GetCmp (&xsl, acceptor) - su);

			  if (sdn >= 0.0 && sdn <= 1.0 && LABS (Cod) > SMALL)
			    sjnCBC = LMIN (1, sdn / Cod);
			  else
			    sjnCBC = sdn;

			  if (sdn >= 0.0 && sdn <= 1.0)
			    sjnUQ = LMIN ((8.0 * V_GetCmp (&Col, donor) * sdn + (1.0 - V_GetCmp (&Col, donor)) * (6.0 * sdn + 3.0)) / 8.0, sjnCBC);
			  else
			    sjnUQ = sdn;

			  if (LABS (l1 * l2) > SMALL)
			    ang = LABS (dot / (l1 * l2));
			  else
			    ang = LABS (dot / SMALL);

			  if (ang > 1.0)
			    ang = 1.0;

			  tetaj = acos (ang);

			  qj = LMIN (parameter.kq * 0.5 * (cos (2 * tetaj) + 1.0), 1.0);

			  sjn = qj * sjnCBC + (1 - qj) * sjnUQ;

			  if (LABS (1.0 - sdn) > SMALL)
			    {
			      betaj = LMIN (LMAX ((sjn - sdn) / (1.0 - sdn), 0.0), 1.0);
			    }

			}

		    }

		}

	    }
	}

      V_SetCmp (&betaf, face, betaj);

    }

  VecGhostRestoreLocalForm (xs, &xsl);
  VecGhostRestoreLocalForm (Co, &Col);
    
  VecAssemblyBegin (betaf);
  VecAssemblyEnd (betaf);

}

void
CorrectBeta (double dt)
{

  int i;

  register unsigned int face, pair;
  register unsigned int element, neighbor;
  register unsigned int donor, acceptor;

  double Cj;

  double cbetaj, betaj;

  double ds, Ep, Em;

  VecGhostGetLocalForm (xs0, &xs0l);
  VecGhostGetLocalForm (xs, &xsl);
  
  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

      cbetaj = 0.0;

      betaj = V_GetCmp (&betaf, face);

      if (betaj < 1E-2) continue;

      if (pair != -1)
	{

	  neighbor = faces[pair].element;

	  if (V_GetCmp (&uf, face) != 0.0)
	    {

	      if (V_GetCmp (&uf, face) > 0.0)
		{

		  acceptor = neighbor;
		  donor = element;

		}
	      else
		{

		  acceptor = element;
		  donor = neighbor;

		}

	      Cj = LMIN (LMAX (-V_GetCmp (&uf, face) * faces[face].Aj * dt /
			       elements[element].Vp, 0.0), 1.0);

	      ds =
		0.5 * (V_GetCmp (&xs0l, acceptor) +
		       V_GetCmp (&xsl, acceptor)) - 0.5 *
		(V_GetCmp (&xs0l, donor) + V_GetCmp (&xsl, donor));

	      if (V_GetCmp (&xsl, donor) < 0.0)
		{
		  Em = LMAX (-V_GetCmp (&xsl, donor), 0.0);

		  // Donor value < 0.0 Ex: sd = -0.1 -> Em = +0.1 
		  if (Em > SMALL && Cj > SMALL)
		    {
		      if (ds > Em)
			{
			  cbetaj = Em * (2 + Cj - 2 * Cj * betaj) / (2 * Cj * (ds - Em));

			  cbetaj = LMIN (cbetaj, betaj);
			}
		    }

		}

	      if (V_GetCmp (&xsl, donor) > 1.0)
		{

		  Ep = LMAX (V_GetCmp (&xsl, donor) - 1.0, 0.0);

		  // Donor value > 1.0 Ex: sd = 1.1 -> Ep = +0.1 
		  if (Ep > SMALL && Cj > SMALL)
		    {
		      if (ds < -Ep)
			{
			  cbetaj = Ep * (2 + Cj - 2 * Cj * betaj) / (2 * Cj * (-ds - Ep));

			  cbetaj = LMIN (cbetaj, betaj);
			}
		    }
		}

	    }

	}
      else
	{

	  if (faces[face].bc == PROCESSOR)
	    {

	      if (V_GetCmp (&uf, face) != 0.0)
		{

		  if (V_GetCmp (&uf, face) > 0.0)
		    {

		      acceptor = faces[face].ghost;
		      donor = element;

		    }
		  else
		    {
		      acceptor = element;
		      donor = faces[face].ghost;

		    }

		  Cj = LMAX (-V_GetCmp (&uf, face) * faces[face].Aj * dt / elements[element].Vp, 0.0);

		  ds = 0.5 * (V_GetCmp (&xs0l, acceptor) + V_GetCmp (&xsl, acceptor)) - 0.5 * (V_GetCmp (&xs0l, donor) + V_GetCmp (&xsl, donor));

		  if (V_GetCmp (&xsl, donor) < 0.0)
		    {
		      Em = LMAX (-V_GetCmp (&xsl, donor), 0.0);

		      // Donor value < 0.0 Ex: sd = -0.1 -> Em = +0.1 
		      if (Em > SMALL && Cj > SMALL)
			{
			  if (ds > Em)
			    {
			      cbetaj = Em * (2 + Cj - 2 * Cj * betaj) / (2 * Cj * (ds - Em));

			      cbetaj = LMIN (cbetaj, betaj);
			    }
			}

		    }


		  if (V_GetCmp (&xsl, donor) > 1.0)
		    {

		      Ep = LMAX (V_GetCmp (&xsl, donor) - 1.0, 0.0);

		      // Donor value > 1.0 Ex: sd = 1.1 -> Ep = +0.1 
		      if (Ep > SMALL && Cj > SMALL)
			{
			  if (ds < -Ep)
			    {
			      cbetaj = Ep * (2 + Cj - 2 * Cj * betaj) / (2 * Cj * (-ds - Ep));

			      cbetaj = LMIN (cbetaj, betaj);
			    }
			}
		    }

		}
	    }

	}

      betaj -= cbetaj;

      betaj = LMAX (betaj, 0.0);

      V_SetCmp (&betaf, face, betaj);

    }

  VecGhostRestoreLocalForm (xs0, &xs0l);
  VecGhostRestoreLocalForm (xs, &xsl);

  VecAssemblyBegin (betaf);
  VecAssemblyEnd (betaf);

}

void
BoundScalar (Vec * x, Vec * xl, double min, double max)
{

  int i;

  int element;

  VecGhostGetLocalForm (*x, xl);

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      V_SetCmp (x, elements[element].index, LMAX (LMIN (V_GetCmp (xl, element), max), min));

    }

  VecGhostRestoreLocalForm (*x, xl);

  VecAssemblyBegin (*x);
  VecAssemblyEnd (*x);

  VecGhostUpdateBegin (*x, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (*x, INSERT_VALUES, SCATTER_FORWARD);
  
}

void
SmoothScalar (Vec * x, Vec * xl, int n)
{

  unsigned int i, j, k;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  double sj;

  //double dNf, dPf;
  double lambda;

  double sum1, sum2;

  double *sa;
  double *sm;

  sa = calloc (nbelements, sizeof (double));
  sm = calloc (nbelements, sizeof (double));

  VecGhostGetLocalForm (*x, xl);

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      sa[i] = V_GetCmp (xl, element);

    }

  for (k = 0; k < n; k++)
    {

      for (i = 0; i < nbelements; i++)
	{

	  element = i;

	  sum1 = 0.0;

	  sum2 = 0.0;

	  for (j = 0; j < elements[element].nbfaces; j++)
	    {

	      face = elements[element].face[j];

	      pair = faces[face].pair;

	      if (pair != -1)
		{
		  neighbor = faces[pair].element;

		  /*
		  dNf = GeoMagVector (GeoSubVectorVector (elements[neighbor].celement, faces[pair].cface));
		  dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement, faces[face].cface));

		  lambda = dPf / (dPf + dNf);
		  */

  	          lambda = 0.5;

		  sj = sa[neighbor] * lambda + sa[element] * (1.0 - lambda);

		}
	      else
		{

		  sj = sa[element];

		}

	      sum1 += sj * faces[face].Aj;

	      sum2 += faces[face].Aj;

	    }

	  sm[i] = sum1 / sum2;

	}

      for (i = 0; i < nbelements; i++)
	{

	  sa[i] = sm[i];

	}

    }

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      V_SetCmp (xl, element, sm[i]);

    }

  VecGhostRestoreLocalForm (*x, xl);

  free (sa);
  free (sm);

  VecAssemblyBegin (*x);
  VecAssemblyEnd (*x);

  VecGhostUpdateBegin (*x, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (*x, INSERT_VALUES, SCATTER_FORWARD);

}

void
CalculateMassFraction ()
{

  unsigned int i;

  register unsigned int element;

  double f[2];
  double vol[2];

  vol[0] = 0.0;
  vol[1] = 0.0;

  VecGhostGetLocalForm (xs, &xsl);
  VecGhostGetLocalForm (dens, &densl);
  
  for (i = 0; i < nbelements; i++)
    {
      element = i;

      f[1] = V_GetCmp (&xsl, element);
      f[0] = 1 - f[1];

      vol[0] += f[0] * V_GetCmp (&densl, element) * elements[element].Vp;
      vol[1] += f[1] * V_GetCmp (&densl, element) * elements[element].Vp;
    }

  PetscPrintf (PETSC_COMM_WORLD, "\nMass of fluid %d: %+E kg\n", 0, vol[0]);
  PetscPrintf (PETSC_COMM_WORLD, "\nMass of fluid %d: %+E kg\n", 1, vol[1]);

  VecGhostRestoreLocalForm (xs, &xsl);
  VecGhostRestoreLocalForm (dens, &densl);
  
}

void
BuildVolumeOfFluidMatrix (double dt)
{

  unsigned int i, j, n;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  msh_element ghost;

  double aip;
  double ain[MAXFACES];
  unsigned int ani[MAXFACES];
  double bip;

  double betaj;

  // MatSetValues
  int row;
  int ncols;
  int col[MAXFACES+1];
  int nvals;
  double val[MAXFACES+1];

  VecGhostGetLocalForm (xs0, &xs0l);
  VecGhostGetLocalForm (xs, &xsl);
  
  // Equation: ds/dt + div(s*U) = 0

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      aip = 0;

      bip = 0;

      n = 0;

      for (j = 0; j < elements[element].nbfaces; j++)
	{

	  face = elements[element].face[j];

	  pair = faces[face].pair;

	  if (parameter.scheme[is] == UDS)
	    {

	      // UDS
	      if (V_GetCmp (&uf, face) > 0.0)
		betaj = 0.0;
	      else
		betaj = 1.0;

	    }
	  else
	    {

	      // CDS
	      if (V_GetCmp (&uf, face) > 0.0)
		betaj = V_GetCmp (&betaf, face);
	      else
		betaj = 1.0 - V_GetCmp (&betaf, face);

	    }

	  if (pair != -1)
	    {

	      neighbor = faces[pair].element;

	      aip += 0.5 * (1.0 - betaj) * V_GetCmp (&uf, face) * faces[face].Aj;

	      ain[n] = 0.5 * betaj * V_GetCmp (&uf, face) * faces[face].Aj;

	      ani[n] = elements[neighbor].index;
	      n++;

	      bip += -0.5 * (1.0 - betaj) * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xsl, element);

	      bip += -0.5 * betaj * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xsl, neighbor);

	    }
	  else
	    {

	      if (faces[face].bc == PROCESSOR)
		{

		  ghost.index = faces[face].physreg;

		  aip += 0.5 * (1.0 - betaj) * V_GetCmp (&uf, face) * faces[face].Aj;

		  ain[n] = 0.5 * betaj * V_GetCmp (&uf, face) * faces[face].Aj;

		  ani[n] = ghost.index;
		  n++;

		  bip += -0.5 * (1.0 - betaj) * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xsl, element);

		  bip += -0.5 * betaj * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xsl, faces[face].ghost);

		}
	      else
		{

		  bip += -1.0 * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xsf, face);
		}
	    }

	}

      if (dt > 0.0)
	{

	  aip += elements[element].Vp / dt;

	  bip += elements[element].Vp / dt * V_GetCmp (&xs0l, element);

	}

      if (aip == 0.0 || aip != aip)
	{
	  PetscPrintf (PETSC_COMM_WORLD, "\nError: Problem setting up volume-of-fluid matrix\n");
	  exit (LOGICAL_ERROR);
	}

      /*
      Q_SetEntry (&As, elements[element].index, elements[element].index, aip);

      for (j = 0; j < n; j++)
	{
	  Q_SetEntry (&As, elements[element].index, ani[j], ain[j]);
	}
       */

      ncols = 0;
      nvals = 0;

      row = elements[element].index;

      col[ncols] = elements[element].index;
      ncols++;

      val[nvals] = aip;
      nvals++;

      for (j = 0; j < n; j++)
	{
		col[ncols] = ani[j];
	        ncols++;

		val[nvals] = ain[j];
	        nvals++;
	}

      Q_SetEntries (&As, 1, &row, ncols, col, val);

      V_SetCmp (&bs, elements[element].index, bip);

    }

  VecGhostRestoreLocalForm (xs0, &xs0l);
  VecGhostRestoreLocalForm (xs, &xsl);
    
  MatAssemblyBegin (As, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd (As, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin (bs);
  VecAssemblyEnd (bs);

}

void
SolveVolumeOfFluidExplicit (double dt)
{

  unsigned int i, j, n;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  msh_element ghost;

  double aip;
  double ain[MAXFACES];
  unsigned int ani[MAXFACES];
  double bip;

  double betaj;

  double sums;

  VecGhostUpdateBegin (xs0, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (xs0, INSERT_VALUES, SCATTER_FORWARD);

  VecGhostUpdateBegin (xs, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (xs, INSERT_VALUES, SCATTER_FORWARD);

  VecGhostGetLocalForm (xs0, &xs0l);
  VecGhostGetLocalForm (xs, &xsl);
  
  // Equation: ds/dt + div(s*U) = 0

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      aip = 0;

      bip = 0;

      n = 0;

      for (j = 0; j < elements[element].nbfaces; j++)
	{

	  face = elements[element].face[j];

	  pair = faces[face].pair;

	  if (parameter.scheme[is] == UDS)
	    {

	      // UDS
	      if (V_GetCmp (&uf, face) > 0.0)
		betaj = 0.0;
	      else
		betaj = 1.0;

	    }
	  else
	    {

	      // CDS
	      if (V_GetCmp (&uf, face) > 0.0)
		betaj = V_GetCmp (&betaf, face);
	      else
		betaj = 1.0 - V_GetCmp (&betaf, face);

	    }

	  if (pair != -1)
	    {

	      neighbor = faces[pair].element;

	      aip += 0.5 * (1.0 - betaj) * V_GetCmp (&uf, face) * faces[face].Aj;

	      ain[n] = 0.5 * betaj * V_GetCmp (&uf, face) * faces[face].Aj;

	      ani[n] = neighbor;
	      n++;

	      bip += -0.5 * (1.0 - betaj) * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xsl, element);

	      bip += -0.5 * betaj * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xsl, neighbor);

	    }
	  else
	    {

	      if (faces[face].bc == PROCESSOR)
		{

		  ghost.index = faces[face].physreg;

		  aip += 0.5 * (1.0 - betaj) * V_GetCmp (&uf, face) * faces[face].Aj;

		  ain[n] = 0.5 * betaj * V_GetCmp (&uf, face) * faces[face].Aj;

		  ani[n] = faces[face].ghost;
		  n++;

		  bip += -0.5 * (1.0 - betaj) * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xsl, element);

		  bip += -0.5 * betaj * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xsl, faces[face].ghost);

		}
	      else
		{

		  bip += -1.0 * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xsf, face);
		}
	    }

	}

      if (dt > 0.0)
	{

	  aip += elements[element].Vp / dt;

	  bip += elements[element].Vp / dt * V_GetCmp (&xs0l, element);

	}

      if (aip == 0.0)
	{
	  PetscPrintf (PETSC_COMM_WORLD,
		       "\nError: Problem setting up volume-of-fluid matrix\n");
	  exit (LOGICAL_ERROR);
	}

      sums = 0.0;

      for (j = 0; j < n; j++)
	{
	  sums += ain[j] * V_GetCmp (&xs0l, ani[j]);
	}

      V_SetCmp (&xs, elements[element].index, (bip - sums) / aip);

    }

  VecGhostRestoreLocalForm (xs0, &xs0l);
  VecGhostRestoreLocalForm (xs, &xsl);
    
  VecAssemblyBegin (xs);
  VecAssemblyEnd (xs);

}

void
CalculateGamma (char *var, int *fiter, double dt, double *maxCp, int verbose, int pchecks)
{

  int i, j;

  int loc;
  
  double mres;
  int miter;
  double mtime;

  if (parameter.calc[is] == LOGICAL_FALSE)
    return;

  CalculateMaxCourantNumber (dt, 1);

  VecMax (Co, &loc, maxCp);
  
  parameter.ncicsamsteps = LMIN ((int) (*maxCp * 2) + 1, 100);
 
  if (verbose == LOGICAL_TRUE)
    PetscPrintf (PETSC_COMM_WORLD, "\nMaximum Courant number at interface: %.3f\n", *maxCp);
  
  V_Constr (&betaf, nbfaces, 1);	// CICSAM interpolation factor

  for (i = 0; i < parameter.ncicsamsteps; i++)
    {

      CalculateMaxCourantNumber (dt / parameter.ncicsamsteps, 0);

      VecMax (Co, &loc, maxCp);

      // Store previous time step values
      VecCopy (xs, xs0);

      fiter[is]++;

      // Predict beta - CICSAM       
      PredictBeta ();

      for (j = 0; j <= parameter.ncicsamcor; j++)
	{

	  if (parameter.timemethod[is] == IMPLICITEULER)
	    {
	      Q_Constr (&As, nbelements, LOGICAL_FALSE);
	      V_Constr (&bs, nbelements, 0);	// Gamma source
	    }

	  if (parameter.timemethod[is] == EXPLICITEULER)
	    {

	      // Matrix free VOF 
	      SolveVolumeOfFluidExplicit (dt / parameter.ncicsamsteps);

	    }

	  if (parameter.timemethod[is] == IMPLICITEULER)
	    {

	      // Build VOF matrix     
	      BuildVolumeOfFluidMatrix (dt / parameter.ncicsamsteps);

	      if (pchecks == LOGICAL_TRUE)
		{
		  if (!CheckIfDiagonalMatrix (&As))
		    {
		      PetscPrintf (PETSC_COMM_WORLD, "\nWarning: Volume-of-fluid matrix is not diagonal dominant\n");
		      //MatView(As, PETSC_VIEWER_STDOUT_WORLD);
		      WriteMatrix (&As);
		      WriteVector (&bs);
		      //exit (LOGICAL_ERROR);
		    }
		}

	      // Solve matrix to get indicator function s
	      SolveMatrix (&As, &xs, &bs, &miter, &mres, &mtime, parameter.msolver[is], parameter.mprecond[is], parameter.miter[is], parameter.mtol[is]);

	      if (verbose == LOGICAL_TRUE)
		PetscPrintf (PETSC_COMM_WORLD,
			     "\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
			     var[is], miter, mres, mtime);

	      if (pchecks == LOGICAL_TRUE)
		{
		  if (mres > parameter.mtol[is]
		      && miter == parameter.miter[is])
		    {
		      PetscPrintf (PETSC_COMM_WORLD,
				   "\nError: Problem solving matrix %c\n",
				   var[is]);
		      exit (LOGICAL_ERROR);
		    }
		}

	    }

	  // Correct beta
	  CorrectBeta (dt / parameter.ncicsamsteps);

	  CorrectFaceS ();

	  if (parameter.timemethod[is] == IMPLICITEULER)
	    {
	      Q_Destr (&As);
	      V_Destr (&bs);
	    }

	}

      if (parameter.fill == LOGICAL_TRUE)
	{
	  // Bound indicator function      
	  BoundScalar (&xs, &xsl, 0.0, 1.0);
	}

      //SmoothScalar (&xs, &xsl, 2);

    }

  if (pchecks == LOGICAL_TRUE)
    {
      // Calculate mass fractions             
      CalculateMassFraction ();
    }

  V_Destr (&betaf);

}
