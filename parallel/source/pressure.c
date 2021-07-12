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

#include "pressure.h"

void
CorrectFaceP () //修正面心压力值
{

  unsigned int i;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  msh_element ghost;

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

  VecGhostGetLocalForm (ap, &apl);

  VecGhostGetLocalForm (xp, &xpl);
  VecGhostGetLocalForm (xs, &xsl);

  VecGhostGetLocalForm (dens, &densl);
  
  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

      if (parameter.orthof != 0.0)
	gradpp = Gradient (&xpl, &xpf, LOGICAL_TRUE, element);

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

	  phij = V_GetCmp (&xpl, neighbor) * lambda +
	    V_GetCmp (&xpl, element) * (1.0 - lambda);

	  V_SetCmp (&xpf, face, phij);

	}
      else
	{

	  if (faces[face].bc == PROCESSOR)
	    {

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

	      phij = V_GetCmp (&xpl, faces[face].ghost) * lambda +
		V_GetCmp (&xpl, element) * (1.0 - lambda);

	      V_SetCmp (&xpf, face, phij);

	    }
	  else
	    {

	      ppl = V_GetCmp (&xpl, element);

	      apj = V_GetCmp (&apl, element);

	      if (parameter.orthof != 0.0)
		{
		  ppl += parameter.orthof *
		    GeoDotVectorVector (gradpp,
					GeoSubVectorVector (faces[face].rpl,
							    elements[element].
							    celement));
		}

	      ghf =
		V_GetCmp (&densl, element) * GeoDotVectorVector (g, faces[face].d);

	      ppl += ghf;

	      if (faces[face].bc == PERMEABLE)
		{

		  V_SetCmp (&xpf, face,
			    V_GetCmp (&xpf,
				      face) * (1.0 - V_GetCmp (&xsl,
							       element)) +
			    ppl * V_GetCmp (&xsl, element));

		}

	      if (faces[face].bc == INLET)
		{

		  V_SetCmp (&xpf, face, ppl - V_GetCmp (&uf, face) * apj * (faces[face].dj + faces[face].kj));

		}

	      if (faces[face].bc == MOVINGWALL ||
		  faces[face].bc == WALL ||
		  faces[face].bc == ADIABATICWALL ||
		  faces[face].bc == SURFACE)
		{

		  V_SetCmp (&xpf, face, ppl);

		}
	    }

	}

    }

  VecGhostRestoreLocalForm (ap, &apl);

  VecGhostRestoreLocalForm (xp, &xpl);
  VecGhostRestoreLocalForm (xs, &xsl);

  VecGhostRestoreLocalForm (dens, &densl);
    
  VecAssemblyBegin (xpf);
  VecAssemblyEnd (xpf);

}

void
BuildContinuityMatrix (double dt) //建立连续性方程的矩阵
{

  int i, j, n;

  register unsigned int  face, pair;
  register unsigned int element, neighbor;

  msh_element ghost;

  double acp;
  double acn[MAXFACES];
  unsigned int ani[MAXFACES];
  double bcp;

  double apj;

  double Huj, Hvj, Hwj;
  double Hf;

  //double dNf, dPf;
  double lambda;

  double dj;

  msh_vector gradpp;
  msh_vector gradpn;

  // MatSetValues
  int row;
  int ncols;
  int col[MAXFACES+1];
  int nvals;
  double val[MAXFACES+1];

  VecGhostGetLocalForm (ap, &apl);
  VecGhostGetLocalForm (hu, &hul);
  VecGhostGetLocalForm (hv, &hvl);
  VecGhostGetLocalForm (hw, &hwl);

  VecGhostGetLocalForm (xu, &xul);
  VecGhostGetLocalForm (xv, &xvl);
  VecGhostGetLocalForm (xw, &xwl);
  VecGhostGetLocalForm (xp, &xpl);
  VecGhostGetLocalForm (xs, &xsl);

  VecGhostGetLocalForm (dens, &densl);
  
  // Equation: div(U) = 0

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      acp = 0.0;

      bcp = 0.0;

      n = 0;

      if (parameter.orthof != 0.0)
	gradpp = Gradient (&xpl, &xpf, LOGICAL_TRUE, element);

      for (j = 0; j < elements[element].nbfaces; j++)
	{

	  face = elements[element].face[j];

	  pair = faces[face].pair;

	  if (pair != -1)
	    {

	      neighbor = faces[pair].element;

	      /*
	      dNf =
		GeoMagVector (GeoSubVectorVector
			      (elements[neighbor].celement,
			       faces[face].cface));
	      dPf =
		GeoMagVector (GeoSubVectorVector
			      (elements[element].celement,
			       faces[face].cface));

	      lambda = dPf / (dPf + dNf);
	      */

	      lambda = 0.5;

	      apj =
		V_GetCmp (&apl, neighbor) * lambda +
		V_GetCmp (&apl, element) * (1.0 - lambda);

	      Huj =
		V_GetCmp (&hul, neighbor) * lambda +
		V_GetCmp (&hul, element) * (1.0 - lambda);

	      Hvj =
		V_GetCmp (&hvl, neighbor) * lambda +
		V_GetCmp (&hvl, element) * (1.0 - lambda);

	      Hwj =
		V_GetCmp (&hwl, neighbor) * lambda +
		V_GetCmp (&hwl, element) * (1.0 - lambda);

	      Hf = Huj * faces[face].n.x +
		Hvj * faces[face].n.y + Hwj * faces[face].n.z;

	      acp +=
		-1.0 / (apj * (faces[face].dj + faces[face].kj)) *
		faces[face].Aj;

	      acn[n] =
		1.0 / (apj * (faces[face].dj + faces[face].kj)) *
		faces[face].Aj;

	      ani[n] = elements[neighbor].index;
	      n++;

	      V_SetCmp (&uf, face, 1.0 / apj * Hf);

	      bcp += V_GetCmp (&uf, face) * faces[face].Aj;

	      if (parameter.orthof != 0.0)
		gradpn = Gradient (&xpl, &xpf, LOGICAL_TRUE, neighbor);

	      if (parameter.orthof != 0.0)
		{
		  // Non-orthogonal correction term
		  bcp +=
		    -1.0 * parameter.orthof / (apj *
					       (faces[face].dj +
						faces[face].kj)) *
		    faces[face].Aj *
		    (GeoDotVectorVector
		     (gradpn,
		      GeoSubVectorVector (faces[face].rnl,
					  elements[neighbor].celement)) -
		     GeoDotVectorVector (gradpp,
					 GeoSubVectorVector (faces[face].rpl,
							     elements
							     [element].
							     celement)));
		}

	    }
	  else
	    {

	      if (faces[face].bc == PROCESSOR)
		{

		  ghost.index = faces[face].physreg;

		  ghost.celement.x = V_GetCmp (&cexl, faces[face].ghost);
		  ghost.celement.y = V_GetCmp (&ceyl, faces[face].ghost);
		  ghost.celement.z = V_GetCmp (&cezl, faces[face].ghost);

		  dj = GeoMagVector (GeoSubVectorVector
				     (ghost.celement,
				      elements[element].celement));

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

		  apj =
		    V_GetCmp (&apl, faces[face].ghost) * lambda +
		    V_GetCmp (&apl, element) * (1.0 - lambda);

		  Huj =
		    V_GetCmp (&hul, faces[face].ghost) * lambda +
		    V_GetCmp (&hul, element) * (1.0 - lambda);

		  Hvj =
		    V_GetCmp (&hvl, faces[face].ghost) * lambda +
		    V_GetCmp (&hvl, element) * (1.0 - lambda);

		  Hwj =
		    V_GetCmp (&hwl, faces[face].ghost) * lambda +
		    V_GetCmp (&hwl, element) * (1.0 - lambda);

		  Hf = Huj * faces[face].n.x +
		    Hvj * faces[face].n.y + Hwj * faces[face].n.z;

		  acp += -1.0 / (apj * dj) * faces[face].Aj;

		  acn[n] = 1.0 / (apj * dj) * faces[face].Aj;

		  ani[n] = ghost.index;
		  n++;

		  V_SetCmp (&uf, face, 1.0 / apj * Hf);

		  bcp += V_GetCmp (&uf, face) * faces[face].Aj;

		  /*
		     if (parameter.orthof != 0.0)                  
		     gradpn = Gradient (&xpl, &xpf, LOGICAL_TRUE,  faces[face].ghost);

		     if (parameter.orthof != 0.0)
		     {
		     // Non-orthogonal correction term
		     bcp += -1.0 * parameter.orthof / (apj * dj) * faces[face].Aj *
		     (GeoDotVectorVector
		     (gradpn,
		     GeoSubVectorVector (faces[face].rnl,
		     ghost.celement)) -
		     GeoDotVectorVector (gradpp,
		     GeoSubVectorVector (faces[face].rpl,
		     elements[element].
		     celement))); 
		     }
		   */

		}
	      else
		{

		  apj = V_GetCmp (&apl, element);

		  Huj = V_GetCmp (&hul, element);
		  Hvj = V_GetCmp (&hvl, element);
		  Hwj = V_GetCmp (&hwl, element);

		  Hf =
		    Huj / apj * faces[face].n.x +
		    Hvj / apj * faces[face].n.y + Hwj / apj * faces[face].n.z;

		  if (faces[face].bc == PERMEABLE)
		    {

		      acp +=
			-1.0 / (apj * (faces[face].dj + faces[face].kj)) *
			faces[face].Aj * (1.0 - V_GetCmp (&xsl, element));

		      bcp +=
			-1.0 / (apj * (faces[face].dj + faces[face].kj)) *
			faces[face].Aj * V_GetCmp (&xpf,
						   face) * (1.0 -
							    V_GetCmp (&xsl,
								      element));

		      V_SetCmp (&uf, face,
				1.0 / apj * Hf * (1.0 -
						  V_GetCmp (&xsl, element)));

		      bcp += V_GetCmp (&uf, face) * faces[face].Aj;

		      if (parameter.orthof != 0.0)
			{
			  // Non-orthogonal correction term               
			  bcp +=
			    1.0 * parameter.orthof / (apj *
						      (faces[face].dj +
						       faces[face].kj)) *
			    (GeoDotVectorVector
			     (gradpp,
			      GeoSubVectorVector (faces[face].rpl,
						  elements[element].
						  celement))) *
			    faces[face].Aj * (1.0 - V_GetCmp (&xsl, element));
			}

		    }

		  if (faces[face].bc == OUTLET)
		    {

		      // velocity gradient = 0
		      // specified pressure

		      acp +=
			-1.0 / (apj * (faces[face].dj + faces[face].kj)) *
			faces[face].Aj;

		      bcp +=
			-1.0 / (apj * (faces[face].dj + faces[face].kj)) *
			faces[face].Aj * V_GetCmp (&xpf, face);

		      V_SetCmp (&uf, face, Hf / apj);

		      bcp += V_GetCmp (&uf, face) * faces[face].Aj;

		      if (parameter.orthof != 0.0)
			{
			  // Non-orthogonal correction term
			  bcp +=
			    1.0 * parameter.orthof / (apj *
						      (faces[face].dj +
						       faces[face].kj)) *
			    (GeoDotVectorVector
			     (gradpp,
			      GeoSubVectorVector (faces[face].rpl,
						  elements[element].
						  celement))) *
			    faces[face].Aj;
			}

		    }

		  if (faces[face].bc == PRESSUREINLET)
		    {

		      // specified pressure 
		      // velocity gradient = 0

		      acp +=
			-1.0 / (apj * (faces[face].dj + faces[face].kj)) *
			faces[face].Aj;

		      bcp +=
			-1.0 / (apj * (faces[face].dj + faces[face].kj)) *
			faces[face].Aj * V_GetCmp (&xpf, face);

		      V_SetCmp (&uf, face, 1.0 / apj * Hf);

		      bcp += V_GetCmp (&uf, face) * faces[face].Aj;

		      if (parameter.orthof != 0.0)
			{
			  // Non-orthogonal correction term
			  bcp +=
			    1.0 * parameter.orthof / (apj *
						      (faces[face].dj +
						       faces[face].kj)) *
			    (GeoDotVectorVector
			     (gradpp,
			      GeoSubVectorVector (faces[face].rpl,
						  elements[element].
						  celement))) *
			    faces[face].Aj;
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
		      V_SetCmp (&uf, face,
				V_GetCmp (&xuf, face) * faces[face].n.x +
				V_GetCmp (&xvf, face) * faces[face].n.y +
				V_GetCmp (&xwf, face) * faces[face].n.z);

		      bcp += V_GetCmp (&uf, face) * faces[face].Aj;

		    }

		}

	    }

	}

      if (acp == 0.0 || acp != acp)
	{
	  PetscPrintf (PETSC_COMM_WORLD, "\nError: Problem setting up continuity matrix\n");
	  exit (LOGICAL_ERROR);
	}
	
      /*
      Q_SetEntry (&Ac, elements[element].index, elements[element].index, acp);

      for (j = 0; j < n; j++)
	{
	  if (ani[j] > elements[element].index)
	    {
	      Q_SetEntry (&Ac, elements[element].index, ani[j], acn[j]);
	    }
	}
      */

      ncols = 0;
      nvals = 0;

      row = elements[element].index;

      col[ncols] = elements[element].index;
      ncols++;

      val[nvals] = acp;
      nvals++;

      for (j = 0; j < n; j++)
	{
	  if (ani[j] > elements[element].index)
	    {
      		col[ncols] = ani[j];
	        ncols++;

      		val[nvals] = acn[j];
	        nvals++;
	    }
	}

      Q_SetEntries (&Ac, 1, &row, ncols, col, val);

      V_SetCmp (&bp, elements[element].index, bcp);

    }

  VecGhostRestoreLocalForm (ap, &apl);
  VecGhostRestoreLocalForm (hu, &hul);
  VecGhostRestoreLocalForm (hv, &hvl);
  VecGhostRestoreLocalForm (hw, &hwl);

  VecGhostRestoreLocalForm (xu, &xul);
  VecGhostRestoreLocalForm (xv, &xvl);
  VecGhostRestoreLocalForm (xw, &xwl);
  VecGhostRestoreLocalForm (xp, &xpl);
  VecGhostRestoreLocalForm (xs, &xsl);

  VecGhostRestoreLocalForm (dens, &densl);

  MatAssemblyBegin (Ac, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd (Ac, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin (bp);
  VecAssemblyEnd (bp);

}

void CalculatePressure (char *var, int *fiter, double dt, double maxCp, int verbose, int pchecks)
{

  int i;

  double mres;
  int miter;
  double mtime;

  double presc;

  if (parameter.calc[ip] == LOGICAL_FALSE)
    return;
    
  V_Constr (&xpp, nbelements, 0);	// Pressure at cell center - previous iteration

  // Store previous time step values
  VecCopy (xp, xp0);

  fiter[ip]++;
   
  for (i = 0; i <= parameter.northocor; i++)
    {

      Q_Constr (&Ac, nbelements, LOGICAL_TRUE);
      V_Constr (&bp, nbelements, 0);	// Continuity source

      // Store previous iteration values
      VecCopy (xp, xpp);

      // Build the continuity matrix (mass conservation)          
      BuildContinuityMatrix (dt);

      if (pchecks == LOGICAL_TRUE)
	{
	  if (!CheckIfDiagonalMatrix (&Ac))
	    {
	      PetscPrintf (PETSC_COMM_WORLD,
			   "\nWarning: Continuity matrix is not diagonal dominant\n");
	      //MatView(Ac, PETSC_VIEWER_STDOUT_WORLD);
	      WriteMatrix (&Ac);
	      WriteVector (&bp);
	      //exit (LOGICAL_ERROR);
	    }
	}

      // Solve matrix to get pressure p
      SolveMatrix (&Ac, &xp, &bp, &miter, &mres, &mtime,
		   parameter.msolver[ip], parameter.mprecond[ip],
		   parameter.miter[ip], parameter.mtol[ip]);

      if (verbose == LOGICAL_TRUE)
	PetscPrintf (PETSC_COMM_WORLD,
		     "\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
		     var[ip], miter, mres, mtime);

      if (pchecks == LOGICAL_TRUE)
	{
	  if (mres > parameter.mtol[ip] && miter == parameter.miter[ip])
	    {
	      PetscPrintf (PETSC_COMM_WORLD,
			   "\nError: Problem solving matrix %c\n", var[ip]);
	      exit (LOGICAL_ERROR);
	    }
	}

      // Calculate pressure convergence
      VecWAXPY (temp1, -1.0, xp, xpp);
      VecNorm (temp1, NORM_2, &presc);

      if (verbose == LOGICAL_TRUE)
	PetscPrintf (PETSC_COMM_WORLD, "\nNon-orthogonality error (continuity): %+E\n", presc);

      if (presc < parameter.mtol[ip])
	break;

      CorrectFaceP ();

      Q_Destr (&Ac);
      V_Destr (&bp);

    }

  V_Destr (&xpp);
    
}
