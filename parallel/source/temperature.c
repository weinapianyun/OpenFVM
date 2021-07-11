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

#include "temperature.h"

void
CorrectFaceT ()
{

  int i;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  msh_element ghost;

  //double dNf, dPf;
  double lambda;

  double Tpl;

  msh_vector gradTp;

  VecGhostGetLocalForm (xT, &xTl);
  
  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

      if (parameter.orthof != 0.0)
	gradTp = Gradient (&xTl, &xTf, LOGICAL_TRUE, element);

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

	  V_SetCmp (&xTf, face,
		    V_GetCmp (&xTl, neighbor) * lambda +
		    V_GetCmp (&xTl, element) * (1.0 - lambda));

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

	      V_SetCmp (&xTf, face,
			V_GetCmp (&xTl, faces[face].ghost) * lambda +
			V_GetCmp (&xTl, element) * (1.0 - lambda));

	    }
	  else
	    {

	      Tpl = V_GetCmp (&xTl, element);

	      if (parameter.orthof != 0.0)
		{
		  // Non-orthogonal correction
		  Tpl += parameter.orthof * GeoDotVectorVector (gradTp,
								GeoSubVectorVector
								(faces
								 [face].
								 rpl,
								 elements
								 [element].
								 celement));
		}

	      if (faces[face].bc == ADIABATICWALL)
		{

		  V_SetCmp (&xTf, face, Tpl);

		}

	    }
	}

    }

  VecGhostRestoreLocalForm (xT, &xTl);
    
  VecAssemblyBegin (xTf);
  VecAssemblyEnd (xTf);

}

void
BuildEnergyMatrix (double dt, double schemefactor)
{

  unsigned int i, j, n;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  msh_element ghost;

  double aep;
  double aen[MAXFACES];
  unsigned int ani[MAXFACES];

  double bep;

  //double dNf, dPf;
  double lambda;

  double dj;

  double xsi;

  msh_vector gradTp;
  msh_vector gradTn;

  //msh_vector gradup, gradvp, gradwp;

  double densp;
  double spheatp;
  double thcondj;

  // MatSetValues
  int row;
  int ncols;
  int col[MAXFACES+1];
  int nvals;
  double val[MAXFACES+1];

  VecGhostGetLocalForm (xT0, &xT0l);
  VecGhostGetLocalForm (xT, &xTl);
	  
  VecGhostGetLocalForm (dens, &densl);
  VecGhostGetLocalForm (visc, &viscl);
  VecGhostGetLocalForm (spheat, &spheatl);
  VecGhostGetLocalForm (thcond, &thcondl);
  
  for (i = 0; i < nbelements; i++)
    {

      element = i;

      aep = 0.0;
      bep = 0.0;

      n = 0;

      if (parameter.orthof != 0.0)
	gradTp = Gradient (&xT0l, &xTf, LOGICAL_TRUE, element);

      densp = V_GetCmp (&densl, element);
      spheatp = V_GetCmp (&spheatl, element);

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

	      thcondj =
		V_GetCmp (&thcondl,
			  element) * (1.0 - lambda) + V_GetCmp (&thcondl,
								neighbor) *
		lambda;

	      // Conduction 
	      aep += schemefactor *
		thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
		elements[element].Vp;
	      aen[n] = schemefactor *
		-thcondj * faces[face].Aj / (faces[face].dj +
					     faces[face].kj) /
		elements[element].Vp;

	      // Convection 
	      if (parameter.scheme[iT] == UDS)
		{
		  // UDS
		  if (V_GetCmp (&uf, face) > 0.0)
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
	      aep += schemefactor *
		(1.0 - xsi) * densp * spheatp *
		V_GetCmp (&uf, face) * faces[face].Aj / elements[element].Vp;
	      aen[n] += schemefactor *
		xsi * densp * spheatp *
		V_GetCmp (&uf, face) * faces[face].Aj / elements[element].Vp;

	      ani[n] = elements[neighbor].index;
	      n++;

	      // Conduction
	      bep += -(1.0 - schemefactor) *
		thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
		elements[element].Vp * V_GetCmp (&xT0l, element);

	      bep += +(1.0 - schemefactor) *
		thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
		elements[element].Vp * V_GetCmp (&xT0l, neighbor);

	      // Convection
	      bep += -(1.0 - schemefactor) *
		(1.0 - xsi) * densp * spheatp * V_GetCmp (&uf, face) *
		faces[face].Aj / elements[element].Vp * V_GetCmp (&xT0l,
								  element);

	      bep += -(1.0 - schemefactor) *
		xsi * densp * spheatp * V_GetCmp (&uf, face) *
		faces[face].Aj / elements[element].Vp * V_GetCmp (&xT0l,
								  neighbor);

	      if (parameter.orthof != 0.0)
		gradTn = Gradient (&xT0l, &xTf, LOGICAL_TRUE, neighbor);

	      if (parameter.orthof != 0.0)
		{
		  // Non-orthogonal correction term
		  bep += parameter.orthof *
		    thcondj * faces[face].Aj / (faces[face].dj +
						faces[face].kj) /
		    elements[element].Vp *
		    (GeoDotVectorVector
		     (gradTn,
		      GeoSubVectorVector (faces[face].rnl,
					  elements[neighbor].celement)) -
		     GeoDotVectorVector (gradTp,
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

		  thcondj =
		    V_GetCmp (&thcondl,
			      element) * (1.0 - lambda) + V_GetCmp (&thcondl,
								    faces
								    [face].
								    ghost) *
		    lambda;

		  // Conduction 
		  aep += schemefactor *
		    thcondj * faces[face].Aj / dj / elements[element].Vp;
		  aen[n] = schemefactor *
		    -thcondj * faces[face].Aj / dj / elements[element].Vp;

		  // Convection 
		  if (parameter.scheme[iT] == UDS)
		    {
		      // UDS
		      if (V_GetCmp (&uf, face) > 0.0)
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
		  aep += schemefactor *
		    (1.0 - xsi) * densp * spheatp *
		    V_GetCmp (&uf, face) *
		    faces[face].Aj / elements[element].Vp;
		  aen[n] += schemefactor *
		    xsi * densp * spheatp *
		    V_GetCmp (&uf, face) * faces[face].Aj /
		    elements[element].Vp;

		  ani[n] = ghost.index;
		  n++;

		  // Conduction
		  bep += -(1.0 - schemefactor) *
		    thcondj * faces[face].Aj / (faces[face].dj +
						faces[face].kj) /
		    elements[element].Vp * V_GetCmp (&xT0l, element);

		  bep += +(1.0 - schemefactor) *
		    thcondj * faces[face].Aj / (faces[face].dj +
						faces[face].kj) /
		    elements[element].Vp * V_GetCmp (&xT0l,
						     faces[face].ghost);

		  // Convection
		  bep += -(1.0 - schemefactor) *
		    (1.0 - xsi) * densp * spheatp * V_GetCmp (&uf, face) *
		    faces[face].Aj / elements[element].Vp * V_GetCmp (&xT0l,
								      element);

		  bep += -(1.0 - schemefactor) *
		    xsi * densp * spheatp * V_GetCmp (&uf, face) *
		    faces[face].Aj / elements[element].Vp * V_GetCmp (&xT0l,
								      faces
								      [face].
								      ghost);

		  /*
		     if (parameter.orthof != 0.0)                   
		     gradTn = Gradient (&xTl, &xTf, LOGICAL_TRUE, faces[face].ghost);

		     if (parameter.orthof != 0.0)
		     {
		     // Non-orthogonal correction term
		     bep += parameter.orthof *
		     thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
		     elements[element].Vp *
		     (GeoDotVectorVector
		     (gradTn,
		     GeoSubVectorVector (faces[face].rnl,
		     ghost.celement)) -
		     GeoDotVectorVector (gradTp,
		     GeoSubVectorVector (faces[face].rpl,
		     elements[element].
		     celement)));
		     }
		   */


		}
	      else
		{

		  thcondj = material.bthcond;

		  if (faces[face].bc != EMPTY
		      && faces[face].bc != ADIABATICWALL)
		    {

		      // Conduction
		      aep += schemefactor *
			thcondj * faces[face].Aj / (faces[face].dj +
						    faces[face].kj) /
			elements[element].Vp;
		      bep +=
			thcondj * faces[face].Aj / (faces[face].dj +
						    faces[face].kj) /
			elements[element].Vp * V_GetCmp (&xTf, face);

		      bep += -(1.0 - schemefactor) *
			thcondj * faces[face].Aj / (faces[face].dj +
						    faces[face].kj) /
			elements[element].Vp * V_GetCmp (&xT0l, element);

		      // Convection
		      if (parameter.scheme[iT] == UDS)
			{
			  // UDS
			  if (V_GetCmp (&uf, face) > 0.0)
			    {
			      aep +=
				densp * spheatp *
				V_GetCmp (&uf, face) *
				faces[face].Aj / elements[element].Vp;
			    }
			  else
			    {
			      bep +=
				-densp * spheatp *
				V_GetCmp (&uf, face) *
				faces[face].Aj / elements[element].Vp *
				V_GetCmp (&xTf, face);
			    }

			}
		      else
			{
			  // CDS
			  bep +=
			    -densp * spheatp *
			    V_GetCmp (&uf, face) * faces[face].Aj /
			    elements[element].Vp * V_GetCmp (&xTf, face);
			}

		      if (parameter.orthof != 0.0)
			{
			  // Non-orthogonal correction term
			  bep += parameter.orthof *
			    thcondj * faces[face].Aj / (faces[face].dj +
							faces[face].kj) /
			    elements[element].Vp * GeoDotVectorVector (gradTp,
								       GeoSubVectorVector
								       (faces
									[face].
									rpl,
									elements
									[element].
									celement));
			}
		    }
		}

	    }

	}

      // Unsteady term
      if (dt > 0)
	{
	  aep += densp * spheatp / dt;

	  bep += densp * spheatp / dt * V_GetCmp (&xT0l, element);
	}

      /*
         gradup = Gradient (&xul, &xuf, LOGICAL_TRUE, element);
         gradvp = Gradient (&xvl, &xvf, LOGICAL_TRUE, element);
         gradwp = Gradient (&xwl, &xwf, LOGICAL_TRUE, element);
       */

      if (aep == 0.0 || aep != aep)
	{
	  PetscPrintf (PETSC_COMM_WORLD, "\nError: Problem setting up energy matrix\n");
	  exit (LOGICAL_ERROR);
	}

      /*
      Q_SetEntry (&Ae, elements[element].index, elements[element].index, aep);

      for (j = 0; j < n; j++)
	{
	  Q_SetEntry (&Ae, elements[element].index, ani[j], aen[j]);
	}
      */
		
      ncols = 0;
      nvals = 0;

      row = elements[element].index;

      col[ncols] = elements[element].index;
      ncols++;

      val[nvals] = aep;
      nvals++;

      for (j = 0; j < n; j++)
	{
		col[ncols] = ani[j];
	        ncols++;

		val[nvals] = aen[j];
		nvals++;
	}

      Q_SetEntries (&Ae, 1, &row, ncols, col, val);

      V_SetCmp (&bT, elements[element].index, bep);

    }

  VecGhostRestoreLocalForm (xT0, &xT0l);
  VecGhostRestoreLocalForm (xT, &xTl);
	  	  
  VecGhostRestoreLocalForm (dens, &densl);
  VecGhostRestoreLocalForm (visc, &viscl);
  VecGhostRestoreLocalForm (spheat, &spheatl);
  VecGhostRestoreLocalForm (thcond, &thcondl);
    
  MatAssemblyBegin (Ae, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd (Ae, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin (bT);
  VecAssemblyEnd (bT);

}

void
SolveEnergyExplicit (double dt)
{

  unsigned int i, j, n;

  register unsigned int face, pair;
  register unsigned element, neighbor;

  msh_element ghost;

  double aep;
  double aen[MAXFACES];
  unsigned int ani[MAXFACES];

  double bep;

  //double dNf, dPf;
  double lambda;

  double dj;

  double xsi;

  msh_vector gradTp;
  msh_vector gradTn;

  //msh_vector gradup, gradvp, gradwp;

  double densp;
  double spheatp;
  double thcondj;

  double sumT;

  VecGhostUpdateBegin (xT0, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (xT0, INSERT_VALUES, SCATTER_FORWARD);

  VecGhostGetLocalForm (xT0, &xT0l);
  VecGhostGetLocalForm (xT, &xTl);
	  
  VecGhostGetLocalForm (dens, &densl);
  VecGhostGetLocalForm (visc, &viscl);
  VecGhostGetLocalForm (spheat, &spheatl);
  VecGhostGetLocalForm (thcond, &thcondl);
  
  for (i = 0; i < nbelements; i++)
    {

      element = i;

      aep = 0.0;
      bep = 0.0;

      n = 0;

      if (parameter.orthof != 0.0)
	gradTp = Gradient (&xTl, &xTf, LOGICAL_TRUE, element);

      densp = V_GetCmp (&densl, element);
      spheatp = V_GetCmp (&spheatl, element);

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

	      thcondj =
		V_GetCmp (&thcondl,
			  element) * (1.0 - lambda) + V_GetCmp (&thcondl,
								neighbor) *
		lambda;

	      // Conduction 
	      aep +=
		thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
		elements[element].Vp;
	      aen[n] =
		-thcondj * faces[face].Aj / (faces[face].dj +
					     faces[face].kj) /
		elements[element].Vp;

	      // Convection 
	      if (parameter.scheme[iT] == UDS)
		{
		  // UDS
		  if (V_GetCmp (&uf, face) > 0.0)
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
	      aep +=
		(1.0 - xsi) * densp * spheatp *
		V_GetCmp (&uf, face) * faces[face].Aj / elements[element].Vp;
	      aen[n] +=
		xsi * densp * spheatp *
		V_GetCmp (&uf, face) * faces[face].Aj / elements[element].Vp;

	      ani[n] = neighbor;
	      n++;

	      if (parameter.orthof != 0.0)
		gradTn = Gradient (&xTl, &xTf, LOGICAL_TRUE, neighbor);

	      if (parameter.orthof != 0.0)
		{
		  // Non-orthogonal correction term
		  bep += parameter.orthof *
		    thcondj * faces[face].Aj / (faces[face].dj +
						faces[face].kj) /
		    elements[element].Vp *
		    (GeoDotVectorVector
		     (gradTn,
		      GeoSubVectorVector (faces[face].rnl,
					  elements[neighbor].celement)) -
		     GeoDotVectorVector (gradTp,
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

		  thcondj =
		    V_GetCmp (&thcondl,
			      element) * (1.0 - lambda) + V_GetCmp (&thcondl,
								    faces
								    [face].
								    ghost) *
		    lambda;

		  // Conduction 
		  aep += thcondj * faces[face].Aj / dj / elements[element].Vp;
		  aen[n] =
		    -thcondj * faces[face].Aj / dj / elements[element].Vp;

		  // Convection 
		  if (parameter.scheme[iT] == UDS)
		    {
		      // UDS
		      if (V_GetCmp (&uf, face) > 0.0)
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
		  aep +=
		    (1.0 - xsi) * densp * spheatp *
		    V_GetCmp (&uf, face) *
		    faces[face].Aj / elements[element].Vp;
		  aen[n] +=
		    xsi * densp * spheatp *
		    V_GetCmp (&uf, face) * faces[face].Aj /
		    elements[element].Vp;

		  ani[n] = faces[face].ghost;
		  n++;

		  /*
		     if (parameter.orthof != 0.0)                   
		     gradTn = Gradient (&xTl, &xTf, LOGICAL_TRUE, faces[face].ghost);

		     if (parameter.orthof != 0.0)
		     {
		     // Non-orthogonal correction term
		     bep += parameter.orthof *
		     thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
		     elements[element].Vp *
		     (GeoDotVectorVector
		     (gradTn,
		     GeoSubVectorVector (faces[face].rnl,
		     ghost.celement)) -
		     GeoDotVectorVector (gradTp,
		     GeoSubVectorVector (faces[face].rpl,
		     elements[element].
		     celement)));
		     }
		   */


		}
	      else
		{

		  thcondj = material.bthcond;

		  if (faces[face].bc != EMPTY
		      && faces[face].bc != ADIABATICWALL)
		    {

		      // Conduction
		      aep +=
			thcondj * faces[face].Aj / (faces[face].dj +
						    faces[face].kj) /
			elements[element].Vp;
		      bep +=
			thcondj * faces[face].Aj / (faces[face].dj +
						    faces[face].kj) /
			elements[element].Vp * V_GetCmp (&xTf, face);

		      // Convection
		      if (parameter.scheme[iT] == UDS)
			{
			  // UDS
			  if (V_GetCmp (&uf, face) > 0.0)
			    {
			      aep +=
				densp * spheatp *
				V_GetCmp (&uf, face) *
				faces[face].Aj / elements[element].Vp;
			    }
			  else
			    {
			      bep +=
				-densp * spheatp *
				V_GetCmp (&uf, face) *
				faces[face].Aj / elements[element].Vp *
				V_GetCmp (&xTf, face);
			    }

			}
		      else
			{
			  // CDS
			  bep +=
			    -densp * spheatp *
			    V_GetCmp (&uf, face) * faces[face].Aj /
			    elements[element].Vp * V_GetCmp (&xTf, face);
			}

		      if (parameter.orthof != 0.0)
			{
			  // Non-orthogonal correction term
			  bep += parameter.orthof *
			    thcondj * faces[face].Aj / (faces[face].dj +
							faces[face].kj) /
			    elements[element].Vp * GeoDotVectorVector (gradTp,
								       GeoSubVectorVector
								       (faces
									[face].
									rpl,
									elements
									[element].
									celement));
			}
		    }
		}

	    }

	}

      // Unsteady term
      if (dt > 0)
	{
	  aep += densp * spheatp / dt;

	  bep += densp * spheatp / dt * V_GetCmp (&xT0l, element);
	}

      /*
         gradup = Gradient (&xul, &xuf, LOGICAL_TRUE, element);
         gradvp = Gradient (&xvl, &xvf, LOGICAL_TRUE, element);
         gradwp = Gradient (&xwl, &xwf, LOGICAL_TRUE, element);
       */

      if (aep == 0.0)
	{
	  PetscPrintf (PETSC_COMM_WORLD,
		       "\nError: Problem setting up energy matrix\n");
	  exit (LOGICAL_ERROR);
	}

      sumT = 0.0;

      for (j = 0; j < n; j++)
	{
	  sumT += aen[j] * V_GetCmp (&xT0l, ani[j]);
	}

      V_SetCmp (&xT, elements[element].index, (bep - sumT) / aep);

    }

  VecGhostRestoreLocalForm (xT0, &xT0l);
  VecGhostRestoreLocalForm (xT, &xTl);
	  	  
  VecGhostRestoreLocalForm (dens, &densl);
  VecGhostRestoreLocalForm (visc, &viscl);
  VecGhostRestoreLocalForm (spheat, &spheatl);
  VecGhostRestoreLocalForm (thcond, &thcondl);
    
  VecAssemblyBegin (xT);
  VecAssemblyEnd (xT);

}

void CalculateTemperature (char *var, int *fiter, double dt, double maxCp, int verbose, int pchecks)
{

  int i;

  double mres;
  int miter;
  double mtime;

  double tempc;

  if (parameter.calc[iT] == LOGICAL_FALSE)
    return;
    
  V_Constr (&xTp, nbelements, 0);	// Temperature at cell center - previous iteration

  // Store previous time step values
  VecCopy (xT, xT0);

  fiter[iT]++;

  for (i = 0; i <= parameter.northocor; i++)
    {

      if (parameter.timemethod[iT] == IMPLICITEULER
	  || parameter.timemethod[iT] == CRANKNICOLSON)
	{
	  Q_Constr (&Ae, nbelements, LOGICAL_FALSE);
	  V_Constr (&bT, nbelements, 0);	// Energy source
	}

      // Store previous iteration values
      VecCopy (xT, xTp);

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
		  PetscPrintf (PETSC_COMM_WORLD,
			       "\nWarning: Energy matrix is not diagonal dominant\n");
		  //MatView(Ae, PETSC_VIEWER_STDOUT_WORLD);
		  WriteMatrix (&Ae);
		  WriteVector (&bT);
		  //exit (LOGICAL_ERROR);
		}
	    }

	  // Solve matrix to get temperature T
	  SolveMatrix (&Ae, &xT, &bT, &miter, &mres, &mtime,
		       parameter.msolver[iT], parameter.mprecond[iT],
		       parameter.miter[iT], parameter.mtol[iT]);

	  if (verbose == LOGICAL_TRUE)
	    PetscPrintf (PETSC_COMM_WORLD,
			 "\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
			 var[iT], miter, mres, mtime);

	  if (pchecks == LOGICAL_TRUE)
	    {
	      if (mres > parameter.mtol[iT] && miter == parameter.miter[iT])
		{
		  PetscPrintf (PETSC_COMM_WORLD,
			       "\nError: Problem solving matrix %c\n",
			       var[iT]);
		  exit (LOGICAL_ERROR);
		}
	    }

	}

      // Calculate temperature convergence
      VecWAXPY (temp1, -1.0, xT, xTp);
      VecNorm (temp1, NORM_2, &tempc);

      if (verbose == LOGICAL_TRUE)
	PetscPrintf (PETSC_COMM_WORLD,
		     "\nNon-orthogonality error (energy): %+E\n", tempc);

      CorrectFaceT ();

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
