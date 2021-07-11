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

#include "velocity.h"

void
CorrectVelocityField ()
{

  unsigned int i, j;

  unsigned int face;

  register unsigned int element;

  msh_vector gradp;

  msh_vector sum1;
  msh_vector sum2;

  double u, v, w;

  VecGhostGetLocalForm (xp, &xpl);

  VecGhostGetLocalForm (ap, &apl);
  VecGhostGetLocalForm (hu, &hul);
  VecGhostGetLocalForm (hv, &hvl);
  VecGhostGetLocalForm (hw, &hwl);
    
  if (parameter.calc[is] == LOGICAL_TRUE)
    {

      for (i = 0; i < nbelements; i++)
	{

	  element = i;

	  sum1.x = 0.0;
	  sum1.y = 0.0;
	  sum1.z = 0.0;

	  sum2.x = 0.0;
	  sum2.y = 0.0;
	  sum2.z = 0.0;

	  for (j = 0; j < elements[element].nbfaces; j++)
	    {

	      face = elements[element].face[j];

	      sum1.x += LABS (faces[face].A.x);
	      sum1.y += LABS (faces[face].A.y);
	      sum1.z += LABS (faces[face].A.z);

	      sum2.x += V_GetCmp (&uf, face) * faces[face].A.x;
	      sum2.y += V_GetCmp (&uf, face) * faces[face].A.y;
	      sum2.z += V_GetCmp (&uf, face) * faces[face].A.z;
	      
	    }

	  u = sum2.x / sum1.x;
	  v = sum2.y / sum1.y;
	  w = sum2.z / sum1.z;
	
	  //printf("u: %f\n", u);
	  //printf("v: %f\n", v);
	  //printf("w: %f\n", w);
	  
	  V_SetCmp (&xu, elements[element].index, u);
	  V_SetCmp (&xv, elements[element].index, v);
	  V_SetCmp (&xw, elements[element].index, w);

	}

    }
  else
    {

      for (i = 0; i < nbelements; i++)
	{

	  element = i;

	  gradp = Gradient (&xpl, &xpf, LOGICAL_TRUE, element);

	  V_SetCmp (&xu, elements[element].index, (V_GetCmp (&hul, element) - gradp.x) / V_GetCmp (&apl, element));
	  V_SetCmp (&xv, elements[element].index, (V_GetCmp (&hvl, element) - gradp.y) / V_GetCmp (&apl, element));
	  V_SetCmp (&xw, elements[element].index, (V_GetCmp (&hwl, element) - gradp.z) / V_GetCmp (&apl, element));

	}

    }

  VecGhostRestoreLocalForm (xp, &xpl);
        
  VecGhostRestoreLocalForm (ap, &apl);
  VecGhostRestoreLocalForm (hu, &hul);
  VecGhostRestoreLocalForm (hv, &hvl);
  VecGhostRestoreLocalForm (hw, &hwl);
    
  VecAssemblyBegin (xu);
  VecAssemblyEnd (xu);
  VecAssemblyBegin (xv);
  VecAssemblyEnd (xv);
  VecAssemblyBegin (xw);
  VecAssemblyEnd (xw);

}

void
CorrectFaceUVW ()
{

  int i;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  msh_element ghost;

  double apj;

  //double dNf, dPf;
  double lambda;

  double dj;

  double ppl;
  double pnl;

  msh_vector gradpp;
  msh_vector gradpn;

  VecGhostGetLocalForm (xu, &xul);
  VecGhostGetLocalForm (xv, &xvl);
  VecGhostGetLocalForm (xw, &xwl);
  VecGhostGetLocalForm (xp, &xpl);
  VecGhostGetLocalForm (xs, &xsl);
  VecGhostGetLocalForm (ap, &apl);

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
	  dNf = GeoMagVector (GeoSubVectorVector (elements[neighbor].celement, faces[face].cface));
	  dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement, faces[face].cface));

	  lambda = dPf / (dPf + dNf);
	  */

	  lambda = 0.5;

	  apj = V_GetCmp (&apl, neighbor) * lambda + V_GetCmp (&apl, element) * (1.0 - lambda);

	  if (parameter.orthof != 0.0)
	    gradpn = Gradient (&xpl, &xpf, LOGICAL_TRUE, neighbor);

	  ppl = V_GetCmp (&xpl, element);

	  if (parameter.orthof != 0.0)
	    {
	      ppl += parameter.orthof * GeoDotVectorVector (gradpp, GeoSubVectorVector (faces[face].rpl, elements[element].celement));
	    }

	  pnl = V_GetCmp (&xpl, neighbor);

	  if (parameter.orthof != 0.0)
	    {

	      pnl += parameter.orthof * GeoDotVectorVector (gradpn, GeoSubVectorVector (faces[face].rnl, elements[neighbor].celement));
	    }

	  V_SetCmp (&uf, face, V_GetCmp (&uf, face) - 1.0 / (apj * (faces[face].dj + faces[face].kj)) * (pnl - ppl));

	  V_SetCmp (&xuf, face, V_GetCmp (&uf, face) * faces[face].n.x);
	  V_SetCmp (&xvf, face, V_GetCmp (&uf, face) * faces[face].n.y);
	  V_SetCmp (&xwf, face, V_GetCmp (&uf, face) * faces[face].n.z);

	  V_SetCmp (&xpf, face, V_GetCmp (&xpl, neighbor) * lambda + V_GetCmp (&xpl, element) * (1.0 - lambda));

	}
      else
	{

	  if (faces[face].bc == PROCESSOR)
	    {

	      ghost.index = faces[face].physreg;

	      ghost.celement.x = V_GetCmp (&cexl, faces[face].ghost);
	      ghost.celement.y = V_GetCmp (&ceyl, faces[face].ghost);
	      ghost.celement.z = V_GetCmp (&cezl, faces[face].ghost);

	      dj = GeoMagVector (GeoSubVectorVector (ghost.celement, elements[element].celement));

	      /*
	      dNf = GeoMagVector (GeoSubVectorVector (ghost.celement, faces[face].cface));
	      dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement, faces[face].cface));

	      lambda = dPf / (dPf + dNf);
	      */

	      lambda = 0.5;

	      apj = V_GetCmp (&apl, faces[face].ghost) * lambda + V_GetCmp (&apl, element) * (1.0 - lambda);

	      /*
	         gradpn = Gradient (&xpl, &xpf, LOGICAL_TRUE, faces[face].ghost);

	         ppl =
	         V_GetCmp (&xpl, element) 
	         + parameter.orthof * GeoDotVectorVector (gradpp,
	         GeoSubVectorVector
	         (faces
	         [face].
	         rpl,
	         elements
	         [element].
	         celement));

	         pnl =
	         V_GetCmp (&xpl, faces[face].ghost) 
	         + parameter.orthof * GeoDotVectorVector (gradpn,                                     
	         GeoSubVectorVector
	         (faces
	         [face].
	         rnl,
	         ghost.
	         celement));
	       */

	      ppl = V_GetCmp (&xpl, element);

	      pnl = V_GetCmp (&xpl, faces[face].ghost);

	      V_SetCmp (&uf, face, V_GetCmp (&uf, face) - 1.0 / (apj * dj) * (pnl - ppl));

	      V_SetCmp (&xuf, face, V_GetCmp (&uf, face) * faces[face].n.x);
	      V_SetCmp (&xvf, face, V_GetCmp (&uf, face) * faces[face].n.y);
	      V_SetCmp (&xwf, face, V_GetCmp (&uf, face) * faces[face].n.z);

	      V_SetCmp (&xpf, face, V_GetCmp (&xpl, faces[face].ghost) * lambda + V_GetCmp (&xpl, element) * (1.0 - lambda));

	    }
	  else
	    {

	      apj = V_GetCmp (&apl, element);

	      ppl = V_GetCmp (&xpl, element);

	      if (parameter.orthof != 0.0)
		{

		  ppl += parameter.orthof * GeoDotVectorVector (gradpp, GeoSubVectorVector (faces[face].rpl, elements[element].celement));
		}

	      if (faces[face].bc == PERMEABLE)
		{

		  V_SetCmp (&uf, face, (V_GetCmp (&uf, face) - 1.0 / (apj * (faces[face].dj + faces[face].kj)) * (V_GetCmp (&xpf, face) - ppl)));

		  V_SetCmp (&xuf, face, V_GetCmp (&uf, face) * faces[face].n.x * (1.0 - V_GetCmp(&xsl, element)));
		  V_SetCmp (&xvf, face, V_GetCmp (&uf, face) * faces[face].n.y * (1.0 - V_GetCmp(&xsl, element)));
		  V_SetCmp (&xwf, face, V_GetCmp (&uf, face) * faces[face].n.z * (1.0 - V_GetCmp(&xsl, element)));

		}

	      if (faces[face].bc == OUTLET)
		{

		  // velocity gradient = 0
		  // specified pressure

		  V_SetCmp (&uf, face, V_GetCmp (&uf, face) - 1.0 / (apj * (faces[face].dj + faces[face].kj)) * (V_GetCmp (&xpf, face) - ppl));

		  V_SetCmp (&xuf, face, V_GetCmp (&uf, face) * faces[face].n.x);
		  V_SetCmp (&xvf, face, V_GetCmp (&uf, face) * faces[face].n.y);
		  V_SetCmp (&xwf, face, V_GetCmp (&uf, face) * faces[face].n.z);

		}

	      if (faces[face].bc == PRESSUREINLET)
		{

		  // velocity gradient = 0
		  // specified pressure

		  V_SetCmp (&uf, face, V_GetCmp (&uf, face) - 1.0 / (apj * (faces[face].dj + faces[face].kj)) * (V_GetCmp (&xpf, face) - ppl));

		  V_SetCmp (&xuf, face, V_GetCmp (&uf, face) * faces[face].n.x);
		  V_SetCmp (&xvf, face, V_GetCmp (&uf, face) * faces[face].n.y);
		  V_SetCmp (&xwf, face, V_GetCmp (&uf, face) * faces[face].n.z);

		}

	      if (faces[face].bc == PROCESSOR)
		{

		  // specified pressure
		  // specified velocity

		  V_SetCmp (&uf, face, V_GetCmp (&uf, face) - 1.0 / (apj * (faces[face].dj + faces[face].kj)) * (V_GetCmp (&xpf, face) - ppl));

		}

	      if (faces[face].bc == INLET ||
		  faces[face].bc == MOVINGWALL ||
		  faces[face].bc == WALL ||
		  faces[face].bc == ADIABATICWALL
		  || faces[face].bc == SURFACE)
		{

		  // pressure gradient = 0
		  // specified velocity

		  V_SetCmp (&uf, face, V_GetCmp (&xuf, face) * faces[face].n.x + V_GetCmp (&xvf, face) * faces[face].n.y + V_GetCmp (&xwf, face) * faces[face].n.z);

		  V_SetCmp (&xpf, face, ppl);

		}

	      if (faces[face].bc == SLIP)
		{

		  // pressure gradient = 0
		  // velocity gradient = 0

		  V_SetCmp (&xuf, face, V_GetCmp (&xul, element));
		  V_SetCmp (&xvf, face, V_GetCmp (&xvl, element));
		  V_SetCmp (&xwf, face, V_GetCmp (&xwl, element));

		  V_SetCmp (&uf, face, 0.0);

		  V_SetCmp (&xpf, face, ppl);

		}
	    }
	}

    }

  VecGhostRestoreLocalForm (xu, &xul);
  VecGhostRestoreLocalForm (xv, &xvl);
  VecGhostRestoreLocalForm (xw, &xwl);
  VecGhostRestoreLocalForm (xp, &xpl);
  VecGhostRestoreLocalForm (xs, &xsl);
  VecGhostRestoreLocalForm (ap, &apl);
    
  VecAssemblyBegin (xuf);
  VecAssemblyEnd (xuf);
  VecAssemblyBegin (xvf);
  VecAssemblyEnd (xvf);
  VecAssemblyBegin (xwf);
  VecAssemblyEnd (xwf);
  VecAssemblyBegin (uf);
  VecAssemblyEnd (uf);

}

void
CalculateCorrectionFactors ()
{

  VecGhostUpdateBegin (ap, INSERT_VALUES, SCATTER_FORWARD);

  MatGetDiagonal (Am, temp2);
  VecPointwiseMult (temp1, temp2, xu);
  VecAXPY (hu, +1.0, temp1);
  MatMult (Am, xu, temp2);
  VecAXPY (hu, -1.0, temp2);

  MatGetDiagonal (Am, temp2);
  VecPointwiseMult (temp1, temp2, xv);
  VecAXPY (hv, +1.0, temp1);
  MatMult (Am, xv, temp2);
  VecAXPY (hv, -1.0, temp2);

  MatGetDiagonal (Am, temp2);
  VecPointwiseMult (temp1, temp2, xw);
  VecAXPY (hw, +1.0, temp1);
  MatMult (Am, xw, temp2);
  VecAXPY (hw, -1.0, temp2);
  
  VecGhostUpdateBegin (hu, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateBegin (hv, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateBegin (hw, INSERT_VALUES, SCATTER_FORWARD);
  
  VecGhostUpdateEnd (hu, INSERT_VALUES, SCATTER_FORWARD);  
  VecGhostUpdateEnd (hv, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (hw, INSERT_VALUES, SCATTER_FORWARD);

  VecGhostUpdateEnd (ap, INSERT_VALUES, SCATTER_FORWARD);
    
}

void
BuildMomentumMatrix (double dt)
{

  unsigned int i, j, n;

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  double app;

  double apn[MAXFACES];
  unsigned int ani[MAXFACES];

  double bpu, bpv, bpw;

  double densp;
  double viscj;

  msh_element ghost;

  //msh_vector gradup, gradvp, gradwp;
  //msh_vector gradun, gradvn, gradwn;
  //msh_vector gradvisc;

  //msh_vector gradp;

  //double dNf, dPf;
  double lambda;
  double xsi;

  double dj;

  msh_vector g;

  // MatSetValues
  int row;
  int ncols;
  int col[MAXFACES+1];
  int nvals;
  double val[MAXFACES+1];

  g.x = parameter.g[0];
  g.y = parameter.g[1];
  g.z = parameter.g[2];

  VecGhostGetLocalForm (xu0, &xu0l);
  VecGhostGetLocalForm (xv0, &xv0l);
  VecGhostGetLocalForm (xw0, &xw0l);
  
  VecGhostGetLocalForm (xu, &xul);
  VecGhostGetLocalForm (xv, &xvl);
  VecGhostGetLocalForm (xw, &xwl);
  VecGhostGetLocalForm (xp, &xpl);
  VecGhostGetLocalForm (xs, &xsl);
  
  VecGhostGetLocalForm (dens, &densl);
  VecGhostGetLocalForm (visc, &viscl);
  
  // Equation: dU/dt + div(rho*U*U) - div(mi*grad(U)) = qU

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      bpu = 0.0;
      bpv = 0.0;
      bpw = 0.0;

      app = 0.0;

      n = 0;

      /*
         gradup = Gradient (&xul, &xuf, LOGICAL_TRUE, element);
         gradvp = Gradient (&xvl, &xvf, LOGICAL_TRUE, element);
         gradwp = Gradient (&xwl, &xwf, LOGICAL_TRUE, element);
       */

      densp = V_GetCmp (&densl, element);

      for (j = 0; j < elements[element].nbfaces; j++)
	{

	  face = elements[element].face[j];

	  pair = faces[face].pair;

	  if (pair != -1)
	    {

	      neighbor = faces[pair].element;

	      /*
	      dNf = GeoMagVector (GeoSubVectorVector(elements[neighbor].celement, faces[face].cface));
	      dPf = GeoMagVector (GeoSubVectorVector(elements[element].celement,  faces[face].cface));

	      lambda = dPf / (dPf + dNf);
	      */

	      lambda = 0.5;

	      viscj = V_GetCmp (&viscl, element) * (1.0 - lambda) + V_GetCmp (&viscl, neighbor) * lambda;

	      // Convection           
	      if (parameter.scheme[iu] == UDS)
		{
		  // UDS
		  if (V_GetCmp (&uf, face) > 0.0)
		    xsi = 0.0;
		  else
		    xsi = 1.0;

		}
	      else
		{
		  // CDS
		  xsi = lambda;
		}

	      // Convection
	      app += (1.0 - xsi) * densp * V_GetCmp (&uf, face) * faces[face].Aj / elements[element].Vp;

	      // Diffusion
	      app += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;

	      // Convection
	      apn[n] = xsi * densp * V_GetCmp (&uf, face) * faces[face].Aj / elements[element].Vp;

	      // Diffusion
	      apn[n] += -viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;

	      ani[n] = elements[neighbor].index;
	      n++;

	      /*
	         // Non-orthogonal correction terms
	         if (parameter.orthof != 0.0)
	         {
	         gradun = Gradient (&xul, &xuf, LOGICAL_TRUE, neighbor);
	         gradvn = Gradient (&xvl, &xvf, LOGICAL_TRUE, neighbor);
	         gradwn = Gradient (&xwl, &xwf, LOGICAL_TRUE, neighbor);

	         bpu += parameter.orthof * 
	         -viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
	         elements[element].Vp *
	         (GeoDotVectorVector
	         (gradun,
	         GeoSubVectorVector (faces[face].rnl,
	         elements[neighbor].celement)) -
	         GeoDotVectorVector (gradup,
	         GeoSubVectorVector (faces[face].rpl,
	         elements[element].
	         celement)));

	         bpv += parameter.orthof * 
	         -viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
	         elements[element].Vp *
	         (GeoDotVectorVector
	         (gradvn,
	         GeoSubVectorVector (faces[face].rnl,
	         elements[neighbor].celement)) -
	         GeoDotVectorVector (gradvp,
	         GeoSubVectorVector (faces[face].rpl,
	         elements[element].
	         celement)));

	         bpw += parameter.orthof * 
	         -viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
	         elements[element].Vp *
	         (GeoDotVectorVector
	         (gradwn,
	         GeoSubVectorVector (faces[face].rnl,
	         elements[neighbor].celement)) -
	         GeoDotVectorVector (gradwp,
	         GeoSubVectorVector (faces[face].rpl,
	         elements[element].
	         celement)));        
	         }
	       */

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

		  dj = GeoMagVector (GeoSubVectorVector
				     (ghost.celement,
				      elements[element].celement));

		  /*
		  dNf = GeoMagVector (GeoSubVectorVector (ghost.celement, faces[face].cface));
		  dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement, faces[face].cface));

		  lambda = dPf / (dPf + dNf);
		  */

		  lambda = 0.5;

		  viscj = V_GetCmp (&viscl, element) * (1.0 - lambda) + V_GetCmp (&viscl, faces[face].ghost) * lambda;

		  if (parameter.scheme[iu] == UDS)
		    {
		      // UDS
		      if (V_GetCmp (&uf, face) > 0.0)
			xsi = 0.0;
		      else
			xsi = 1.0;

		    }
		  else
		    {
		      // CDS
		      xsi = lambda;
		    }

		  // Convection
		  app += (1.0 - xsi) * densp * V_GetCmp (&uf, face) * faces[face].Aj / elements[element].Vp;

		  // Diffusion
		  app += viscj * faces[face].Aj / dj / elements[element].Vp;

		  // Convection
		  apn[n] = xsi * densp * V_GetCmp (&uf, face) * faces[face].Aj / elements[element].Vp;

		  // Diffusion
		  apn[n] += -viscj * faces[face].Aj / dj / elements[element].Vp;

		  ani[n] = ghost.index;
		  n++;

		  /*
		     // Non-orthogonal correction terms
		     if (parameter.orthof != 0.0)
		     {
		     gradun = Gradient (&xul, &xuf, LOGICAL_TRUE, faces[face].ghost);    
		     gradvn = Gradient (&xvl, &xvf, LOGICAL_TRUE, faces[face].ghost);
		     gradwn = Gradient (&xwl, &xwf, LOGICAL_TRUE, faces[face].ghost);

		     bpu += parameter.orthof * 
		     -viscj * faces[face].Aj / dj /
		     elements[element].Vp *
		     (GeoDotVectorVector
		     (gradun,
		     GeoSubVectorVector (faces[face].rnl,
		     ghost.celement)) -
		     GeoDotVectorVector (gradup,
		     GeoSubVectorVector (faces[face].rpl,
		     elements
		     [element].
		     celement)));

		     bpv += parameter.orthof * 
		     -viscj * faces[face].Aj / dj /
		     elements[element].Vp *
		     (GeoDotVectorVector
		     (gradvn,
		     GeoSubVectorVector (faces[face].rnl,
		     ghost.celement)) -
		     GeoDotVectorVector (gradvp,
		     GeoSubVectorVector (faces[face].rpl,
		     elements
		     [element].
		     celement)));

		     bpw += parameter.orthof * 
		     -viscj * faces[face].Aj / dj /
		     elements[element].Vp *
		     (GeoDotVectorVector
		     (gradwn,
		     GeoSubVectorVector (faces[face].rnl,
		     ghost.celement)) -
		     GeoDotVectorVector (gradwp,
		     GeoSubVectorVector (faces[face].rpl,
		     elements
		     [element].
		     celement)));
		     }
		   */

		}
	      else
		{

		  if (faces[face].bc != EMPTY)
		    {

		      viscj = V_GetCmp (&viscl, element);

		      if ((faces[face].bc != PERMEABLE || V_GetCmp (&xsl, element) > 0.5) && faces[face].bc != SLIP)
		        {
		           // Diffusion
		           app += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;

		           bpu += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) * V_GetCmp (&xuf, face) / elements[element].Vp;
		           bpv += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) * V_GetCmp (&xvf, face) / elements[element].Vp;
		           bpw += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) * V_GetCmp (&xwf, face) / elements[element].Vp;
 		        }

		      // Convection
		      if (parameter.scheme[iu] == UDS)
			{

			  // UDS
			  if (V_GetCmp (&uf, face) > 0.0)
			    {

			      // Convection
			      app += densp * V_GetCmp (&uf, face) * faces[face].Aj / elements[element].Vp;

			    }
			  else
			    {

			      // Convection
			      bpu += -densp * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xuf, face) / elements[element].Vp;
			      bpv += -densp * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xvf, face) / elements[element].Vp;
			      bpw += -densp * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xwf, face) / elements[element].Vp;

			    }

			}
		      else
			{

			  // CDS
			  bpu += -densp * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xuf, face) / elements[element].Vp;
			  bpv += -densp * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xvf, face) / elements[element].Vp;
			  bpw += -densp * V_GetCmp (&uf, face) * faces[face].Aj * V_GetCmp (&xwf, face) / elements[element].Vp;

			}

		      /*
		         // Non-orthogonal correction terms
		         if (parameter.orthof != 0.0)
		         {
		         bpu += parameter.orthof * 
		         viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
		         elements[element].Vp * GeoDotVectorVector (gradup,
		         GeoSubVectorVector
		         (faces[face].
		         rpl,
		         elements
		         [element].
		         celement));

		         bpv += parameter.orthof * 
		         viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
		         elements[element].Vp * GeoDotVectorVector (gradvp,
		         GeoSubVectorVector
		         (faces[face].
		         rpl,
		         elements
		         [element].
		         celement));

		         bpw += parameter.orthof * 
		         viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
		         elements[element].Vp * GeoDotVectorVector (gradwp,
		         GeoSubVectorVector
		         (faces[face].
		         rpl,
		         elements
		         [element].               
		         celement));

		         }
		       */
		    }
		}
	    }
	}

      if (dt > 0)
	{

	  // Unsteady term - Euler

	  app += densp / dt * parameter.st;

	  bpu += densp / dt * V_GetCmp (&xu0l, element);
	  bpv += densp / dt * V_GetCmp (&xv0l, element);
	  bpw += densp / dt * V_GetCmp (&xw0l, element);

	}

      // Source - viscous term
      /*
         gradvisc = Gradient (&viscl, NULL, LOGICAL_FALSE, element);

         bpu += GeoDotVectorVector (gradup, gradvisc);
         bpv += GeoDotVectorVector (gradvp, gradvisc);
         bpw += GeoDotVectorVector (gradwp, gradvisc);
       */

      // Source - gravity
      bpu += densp * g.x;
      bpv += densp * g.y;
      bpw += densp * g.z;

      // Initialize H with source contribution without pressure 
      V_SetCmp (&hu, elements[element].index, bpu);
      V_SetCmp (&hv, elements[element].index, bpv);
      V_SetCmp (&hw, elements[element].index, bpw);

      // Source - pressure 
      /*
         gradp = Gradient (&xpl, &xpf, LOGICAL_TRUE, element);

         bpu += -gradp.x;
         bpv += -gradp.y;
         bpw += -gradp.z;
       */

      if (app == 0.0 || app != app)
	{
	  PetscPrintf (PETSC_COMM_WORLD, "\nError: Problem setting up momentum matrix\n");
	  exit (LOGICAL_ERROR);
	}

      V_SetCmp (&ap, elements[element].index, app);

      /*
      Q_SetEntry (&Am, elements[element].index, elements[element].index, app);

      for (j = 0; j < n; j++)
	{
	  Q_SetEntry (&Am, elements[element].index, ani[j], apn[j]);
	}
      */

      ncols = 0;
      nvals = 0;

      row = elements[element].index;

      col[ncols] = elements[element].index;
      ncols++;

      val[nvals] = app;
      nvals++;

      for (j = 0; j < n; j++)
	{
		col[ncols] = ani[j];
	        ncols++;

		val[nvals] = apn[j];
		nvals++;
	}

      Q_SetEntries (&Am, 1, &row, ncols, col, val);

      if (parameter.calc[iu] == LOGICAL_TRUE)
	V_SetCmp (&bu, elements[element].index, bpu);

      if (parameter.calc[iv] == LOGICAL_TRUE)
	V_SetCmp (&bv, elements[element].index, bpv);

      if (parameter.calc[iw] == LOGICAL_TRUE)
	V_SetCmp (&bw, elements[element].index, bpw);

    }

  VecGhostRestoreLocalForm (xu0, &xu0l);
  VecGhostRestoreLocalForm (xv0, &xv0l);
  VecGhostRestoreLocalForm (xw0, &xw0l);

  VecGhostRestoreLocalForm (xu, &xul);
  VecGhostRestoreLocalForm (xv, &xvl);
  VecGhostRestoreLocalForm (xw, &xwl);
  VecGhostRestoreLocalForm (xp, &xpl);
  VecGhostRestoreLocalForm (xs, &xsl);

  VecGhostRestoreLocalForm (dens, &densl);
  VecGhostRestoreLocalForm (visc, &viscl);

  MatAssemblyBegin (Am, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd (Am, MAT_FINAL_ASSEMBLY);

  if (parameter.calc[iu] == LOGICAL_TRUE)
    {
      VecAssemblyBegin (bu);
      VecAssemblyEnd (bu);
    }

  if (parameter.calc[iv] == LOGICAL_TRUE)
    {
      VecAssemblyBegin (bv);
      VecAssemblyEnd (bv);
    }

  if (parameter.calc[iw] == LOGICAL_TRUE)
    {
      VecAssemblyBegin (bw);
      VecAssemblyEnd (bw);
    }

  VecAssemblyBegin (ap);
  VecAssemblyEnd (ap);
  VecAssemblyBegin (hu);
  VecAssemblyEnd (hu);
  VecAssemblyBegin (hv);
  VecAssemblyEnd (hv);
  VecAssemblyBegin (hw);
  VecAssemblyEnd (hw);

}

void CorrectVelocity (char *var, int *fiter, double dt, double maxCp, int verbose, int pchecks)
{

  if (parameter.calc[ip] == LOGICAL_TRUE)
    {
  
      // Correct face values      
      CorrectFaceUVW ();

      // Correct cell center   
      CorrectVelocityField ();

  }
  
}

void CalculateVelocity (char *var, int *fiter, double dt, double maxCp, int verbose, int pchecks)
{

  double mres;
  int miter;
  double mtime;
      
  Q_Constr (&Am, nbelements, LOGICAL_FALSE);
  V_Constr (&bu, nbelements, 0);	// Momentum source x-component
  V_Constr (&bv, nbelements, 0);	// Momentum source y-component
  V_Constr (&bw, nbelements, 0);	// Momentum source z-component

  // Store previous time step values
  VecCopy (xu, xu0);
  VecCopy (xv, xv0);
  VecCopy (xw, xw0);

  // Build three momentum matrices for u, v, w velocity components
  if (parameter.calc[ip] == LOGICAL_TRUE)
    {
      BuildMomentumMatrix (dt);
    }

  if (pchecks == LOGICAL_TRUE)
    {
      if (!CheckIfDiagonalMatrix (&Am))
	{
	  PetscPrintf (PETSC_COMM_WORLD, "\nWarning: Momentum matrix is not diagonal dominant\n");
	  //MatView(Am, PETSC_VIEWER_STDOUT_WORLD);       
	  WriteMatrix (&Am);
	  WriteVector (&xu);
	  WriteVector (&xv);
	  WriteVector (&xw);
	  //exit (LOGICAL_ERROR);
	}
    }

  if (parameter.calc[iu] == LOGICAL_TRUE)
    {

      fiter[iu]++;

      // Solve matrix for u velocity component
      SolveMatrix (&Am, &xu, &bu, &miter, &mres, &mtime, parameter.msolver[iu], parameter.mprecond[iu], parameter.miter[iu], parameter.mtol[iu]);

      if (verbose == LOGICAL_TRUE)
	PetscPrintf (PETSC_COMM_WORLD, "\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n", var[iu], miter, mres, mtime);

      if (pchecks == LOGICAL_TRUE)
	{
	  if (mres > parameter.mtol[iu] && miter == parameter.miter[iu])
	    {
	      PetscPrintf (PETSC_COMM_WORLD,
			   "\nError: Problem solving matrix %c\n", var[iu]);
	      exit (LOGICAL_ERROR);
	    }
	}
    }

  if (parameter.calc[iv] == LOGICAL_TRUE)
    {

      fiter[iv]++;

      // Solve matrix for v velocity component
      SolveMatrix (&Am, &xv, &bv, &miter, &mres, &mtime, parameter.msolver[iv], parameter.mprecond[iv], parameter.miter[iv], parameter.mtol[iv]);

      if (verbose == LOGICAL_TRUE)
	PetscPrintf (PETSC_COMM_WORLD, "\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n", var[iv], miter, mres, mtime);

      if (pchecks == LOGICAL_TRUE)
	{
	  if (mres > parameter.mtol[iv] && miter == parameter.miter[iv])
	    {
	      PetscPrintf (PETSC_COMM_WORLD,
			   "\nError: Problem solving matrix %c\n", var[iv]);
	      exit (LOGICAL_ERROR);
	    }
	}

    }

  if (parameter.calc[iw] == LOGICAL_TRUE)
    {

      fiter[iw]++;

      // Solve matrix for w velocity component
      SolveMatrix (&Am, &xw, &bw, &miter, &mres, &mtime, parameter.msolver[iw], parameter.mprecond[iw], parameter.miter[iw], parameter.mtol[iw]);

      if (verbose == LOGICAL_TRUE)
	PetscPrintf (PETSC_COMM_WORLD, "\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n", var[iw], miter, mres, mtime);

      if (pchecks == LOGICAL_TRUE)
	{
	  if (mres > parameter.mtol[iw] && miter == parameter.miter[iw])
	    {
	      PetscPrintf (PETSC_COMM_WORLD, "\nProblem solving matrix %c\n", var[iw]);
	      exit (LOGICAL_ERROR);
	    }
	}
    }

  // Calculate correction factors
  if (parameter.calc[ip] == LOGICAL_TRUE)
    {
      CalculateCorrectionFactors ();
    }

  Q_Destr (&Am);

  V_Destr (&bu);
  V_Destr (&bv);
  V_Destr (&bw);
  
}
