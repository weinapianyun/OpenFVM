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

#include "velocity.h"

void
CorrectVelocityField ()
{

  unsigned int i, j;

  register unsigned int face;
  register unsigned int element;

  msh_vector gradp;

  msh_vector sum1;
  msh_vector sum2;

  double u, v, w;

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

	      sum2.x += V_GetCmp (&uf, face + 1) * faces[face].A.x;
	      sum2.y += V_GetCmp (&uf, face + 1) * faces[face].A.y;
	      sum2.z += V_GetCmp (&uf, face + 1) * faces[face].A.z;

	    }

	  u = sum2.x / sum1.x;
	  v = sum2.y / sum1.y;
	  w = sum2.z / sum1.z;

	  //printf("u: %f\n", u);
	  //printf("v: %f\n", v);
	  //printf("w: %f\n", w);

	  V_SetCmp (&xu, element + 1, u);
	  V_SetCmp (&xv, element + 1, v);
	  V_SetCmp (&xw, element + 1, w);

	}

    }
  else
    {

      for (i = 0; i < nbelements; i++)
	{

	  element = i;

	  gradp = Gradient (&xp, &xpf, LOGICAL_TRUE, element);

	  V_SetCmp (&xu, element + 1, (V_GetCmp (&hu, element + 1) - gradp.x) / V_GetCmp (&ap, element + 1));
	  V_SetCmp (&xv, element + 1, (V_GetCmp (&hv, element + 1) - gradp.y) / V_GetCmp (&ap, element + 1));
	  V_SetCmp (&xw, element + 1, (V_GetCmp (&hw, element + 1) - gradp.z) / V_GetCmp (&ap, element + 1));

	}

    }

}

void
CorrectFaceUVW ()
{

  int i;

  unsigned int face, pair;

  register unsigned int element;
  unsigned int neighbor;

  double apj;

  //double dNf, dPf;
  double lambda;

  double ppl;
  double pnl;

  msh_vector gradpp;
  msh_vector gradpn;

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

	  /*
	  dNf = GeoMagVector (GeoSubVectorVector (elements[neighbor].celement, faces[face].cface));
	  dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement, faces[face].cface));

	  lambda = dPf / (dPf + dNf);
	  */

	  lambda = 0.5;

	  apj = V_GetCmp (&ap, neighbor + 1) * lambda + V_GetCmp (&ap, element + 1) * (1.0 - lambda);

	  if (parameter.orthof != 0.0)
	    gradpn = Gradient (&xp, &xpf, LOGICAL_TRUE, neighbor);

	  ppl = V_GetCmp (&xp, element + 1);

	  if (parameter.orthof != 0.0)
	    {
	      ppl += parameter.orthof * GeoDotVectorVector (gradpp, GeoSubVectorVector (faces[face].rpl, elements[element].celement));
	    }

	  pnl = V_GetCmp (&xp, neighbor + 1);

	  if (parameter.orthof != 0.0)
	    {
	      pnl += parameter.orthof * GeoDotVectorVector (gradpn, GeoSubVectorVector (faces[face].rnl, elements[neighbor].celement));
	    }

	  V_SetCmp (&uf, face + 1, V_GetCmp (&uf, face + 1) - 1.0 / apj * (pnl - ppl) / (faces[face].dj + faces[face].kj));

	  V_SetCmp (&xuf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.x);
	  V_SetCmp (&xvf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.y);
	  V_SetCmp (&xwf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.z);

	}
      else
	{

	  apj = V_GetCmp (&ap, element + 1);

	  ppl = V_GetCmp (&xp, element + 1);

	  if (parameter.orthof != 0.0)
	    {
	      ppl += parameter.orthof * GeoDotVectorVector (gradpp, GeoSubVectorVector (faces[face].rpl, elements[element].celement));
	    }

	  if (faces[face].bc == PERMEABLE)
	    {

	      V_SetCmp (&uf, face + 1, (V_GetCmp (&uf, face + 1) - 1.0 / (apj * (faces[face].dj + faces[face].kj)) * (V_GetCmp (&xpf, face + 1) - ppl)));

	      // Wall
	      /*
	      V_SetCmp (&xuf, face + 1, V_GetCmp (&uf,  face + 1) * faces[face].n.x * (1.0 - V_GetCmp(&xs, element + 1)));
	      V_SetCmp (&xvf, face + 1, V_GetCmp (&uf,  face + 1) * faces[face].n.y * (1.0 - V_GetCmp(&xs, element + 1)));
	      V_SetCmp (&xwf, face + 1, V_GetCmp (&uf,  face + 1) * faces[face].n.z * (1.0 - V_GetCmp(&xs, element + 1)));
	      */
		
	      // Slip
	      V_SetCmp (&xuf, face + 1, V_GetCmp (&xu, element + 1) * (1.0 - V_GetCmp(&xs, element + 1)));
	      V_SetCmp (&xvf, face + 1, V_GetCmp (&xv, element + 1) * (1.0 - V_GetCmp(&xs, element + 1)));
	      V_SetCmp (&xwf, face + 1, V_GetCmp (&xw, element + 1) * (1.0 - V_GetCmp(&xs, element + 1)));

	    }

	  if (faces[face].bc == OUTLET)
	    {

	      // velocity gradient = 0
	      // specified pressure

	      V_SetCmp (&uf, face + 1, V_GetCmp (&uf, face + 1) - 1.0 / (apj * (faces[face].dj + faces[face].kj)) * (V_GetCmp (&xpf, face + 1) - ppl));

	      V_SetCmp (&xuf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.x); 
	      V_SetCmp (&xvf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.y);
	      V_SetCmp (&xwf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.z);

	    }

	  if (faces[face].bc == PRESSUREINLET)
	    {

	      // velocity gradient = 0
	      // specified pressure

	      V_SetCmp (&uf, face + 1, V_GetCmp (&uf, face + 1) - 1.0 / (apj * (faces[face].dj + faces[face].kj)) * (V_GetCmp (&xpf, face + 1) - ppl));

	      V_SetCmp (&xuf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.x);
	      V_SetCmp (&xvf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.y);
	      V_SetCmp (&xwf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.z);

	    }

	  if (faces[face].bc == INLET ||
	      faces[face].bc == MOVINGWALL ||
	      faces[face].bc == WALL ||
	      faces[face].bc == ADIABATICWALL || faces[face].bc == SURFACE)
	    {

	      // pressure gradient = 0
	      // specified velocity

	      V_SetCmp (&uf, face + 1, V_GetCmp (&xuf, face + 1) * faces[face].n.x + V_GetCmp (&xvf, face + 1) * faces[face].n.y + V_GetCmp (&xwf, face + 1) * faces[face].n.z);

              V_SetCmp (&xpf, face + 1, ppl);

	    }

	  if (faces[face].bc == SLIP)
	    {

	      // pressure gradient = 0
	      // velocity gradient = 0

	      V_SetCmp (&xuf, face + 1, V_GetCmp (&xu, element + 1));
	      V_SetCmp (&xvf, face + 1, V_GetCmp (&xv, element + 1));
	      V_SetCmp (&xwf, face + 1, V_GetCmp (&xw, element + 1));

	      V_SetCmp (&uf, face + 1, 0.0);

	    }

	}

    }

}

void
CalculateCorrectionFactors ()
{

  AddAsgn_VV (&hu, Mul_QV (Sub_QQ (Diag_Q (&Am), &Am), &xu));
  AddAsgn_VV (&hv, Mul_QV (Sub_QQ (Diag_Q (&Am), &Am), &xv));
  AddAsgn_VV (&hw, Mul_QV (Sub_QQ (Diag_Q (&Am), &Am), &xw));

  //PrintVector(&xu);
  //PrintVector(&hu);

}

void
BuildMomentumMatrix (double dt)
{

  int i, j, n;

  unsigned int face, pair;
 
  register unsigned int element;
  unsigned int neighbor;

  double app;

  double apn[MAXFACES];
  unsigned int ani[MAXFACES];

  double bpu, bpv, bpw;

  double densp;
  double viscj;

  //msh_vector gradup, gradvp, gradwp;
  //msh_vector gradun, gradvn, gradwn;
  //msh_vector gradvisc;

  msh_vector gradp;

  //double dNf, dPf;
  double lambda;
  double xsi;

  msh_vector g;

  g.x = parameter.g[0];
  g.y = parameter.g[1];
  g.z = parameter.g[2];

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
         gradup = Gradient (&xu, &xuf, LOGICAL_TRUE, element);
         gradvp = Gradient (&xv, &xvf, LOGICAL_TRUE, element);
         gradwp = Gradient (&xw, &xwf, LOGICAL_TRUE, element);
       */

      densp = V_GetCmp (&dens, element + 1);

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

	      viscj = V_GetCmp (&visc, element + 1) * (1.0 - lambda) + V_GetCmp (&visc, neighbor + 1) * lambda;

	      // Convection
	      if (parameter.scheme[iu] == UDS)
		{

		  // UDS
		  if (V_GetCmp (&uf, face + 1) > 0.0)
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
	      app += (1.0 - xsi) * densp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp;

	      // Diffusion
	      app += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;

	      // Convection
	      apn[n] = xsi * densp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp;

	      // Diffusion
	      apn[n] += -viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;

	      ani[n] = neighbor;
	      n++;

	      /*
	         if (parameter.orthof != 0.0)
	         {
	         gradun = Gradient (&xu, &xuf, LOGICAL_TRUE, neighbor);
	         gradvn = Gradient (&xv, &xvf, LOGICAL_TRUE, neighbor);
	         gradwn = Gradient (&xw, &xwf, LOGICAL_TRUE, neighbor);

	         // Non-orthogonal correction terms           
	         bpu += parameter.orthof * 
	         -viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp *
	         (GeoDotVectorVector (gradun, GeoSubVectorVector (faces[face].rnl, elements[neighbor].celement)) -
	         GeoDotVectorVector (gradup, GeoSubVectorVector (faces[face].rpl, elements[element].celement)));

	         bpv += parameter.orthof * 
                -viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp *
	         (GeoDotVectorVector (gradvn, GeoSubVectorVector (faces[face].rnl, elements[neighbor].celement)) -
	         GeoDotVectorVector (gradvp, GeoSubVectorVector (faces[face].rpl, elements[element].celement)));

	         bpw += parameter.orthof * 
	         -viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp *
	         (GeoDotVectorVector (gradwn, GeoSubVectorVector (faces[face].rnl, elements[neighbor].celement)) -
		 GeoDotVectorVector (gradwp, GeoSubVectorVector (faces[face].rpl, elements[element].celement)));
	         }
	       */

	    }
	  else
	    {

	      if (faces[face].bc != EMPTY)
		{

		  viscj = V_GetCmp (&visc, element + 1);

		  // Diffusion
		  if ((faces[face].bc != PERMEABLE || V_GetCmp (&xs, element + 1) > 0.5) && faces[face].bc != SLIP)
		    {	
		      app += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;

		      bpu += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) * V_GetCmp (&xuf, face + 1) / elements[element].Vp;
		      bpv += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) * V_GetCmp (&xvf, face + 1) / elements[element].Vp;
		      bpw += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) * V_GetCmp (&xwf, face + 1) / elements[element].Vp;
		    }

		  // Convection
		  if (parameter.scheme[iu] == UDS)
		    {

		      // UDS
		      if (V_GetCmp (&uf, face + 1) > 0.0)
			{

			  app += densp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp;

			}
		      else
			{

			  // Convection               
			  bpu += -densp * V_GetCmp (&uf, face + 1) * faces[face].Aj * V_GetCmp (&xuf, face + 1) / elements[element].Vp;
			  bpv += -densp * V_GetCmp (&uf, face + 1) * faces[face].Aj * V_GetCmp (&xvf, face + 1) / elements[element].Vp;
			  bpw += -densp * V_GetCmp (&uf, face + 1) * faces[face].Aj * V_GetCmp (&xwf, face + 1) / elements[element].Vp;

			}

		    }
		  else
		    {

		      // CDS
		      bpu += -densp * V_GetCmp (&uf, face + 1) * faces[face].Aj * V_GetCmp (&xuf, face + 1) / elements[element].Vp;
		      bpv += -densp * V_GetCmp (&uf, face + 1) * faces[face].Aj * V_GetCmp (&xvf, face + 1) / elements[element].Vp;
		      bpw += -densp * V_GetCmp (&uf, face + 1) * faces[face].Aj * V_GetCmp (&xwf, face + 1) / elements[element].Vp;

		    }

		  // Non-orthogonal correction terms          
		  /*
		     bpu +=
		     viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp * 
		     GeoDotVectorVector (gradup, GeoSubVectorVector(faces[face].rpl, elements[element].celement));

		     bpv +=
		     viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp * 
		     GeoDotVectorVector (gradvp, GeoSubVectorVector(faces[face].rpl, elements[element].celement));

		     bpw +=
		     viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp * 
		     GeoDotVectorVector (gradwp, GeoSubVectorVector(faces[face].rpl, elements[element].celement));
		   */

		}

	    }

	}

      if (dt > 0)
	{

	  // Unsteady term - Euler

	  app += densp / dt * parameter.st;

	  bpu += densp / dt * V_GetCmp (&xu0, element + 1);
	  bpv += densp / dt * V_GetCmp (&xv0, element + 1);
	  bpw += densp / dt * V_GetCmp (&xw0, element + 1);

	}

      // Source - viscous term
      /*
         gradvisc = Gradient (&visc, NULL, LOGICAL_FALSE, element);

         bpu += GeoDotVectorVector (gradup, gradvisc);
         bpv += GeoDotVectorVector (gradvp, gradvisc);
         bpw += GeoDotVectorVector (gradwp, gradvisc);
       */

      // Source - gravity
      bpu += densp * g.x;
      bpv += densp * g.y;
      bpw += densp * g.z;

      // Initialize H with source contribution without pressure 
      V_SetCmp (&hu, element + 1, bpu);
      V_SetCmp (&hv, element + 1, bpv);
      V_SetCmp (&hw, element + 1, bpw);

      // Source - pressure 
      /* 	
      gradp = Gradient (&xp, &xpf, LOGICAL_TRUE, element);
      bpu += -gradp.x;
      bpv += -gradp.y;
      bpw += -gradp.z;
      */

      if (app == 0.0 || app != app)
	{
	  printf ("\nError: Problem setting up momentum matrix\n");
	  exit (LOGICAL_ERROR);
	}

      V_SetCmp (&ap, element + 1, app);

      Q_SetLen (&Am, element + 1, n + 1);

      Q_SetEntry (&Am, element + 1, 0, element + 1, app);

      for (j = 0; j < n; j++)
	{
	  Q_SetEntry (&Am, element + 1, j + 1, ani[j] + 1, apn[j]);
	}

      if (parameter.calc[iu] == LOGICAL_TRUE)
	{
	  V_SetCmp (&bu, element + 1, bpu);
	}

      if (parameter.calc[iv] == LOGICAL_TRUE)
	{
	  V_SetCmp (&bv, element + 1, bpv);
	}

      if (parameter.calc[iw] == LOGICAL_TRUE)
	{
	  V_SetCmp (&bw, element + 1, bpw);
	}

    }

}

void
CorrectVelocity (char *var, int *fiter, double dt, double maxCp, int verbose,
		 int pchecks)
{

  if (parameter.calc[ip] == LOGICAL_TRUE) 
    {

      // Correct face values
      CorrectFaceUVW ();

      // Correct cell center
      CorrectVelocityField ();

   }
   
}

void
CalculateVelocity (char *var, int *fiter, double dt, double maxCp,
		   int verbose, int pchecks)
{

  double mres;
  int miter;
  double mtime;

  Q_Constr (&Am, "Momentum matrix", nbelements, False, Rowws, Normal, True);
  V_Constr (&bu, "Momentum source x-component", nbelements, Normal, True);
  V_Constr (&bv, "Momentum source y-component", nbelements, Normal, True);
  V_Constr (&bw, "Momentum source z-component", nbelements, Normal, True);

  // Store previous time step values
  Asgn_VV (&xu0, &xu);
  Asgn_VV (&xv0, &xv);
  Asgn_VV (&xw0, &xw);

  // Build three momentum matrices for u, v, w velocity components
  if (parameter.calc[ip] == LOGICAL_TRUE)
    {
      BuildMomentumMatrix (dt);
    }

  if (pchecks == LOGICAL_TRUE)
    {
      if (!CheckIfDiagonalMatrix (&Am))
	{
	  printf ("\nWarning: Momentum matrix is not diagonal dominant\n");
	  WriteMatrix (&Am, LOGICAL_FALSE);
	  WriteVector (&bu);
	  WriteVector (&bv);
	  WriteVector (&bw);
	  //exit (LOGICAL_ERROR);
	}
    }

  if (parameter.calc[iu] == LOGICAL_TRUE)
    {

      fiter[iu]++;

      // Set matrix solution accuracy
      SetRTCAccuracy (parameter.mtol[iu]);

      // Solve matrix for u velocity component
      SolveMatrix (&Am, &xu, &bu, &miter, &mres, &mtime, parameter.msolver[iu], parameter.mprecond[iu], parameter.miter[iu]);

      if (verbose == LOGICAL_TRUE)
	printf
	  ("\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
	   var[iu], miter, mres, mtime);

      if ((mres > parameter.mtol[iu] && miter == parameter.miter[iu])
	  || LASResult () != LASOK)
	{
	  printf ("\nError: Problem solving matrix %c\n", var[iu]);
	  exit (LOGICAL_ERROR);
	}

    }

  if (parameter.calc[iv] == LOGICAL_TRUE)
    {

      fiter[iv]++;

      // Set matrix solution accuracy
      SetRTCAccuracy (parameter.mtol[iv]);

      // Solve matrix for v velocity component
      SolveMatrix (&Am, &xv, &bv, &miter, &mres, &mtime, parameter.msolver[iv], parameter.mprecond[iv], parameter.miter[iv]);

      if (verbose == LOGICAL_TRUE)
	printf
	  ("\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
	   var[iv], miter, mres, mtime);

      if (mres > parameter.mtol[iv] && miter == parameter.miter[iv])
	{
	  printf ("\nError: Problem solving matrix %c\n", var[iv]);
	  exit (LOGICAL_ERROR);
	}

    }

  if (parameter.calc[iw] == LOGICAL_TRUE)
    {

      fiter[iw]++;

      // Set matrix solution accuracy
      SetRTCAccuracy (parameter.mtol[iw]);

      // Solve matrix for w velocity component
      SolveMatrix (&Am, &xw, &bw, &miter, &mres, &mtime, parameter.msolver[iw], parameter.mprecond[iw], parameter.miter[iw]);

      if (verbose == LOGICAL_TRUE)
	printf
	  ("\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
	   var[iw], miter, mres, mtime);

      if ((mres > parameter.mtol[iw] && miter == parameter.miter[iw])
	  || LASResult () != LASOK)
	{
	  printf ("\nProblem solving matrix %c\n", var[iw]);
	  exit (LOGICAL_ERROR);
	}

    }

  // Calculate correction factors
  if (parameter.calc[ip] == LOGICAL_TRUE)
    {
      CalculateCorrectionFactors (&Am, &xu, &xv, &xw, &hu, &hv, &hw);
    }

  Q_Destr (&Am);
  V_Destr (&bu);
  V_Destr (&bv);
  V_Destr (&bw);

}
