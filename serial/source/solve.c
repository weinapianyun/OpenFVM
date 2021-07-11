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

#include "gamma.h"
#include "velocity.h"
#include "pressure.h"
#include "temperature.h"

#include "solve.h"

void
AllocateMemory ()
{

  V_Constr (&Co, "Courant number", nbelements, Normal, True);
  V_Constr (&uf, "Face flux velocity", nbfaces, Normal, True);

  V_Constr (&dens, "Density", nbelements, Normal, True);
  V_Constr (&visc, "Dynamic viscosity", nbelements, Normal, True);
  V_Constr (&thcond, "Thermal conductivity", nbelements, Normal, True);
  V_Constr (&spheat, "Specific heat", nbelements, Normal, True);

  V_Constr (&xu0, "Velocity x-component at cell center (previous time step)", nbelements, Normal, True);
  V_Constr (&xv0, "Velocity y-component at cell center (previous time step)", nbelements, Normal, True);
  V_Constr (&xw0, "Velocity z-component at cell center (previous time step)", nbelements, Normal, True);
  V_Constr (&xp0, "Pressure at cell center (previous time step)", nbelements, Normal, True);
  V_Constr (&xT0, "Temperature at cell center (previous time step)", nbelements, Normal, True);
  V_Constr (&xs0, "Gamma at cell center (previous time step)", nbelements, Normal, True);

  V_Constr (&xu, "Velocity x-component at cell center", nbelements, Normal, True);
  V_Constr (&xv, "Velocity y-component at cell center", nbelements, Normal, True);
  V_Constr (&xw, "Velocity z-component at cell center", nbelements, Normal, True);
  V_Constr (&xp, "Pressure at cell center", nbelements, Normal, True);
  V_Constr (&xT, "Temperature at cell center", nbelements, Normal, True);
  V_Constr (&xs, "Gamma at cell center", nbelements, Normal, True);

  V_Constr (&xuf, "Velocity x-component at face center", nbfaces, Normal, True);
  V_Constr (&xvf, "Velocity y-component at face center", nbfaces, Normal, True);
  V_Constr (&xwf, "Velocity z-component at face center", nbfaces, Normal, True);
  V_Constr (&xpf, "Pressure at face center", nbfaces, Normal, True);
  V_Constr (&xTf, "Temperature at face center", nbfaces, Normal, True);
  V_Constr (&xsf, "Gamma at face center", nbfaces, Normal, True);

  V_Constr (&ap, "Momentum matrix diagonal", nbelements, Normal, True);
  V_Constr (&hu, "Momentum matrix source x-component without pressure", nbelements, Normal, True);
  V_Constr (&hv, "Momentum matrix source y-component without pressure", nbelements, Normal, True);
  V_Constr (&hw, "Momentum matrix source z-component without pressure", nbelements, Normal, True);

}

void
DeallocateMemory ()
{

  V_Destr (&Co);
  V_Destr (&uf);

  V_Destr (&dens);
  V_Destr (&visc);
  V_Destr (&thcond);
  V_Destr (&spheat);

  V_Destr (&xu0);
  V_Destr (&xv0);
  V_Destr (&xw0);
  V_Destr (&xp0);
  V_Destr (&xT0);
  V_Destr (&xs0);

  V_Destr (&xu);
  V_Destr (&xv);
  V_Destr (&xw);
  V_Destr (&xp);
  V_Destr (&xT);
  V_Destr (&xs);

  V_Destr (&xuf);
  V_Destr (&xvf);
  V_Destr (&xwf);
  V_Destr (&xpf);
  V_Destr (&xTf);
  V_Destr (&xsf);

  V_Destr (&ap);
  V_Destr (&hu);
  V_Destr (&hv);
  V_Destr (&hw);

}

void
CheckMassConservationError (double dt)
{

  int i, j;

  int face, pair;

  int element;

  double mcp;

  double sum;

  sum = 0.0;

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      mcp = 0.0;

      for (j = 0; j < elements[element].nbfaces; j++)
	{

	  face = elements[element].face[j];

	  pair = faces[face].pair;

	  mcp += V_GetCmp (&uf, face + 1) * faces[face].Aj;

	}

      sum += LABS (mcp);

    }

  printf ("\nMass conservation error: %+E kg\n", sum);

}

void Solve (char *var, int *fiter, double dt, double *maxCp, int verbose, int pchecks)
{

  int i;

  CalculateGamma (var, fiter, dt, maxCp, verbose, pchecks);

  // Set material properties 
  SetMaterialProperties ();

  V_SetAllCmp (&ap, 1.0);

  V_SetAllCmp (&hu, 0.0);
  V_SetAllCmp (&hv, 0.0);
  V_SetAllCmp (&hw, 0.0);

  CalculateVelocity (var, fiter, dt, *maxCp, verbose, pchecks);

  CalculatePressure (var, fiter, dt, *maxCp, verbose, pchecks);

  CorrectVelocity (var, fiter, dt, *maxCp, verbose, pchecks);

  if (pchecks == LOGICAL_TRUE)
    {
      // Check mass conservation
      CheckMassConservationError (dt);
    }

  CalculateTemperature (var, fiter, dt, *maxCp, verbose, pchecks);

}
