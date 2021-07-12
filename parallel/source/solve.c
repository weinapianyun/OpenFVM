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

#include <stdio.h>
#include <stdlib.h>

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

#include "gamma.h"
#include "velocity.h"
#include "pressure.h"
#include "temperature.h"

#include "solve.h"

void
AllocateMemory ()
{

  V_Constr (&cex, nbelements, 0);	// Cell centers x-component
  V_Constr (&cey, nbelements, 0);	// Cell centers y-component
  V_Constr (&cez, nbelements, 0);	// Cell centers z-component

  V_Constr (&Co, nbelements, 0);	// Courant number
  V_Constr (&uf, nbfaces, 1);	// Face flux velocity

  V_Constr (&dens, nbelements, 0);	// Density      
  V_Constr (&visc, nbelements, 0);	// Dynamic viscosity
  V_Constr (&thcond, nbelements, 0);	// Thermal conductivity
  V_Constr (&spheat, nbelements, 0);	// Specific heat

  V_Constr (&xu0, nbelements, 0);	// Velocity x-component at cell center (previous time step)
  V_Constr (&xv0, nbelements, 0);	// Velocity x-component at cell center (previous time step)
  V_Constr (&xw0, nbelements, 0);	// Velocity z-component at cell center (previous time step)
  V_Constr (&xp0, nbelements, 0);	// Pressure at cell center (previous time step)
  V_Constr (&xT0, nbelements, 0);	// Temperature at cell center (previous time step)
  V_Constr (&xs0, nbelements, 0);	// Gamma at cell center (previous time step)

  V_Constr (&xu, nbelements, 0);	// Velocity x-component at cell center
  V_Constr (&xv, nbelements, 0);	// Velocity y-component at cell center (previous time step)
  V_Constr (&xw, nbelements, 0);	// Velocity z-component at cell center
  V_Constr (&xp, nbelements, 0);	// Pressure at cell center
  V_Constr (&xT, nbelements, 0);	// Temperature at cell center   
  V_Constr (&xs, nbelements, 0);	// Gamma at cell center

  V_Constr (&xuf, nbfaces, 1);	// Velocity x-component at face center
  V_Constr (&xvf, nbfaces, 1);	// Velocity x-component at face center
  V_Constr (&xwf, nbfaces, 1);	// Velocity z-component at face center
  V_Constr (&xpf, nbfaces, 1);	// Pressure at face center
  V_Constr (&xTf, nbfaces, 1);	// Temperature at face center
  V_Constr (&xsf, nbfaces, 1);	// Gamma at face center

  V_Constr (&ap, nbelements, 0);	// Momentum matrix diagonal
  V_Constr (&hu, nbelements, 0);	// Momentum matrix source x-component without pressure
  V_Constr (&hv, nbelements, 0);	// Momentum matrix source y-component without pressure
  V_Constr (&hw, nbelements, 0);	// Momentum matrix source z-component without pressure

  V_Constr (&temp1, nbelements, 0);	// Temporary vector 1
  V_Constr (&temp2, nbelements, 0);	// Temporary vector 2

}

void
DeallocateMemory ()
{

  V_Destr (&cex);
  V_Destr (&cey);
  V_Destr (&cez);

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

  V_Destr (&temp1);
  V_Destr (&temp2);

}

void
CheckMassConservationError (double dt) //检查质量守恒误差
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

	  mcp += V_GetCmp (&uf, face) * faces[face].Aj;

	}

      sum += LABS (mcp);

    }

  PetscPrintf (PETSC_COMM_WORLD, "\nMass conservation error: %+E kg\n", sum);


}

void Solve (char *var, int *fiter, double dt, double *maxCp, int verbose, int pchecks)
{


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
