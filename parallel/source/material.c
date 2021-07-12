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
#include <string.h>
#include <math.h>

#include "petsc.h"

#include "ioutils.h"
#include "globals.h"
#include "param.h"
#include "material.h"

int
MtlImportMTL (char *file)
{

  int i, n;

  int inull;

  int tcode;
  int nbmat;

  FILE *fp;
  char descr[512];

  fp = fopen (file, "r");

  if (fp == NULL)
    {
      PetscPrintf (PETSC_COMM_WORLD, "\nError: Material file not found!\n");
      PetscPrintf (PETSC_COMM_WORLD, "%s\n\n", file);
      exit (LOGICAL_ERROR);
    }

  PetscPrintf (PETSC_COMM_WORLD, "\nReading material file: %s ...\n", file);

  do
    {

      do
	{
	  fscanf (fp, "%s", descr);

	  if (strcmp (descr, "$Material") == 0)
	    break;

	  if (strcmp (descr, "$EndFile") == 0)
	    break;

	}
      while (!feof (fp));

      if (strcmp (descr, "$EndFile") == 0)
	break;

      if (strcmp (descr, "$Material") == 0)
	{

	  fscanf (fp, "%d %d", &inull, &nbmat);

	  GetLine (fp);

	  for (i = 0; i < nbmat; i++)
	    {

	      fscanf (fp, "%d %d", &tcode, &n);

	      switch (tcode)
		{
		case 20012:

		  // Compressibility of fluid 0
		  GetLine (fp);
		  fscanf (fp, "%f", &material.psi[0].cpsi);

		  break;

		case 20015:

		  // Density of fluid 0
		  GetLine (fp);
		  fscanf (fp, "%f", &material.dens[0].cdens);

		  break;

		case 20021:

		  // Viscosity of fluid 0
		  GetLine (fp);
		  fscanf (fp, "%f", &material.visc[0].cvisc);

		  break;

		case 20022:

		  // Specific heat of fluid 0
		  GetLine (fp);
		  fscanf (fp, "%f", &material.therm[0].cspheat);

		  break;

		case 20023:

		  // Thermal conductivity of fluid 0
		  GetLine (fp);
		  fscanf (fp, "%f", &material.therm[0].cthcond);

		  break;

		case 21012:

		  // Compressibility of fluid 1
		  GetLine (fp);
		  fscanf (fp, "%f", &material.psi[1].cpsi);

		  break;

		case 21015:

		  // Density of fluid 1
		  GetLine (fp);
		  fscanf (fp, "%f", &material.dens[1].cdens);

		  break;

		case 21021:

		  // Viscosity of fluid 1
		  GetLine (fp);
		  fscanf (fp, "%f", &material.visc[1].cvisc);

		  break;

		case 21022:

		  // Specific heat of fluid 1
		  GetLine (fp);
		  fscanf (fp, "%f", &material.therm[1].cspheat);

		  break;

		case 21023:

		  // Thermal conductivity of fluid 1
		  GetLine (fp);
		  fscanf (fp, "%f", &material.therm[1].cthcond);

		  break;

		case 26050:

		  // Constant surface tension (fluid 0 / fluid 1)
		  GetLine (fp);
		  fscanf (fp, "%f", &material.tens);

		  break;

		case 26060:

		  // Thermal conductivity (boundary)
		  GetLine (fp);
		  fscanf (fp, "%f", &material.bthcond);

		  break;

		default:

		  PetscPrintf (PETSC_COMM_WORLD,
			       "\nError: Unknown material code.\n");
		  exit (LOGICAL_ERROR);
		  break;

		}

	    }

	}

    }
  while (!feof (fp));

  fclose (fp);

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "- Fluid 0\n");
  PetscPrintf (PETSC_COMM_WORLD, "Density: \t\t\t\t\t%+.3E %s/%s^3\n",
	       material.dens[0].cdens, parameter.umass, parameter.ulength);
  PetscPrintf (PETSC_COMM_WORLD, "Viscosity: \t\t\t\t\t%+.3E %s/(%s.%s)\n",
	       material.visc[0].cvisc, parameter.umass, parameter.ulength,
	       parameter.utime);
  PetscPrintf (PETSC_COMM_WORLD,
	       "Specific heat: \t\t\t\t\t%+.3E %s/(%s.%s)\n",
	       material.therm[0].cspheat, parameter.uenergy, parameter.umass,
	       parameter.utemperature);
  PetscPrintf (PETSC_COMM_WORLD,
	       "Thermal conductivity: \t\t\t\t%+.3E %s/(%s.%s.%s)\n",
	       material.therm[0].cthcond, parameter.uenergy, parameter.utime,
	       parameter.umass, parameter.utemperature);

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "- Fluid 1\n");
  PetscPrintf (PETSC_COMM_WORLD, "Density: \t\t\t\t\t%+.3E %s/%s^3\n",
	       material.dens[1].cdens, parameter.umass, parameter.ulength);
  PetscPrintf (PETSC_COMM_WORLD, "Viscosity: \t\t\t\t\t%+.3E %s/(%s.%s)\n",
	       material.visc[1].cvisc, parameter.umass, parameter.ulength,
	       parameter.utime);
  PetscPrintf (PETSC_COMM_WORLD,
	       "Specific heat: \t\t\t\t\t%+.3E %s/(%s.%s)\n",
	       material.therm[1].cspheat, parameter.uenergy, parameter.umass,
	       parameter.utemperature);
  PetscPrintf (PETSC_COMM_WORLD,
	       "Thermal conductivity: \t\t\t\t%+.3E %s/(%s.%s.%s)\n",
	       material.therm[1].cthcond, parameter.uenergy, parameter.utime,
	       parameter.umass, parameter.utemperature);

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "- Fluid 0 / Fluid 1\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "Surace tension: \t\t\t\t%+.3E %s/(%s.%s^2)\n", material.tens,
	       parameter.umass, parameter.ulength, parameter.utime);

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "- Boundary\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "Thermal conductivity: \t\t\t\t%+.3E %s/(%s.%s.%s)\n",
	       material.bthcond, parameter.uenergy, parameter.utime,
	       parameter.umass, parameter.utemperature);

  PetscPrintf (PETSC_COMM_WORLD, "Done.\n");

  return LOGICAL_TRUE;

}
