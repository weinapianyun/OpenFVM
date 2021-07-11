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

#include "petsc.h"

#include "ioutils.h"
#include "globals.h"
#include "param.h"

void
SetDefaults ()
{

  strcpy (parameter.ulength, "m");
  strcpy (parameter.umass, "kg");
  strcpy (parameter.utime, "s");
  strcpy (parameter.uenergy, "J");
  strcpy (parameter.utemperature, "K");

  parameter.inertia = 1;

  parameter.ef[0] = 1.0;
  parameter.ef[1] = 1.0;
  parameter.ef[2] = 1.0;
  parameter.ef[3] = 1.0;
  parameter.ef[4] = 1.0;
  parameter.ef[5] = 1.0;

  parameter.dfactor = 1.0;

  parameter.st = 1.0;

  parameter.timemethod[0] = 1;
  parameter.timemethod[1] = 1;
  parameter.timemethod[2] = 1;
  parameter.timemethod[3] = 1;
  parameter.timemethod[4] = 1;
  parameter.timemethod[5] = 1;

  parameter.scheme[0] = 1;
  parameter.scheme[1] = 1;
  parameter.scheme[2] = 1;
  parameter.scheme[3] = 1;
  parameter.scheme[4] = 1;
  parameter.scheme[5] = 1;

  parameter.steady = 0;

  parameter.ftol[0] = 1E-6;
  parameter.ftol[1] = 1E-6;
  parameter.ftol[2] = 1E-6;
  parameter.ftol[3] = 1E-6;
  parameter.ftol[4] = 1E-6;
  parameter.ftol[5] = 1E-6;

  parameter.wbinary = 0;
  parameter.nsav = 1;

  parameter.calc[0] = 0;
  parameter.calc[1] = 0;
  parameter.calc[2] = 0;
  parameter.calc[3] = 0;
  parameter.calc[4] = 0;
  parameter.calc[5] = 0;

  parameter.savflux = 0;

  parameter.fsav[0] = 0;
  parameter.fsav[1] = 0;
  parameter.fsav[2] = 0;
  parameter.fsav[3] = 0;
  parameter.fsav[4] = 0;
  parameter.fsav[5] = 0;

  parameter.csav[0] = 0;
  parameter.csav[1] = 0;
  parameter.csav[2] = 0;
  parameter.csav[3] = 0;
  parameter.csav[4] = 0;
  parameter.csav[5] = 0;

  parameter.probe[0] = 0;
  parameter.probe[1] = 0;
  parameter.probe[2] = 0;
  parameter.probe[3] = 0;
  parameter.probe[4] = 0;
  parameter.probe[5] = 0;

  parameter.smooth = 1;

  parameter.vortex[0] = 0;
  parameter.vortex[1] = 0;
  parameter.vortex[2] = 0;

  parameter.streamf = 0;

  parameter.fvec = 0;
  parameter.cvec = 0;

  parameter.kq = 2.0;
  parameter.ncicsamcor = 2;

  parameter.g[0] = 0.0;
  parameter.g[1] = 0.0;
  parameter.g[2] = 0.0;

  parameter.msolver[0] = 8;
  parameter.msolver[1] = 8;
  parameter.msolver[2] = 8;
  parameter.msolver[3] = 8;
  parameter.msolver[4] = 8;
  parameter.msolver[5] = 3;

  parameter.mprecond[0] = 4;
  parameter.mprecond[1] = 4;
  parameter.mprecond[2] = 4;
  parameter.mprecond[3] = 4;
  parameter.mprecond[4] = 4;
  parameter.mprecond[5] = 4;

  parameter.northocor = 10;
  parameter.orthof = 1.0;

  parameter.mtol[0] = 1E-8;
  parameter.mtol[1] = 1E-8;
  parameter.mtol[2] = 1E-8;
  parameter.mtol[3] = 1E-8;
  parameter.mtol[4] = 1E-8;
  parameter.mtol[5] = 1E-8;

  parameter.miter[0] = 500;
  parameter.miter[1] = 500;
  parameter.miter[2] = 500;
  parameter.miter[3] = 500;
  parameter.miter[4] = 500;
  parameter.miter[5] = 500;

  parameter.restart = 10000;
  parameter.adjdt = 0;
  parameter.maxCp = 0.25;

  parameter.t0 = 0.0;
  parameter.t1 = 0.001;
  parameter.dt = 0.001;

  parameter.fill = 0;
  parameter.pf = 99.5;

}

int
ParImportPAR (char *file)
{

  int i, j, n;

  int inull;

  int tcode;
  int nbpar;

  FILE *fp;

  char descr[512];

  //Set default values
  SetDefaults ();

  fp = fopen (file, "r");

  if (fp == NULL)
    {
      PetscPrintf (PETSC_COMM_WORLD, "\nError: Parameter file not found!\n");
      PetscPrintf (PETSC_COMM_WORLD, "%s\n\n", file);
      exit (LOGICAL_ERROR);
    }

  PetscPrintf (PETSC_COMM_WORLD, "\nReading parameters file: %s ...\n", file);

  PetscPrintf (PETSC_COMM_WORLD, "\n");

  do
    {

      do
	{
	  fscanf (fp, "%s", descr);

	  if (strcmp (descr, "$Parameter") == 0)
	    break;

	  if (strcmp (descr, "$EndFile") == 0)
	    break;

	}
      while (!feof (fp));

      if (strcmp (descr, "$EndFile") == 0)
	break;

      if (strcmp (descr, "$Parameter") == 0)
	{

	  fscanf (fp, "%d %d ", &inull, &nbpar);

	  GetLine (fp);

	  for (i = 0; i < nbpar; i++)
	    {

	      fscanf (fp, "%d %d", &tcode, &n);

	      switch (tcode)
		{
		case 30001:

		  // Units (length, mass, time)
		  GetLine (fp);
		  fscanf (fp, "%s", parameter.ulength);
		  fscanf (fp, "%s", parameter.umass);
		  fscanf (fp, "%s", parameter.utime);
		  fscanf (fp, "%s", parameter.uenergy);
		  fscanf (fp, "%s", parameter.utemperature);

		  PetscPrintf (PETSC_COMM_WORLD, "Units:\n");
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Length:      \t\t\t\t\t[%s]\n",
			       parameter.ulength);
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Mass:        \t\t\t\t\t[%s]\n",
			       parameter.umass);
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Time:        \t\t\t\t\t[%s]\n",
			       parameter.utime);
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Time:        \t\t\t\t\t[%s]\n",
			       parameter.uenergy);
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Temperature: \t\t\t\t\t[%s]\n",
			       parameter.utemperature);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30003:

		  // Inertia
		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.inertia);

		  if (parameter.inertia == 1)
		    PetscPrintf (PETSC_COMM_WORLD, "Solve momentum\n");
		  else
		    PetscPrintf (PETSC_COMM_WORLD, "Do not solve momentum\n");

		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30005:

		  // Convection interpolation scheme
		  GetLine (fp);
		  for (j = 0; j < nphi; j++)
		    fscanf (fp, "%d", &parameter.scheme[j]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Variable: \t\t\t\t\t[ u v w p T s]\n");
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Interpolation scheme a): \t\t\t[");
		  for (j = 0; j < nphi; j++)
		    PetscPrintf (PETSC_COMM_WORLD, " %d",
				 parameter.scheme[j]);
		  PetscPrintf (PETSC_COMM_WORLD, "]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30009:

		  // Diffusion factor
		  GetLine (fp);
		  fscanf (fp, "%f", &parameter.dfactor);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Diffusion factor: \t\t\t\t%.3f\n",
			       parameter.dfactor);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30010:

		  // Restart files
		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.restart);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Save restart file every %d iterations\n",
			       parameter.restart);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30020:

		  // Binary output
		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.wbinary);

		  if (parameter.wbinary == LOGICAL_TRUE)
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Save results in binary format.\n");
		  else
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Save results in ascii format.\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30030:

		  // Explicit method 
		  GetLine (fp);
		  for (j = 0; j < nphi; j++)
		    fscanf (fp, "%d", &parameter.timemethod[j]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Variable: \t\t\t\t\t[ u v w p T s]\n");
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Time advancement method: b) \t\t\t[");
		  for (j = 0; j < nphi; j++)
		    PetscPrintf (PETSC_COMM_WORLD, " %d",
				 parameter.timemethod[j]);
		  PetscPrintf (PETSC_COMM_WORLD, "]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30031:

		  // Explicit factor
		  GetLine (fp);
		  for (j = 0; j < nphi; j++)
		    fscanf (fp, "%f", &parameter.ef[j]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Variable: \t\t\t\t\t[ u v w p T s]\n");
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Explicit factor: \t\t\t\t[");
		  for (j = 0; j < nphi; j++)
		    PetscPrintf (PETSC_COMM_WORLD, " %.1f", parameter.ef[j]);
		  PetscPrintf (PETSC_COMM_WORLD, "]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30040:

		  // Calculate nth variable 
		  GetLine (fp);
		  for (j = 0; j < nphi; j++)
		    fscanf (fp, "%d", &parameter.calc[j]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Variable: \t\t\t\t\t[ u v w p T s]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "Calculate: \t\t\t\t\t[");
		  for (j = 0; j < nphi; j++)
		    PetscPrintf (PETSC_COMM_WORLD, " %d", parameter.calc[j]);
		  PetscPrintf (PETSC_COMM_WORLD, "]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30090:

		  // Stability factor 
		  GetLine (fp);
		  fscanf (fp, "%f", &parameter.st);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Stability factor: \t\t\t\t%.3f\n",
			       parameter.st);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30100:

		  // Steady state
		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.steady);

		  if (parameter.steady == LOGICAL_TRUE)
		    {
		      PetscPrintf (PETSC_COMM_WORLD,
				   "Calculate until steady state is reached.\n");
		    }
		  else
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Stop simulation at end time.\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30105:

		  // Convergence criterion for steady state    
		  GetLine (fp);
		  for (j = 0; j < nphi; j++)
		    fscanf (fp, "%f", &parameter.ftol[j]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Variable: \t\t\t\t\t[ u v w p T s]\n");
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Convergence criterion for steady state: \t[");
		  for (j = 0; j < nphi; j++)
		    PetscPrintf (PETSC_COMM_WORLD, " %+.3E",
				 parameter.ftol[j]);
		  PetscPrintf (PETSC_COMM_WORLD, "]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30200:

		  // Adjust time interval

		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.adjdt);

		  if (parameter.adjdt == LOGICAL_TRUE)
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Adjust time intervals.\n");
		  else
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Do not adjust time intervals.\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30201:

		  // Maximum Courant number
		  GetLine (fp);
		  fscanf (fp, "%f", &parameter.maxCp);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Maximum Courant number: \t\t\t%.3f\n",
			       parameter.maxCp);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30400:

		  // Number of saves
		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.nsav);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Number of saves: \t\t\t\t%d\n",
			       parameter.nsav);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30440:

		  // Save face flux
		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.savflux);

		  PetscPrintf (PETSC_COMM_WORLD, "Save flux: \t\t\t\t\t%d\n",
			       parameter.savflux);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30450:

		  // Write face scalars (u v w p T s)
		  GetLine (fp);
		  for (j = 0; j < nphi; j++)
		    fscanf (fp, "%d", &parameter.fsav[j]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Variable: \t\t\t\t\t[ u v w p T s]\n");
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Save scalars on face: \t\t\t\t[");
		  for (j = 0; j < nphi; j++)
		    PetscPrintf (PETSC_COMM_WORLD, " %d", parameter.fsav[j]);
		  PetscPrintf (PETSC_COMM_WORLD, "]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30455:

		  //Write face vectors (uvw)
		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.fvec);

		  if (parameter.fvec == LOGICAL_TRUE)
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Save face - vector magnitude: \t\t\t[yes]\n");
		  else
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Save face - vector magnitude: \t\t\t[no]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30460:

		  // Write element scalars (u v w p T s)
		  GetLine (fp);
		  for (j = 0; j < nphi; j++)
		    fscanf (fp, "%d", &parameter.csav[j]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Variable: \t\t\t\t\t[ u v w p T s]\n");
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Save scalars in cell: \t\t\t\t[");
		  for (j = 0; j < nphi; j++)
		    PetscPrintf (PETSC_COMM_WORLD, " %d", parameter.csav[j]);
		  PetscPrintf (PETSC_COMM_WORLD, "]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30465:

		  // Write element vectors (uvw)
		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.cvec);

		  if (parameter.cvec == LOGICAL_TRUE)
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Save cell center - vector: \t\t\t[yes]\n");
		  else
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Save cell center - vector: \t\t\t[no]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30470:

		  // Write vorticity (x y z)
		  GetLine (fp);
		  for (j = 0; j < 3; j++)
		    fscanf (fp, "%d", &parameter.vortex[j]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Axis: \t\t\t\t\t\t[ x y z]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "Save vorticity: \t\t\t\t[");
		  for (j = 0; j < 3; j++)
		    PetscPrintf (PETSC_COMM_WORLD, " %d",
				 parameter.vortex[j]);
		  PetscPrintf (PETSC_COMM_WORLD, "]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30475:

		  // Write stream function (xy)
		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.streamf);

		  if (parameter.streamf == LOGICAL_TRUE)
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Save stream function in xy plane\n");
		  else
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Do not save stream function in xy plane\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30485:

		  // Probe (u v w p T s)
		  GetLine (fp);
		  for (j = 0; j < nphi; j++)
		    fscanf (fp, "%d", &parameter.probe[j]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Variable: \t\t\t\t\t[ u v w p T s]\n");
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Probe options: \t\t\t\t\t[");
		  for (j = 0; j < nphi; j++)
		    PetscPrintf (PETSC_COMM_WORLD, " %d", parameter.probe[j]);
		  PetscPrintf (PETSC_COMM_WORLD, "]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30500:

		  // Smooth
		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.smooth);

		  if (parameter.smooth == LOGICAL_TRUE)
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Smooth probe values: \t\t\t\t[yes]\n");
		  else
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Smooth probe values: \t\t\t\t[no]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30550:

		  // Maximum number of non-othorgonal corrections
		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.northocor);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Maximum number of non-orthogonal corrections: \t%d\n",
			       parameter.northocor);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30551:

		  // Orthogonal factor
		  GetLine (fp);
		  fscanf (fp, "%f", &parameter.orthof);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Orthogonal factor: \t\t\t\t%E\n",
			       parameter.orthof);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30600:

		  // Convergence criterion (matrix solution) 
		  GetLine (fp);
		  for (j = 0; j < nphi; j++)
		    fscanf (fp, "%f", &parameter.mtol[j]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Variable: \t\t\t\t\t[ u v w p T s]\n");
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Matrix solution tolerance: \t\t\t[");
		  for (j = 0; j < nphi; j++)
		    PetscPrintf (PETSC_COMM_WORLD, " %+.3E",
				 parameter.mtol[j]);
		  PetscPrintf (PETSC_COMM_WORLD, "]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30601:

		  // Maximum number of iterations (matrix solution)
		  GetLine (fp);
		  for (j = 0; j < nphi; j++)
		    fscanf (fp, "%d", &parameter.miter[j]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Variable: \t\t\t\t\t[ u v w p T s]\n");
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Matrix solution iterations: \t\t\t[");
		  for (j = 0; j < nphi; j++)
		    PetscPrintf (PETSC_COMM_WORLD, " %d", parameter.miter[j]);
		  PetscPrintf (PETSC_COMM_WORLD, "]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30650:

		  // Matrix solver (u v w p T s) (0-NONE, 1-JACOBI, 2-SOR, 3-QMR, 4-GMRES, 5-CG, 6-CGN, 7-CGS, 8-BICG, 9-BICGS)

		  GetLine (fp);
		  for (j = 0; j < nphi; j++)
		    fscanf (fp, "%d", &parameter.msolver[j]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Variable: \t\t\t\t\t[ u v w p T s]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "Solver: c) \t\t\t\t\t[");
		  for (j = 0; j < nphi; j++)
		    PetscPrintf (PETSC_COMM_WORLD, " %d",
				 parameter.msolver[j]);
		  PetscPrintf (PETSC_COMM_WORLD, "]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30651:

		  // Matrix preconditioner (0-NONE, 1-JACOBI, 2-SOR, 3-ILU)
		  GetLine (fp);
		  for (j = 0; j < nphi; j++)
		    fscanf (fp, "%d", &parameter.mprecond[j]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Variable: \t\t\t\t\t[ u v w p T s]\n");
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Pre-conditioner: d) \t\t\t\t[");
		  for (j = 0; j < nphi; j++)
		    PetscPrintf (PETSC_COMM_WORLD, " %d",
				 parameter.mprecond[j]);
		  PetscPrintf (PETSC_COMM_WORLD, "]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30800:

		  // Interface scheme factor - CICSAM
		  GetLine (fp);
		  fscanf (fp, "%f", &parameter.kq);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Interface scheme factor: \t\t\t%.3f\n",
			       parameter.kq);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 30900:

		  // Maximum number of CICSAM corrections

		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.ncicsamcor);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Number of CICSAM corrections: \t\t\t%d\n",
			       parameter.ncicsamcor);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 32000:

		  // Start time
		  GetLine (fp);
		  fscanf (fp, "%f", &parameter.t0);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Start time: \t\t\t\t\t%+.3E s\n",
			       parameter.t0);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 32001:

		  // End time
		  GetLine (fp);
		  fscanf (fp, "%f", &parameter.t1);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "End time: \t\t\t\t\t%+.3E s\n", parameter.t1);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 32002:

		  // Time interval
		  GetLine (fp);
		  fscanf (fp, "%f", &parameter.dt);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Time interval: \t\t\t\t\t%+.3E s\n",
			       parameter.dt);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 33050:

		  // Simulate filling
		  GetLine (fp);
		  fscanf (fp, "%d", &parameter.fill);

		  if (parameter.fill == LOGICAL_TRUE)
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Simulate filling: \t\t\t\t[yes]\n");
		  else
		    PetscPrintf (PETSC_COMM_WORLD,
				 "Simulate filling: \t\t\t\t[no]\n");
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 33051:

		  // Completely filled percentage
		  GetLine (fp);
		  fscanf (fp, "%f", &parameter.pf);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Stop filling at: \t\t\t\t%.2f%%\n",
			       parameter.pf);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		case 34000:

		  // Gravity vector
		  GetLine (fp);
		  fscanf (fp, "%f %f %f", &parameter.g[0], &parameter.g[1],
			  &parameter.g[2]);

		  PetscPrintf (PETSC_COMM_WORLD,
			       "Gravity vector (x): \t\t\t\t%+.3E m/s^2\n",
			       parameter.g[0]);
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Gravity vector (y): \t\t\t\t%+.3E m/s^2\n",
			       parameter.g[1]);
		  PetscPrintf (PETSC_COMM_WORLD,
			       "Gravity vector (z): \t\t\t\t%+.3E m/s^2\n",
			       parameter.g[2]);
		  PetscPrintf (PETSC_COMM_WORLD, "\n");

		  break;

		default:

		  PetscPrintf (PETSC_COMM_WORLD,
			       "\nError: Unknown parameter code (%d).\n",
			       tcode);
		  exit (LOGICAL_ERROR);
		  break;

		}

	    }

	}

    }
  while (!feof (fp));

  fclose (fp);

  PetscPrintf (PETSC_COMM_WORLD, "a) Interpolation scheme list\n");
  PetscPrintf (PETSC_COMM_WORLD, "%d-UDS\n", UDS);
  PetscPrintf (PETSC_COMM_WORLD, "%d-CDS\n", CDS);
  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "b) Time advancement method list\n");
  PetscPrintf (PETSC_COMM_WORLD, "%d-Explicit\n", EXPLICITEULER);
  PetscPrintf (PETSC_COMM_WORLD, "%d-Implicit\n", IMPLICITEULER);
  PetscPrintf (PETSC_COMM_WORLD, "%d-Crank Nicolson\n", CRANKNICOLSON);
  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "c) Solvers list\n");
  //PetscPrintf (PETSC_COMM_WORLD, "%d-Jacobi\n", sJACOBI);
  //PetscPrintf (PETSC_COMM_WORLD, "%d-SOR\n", sSOR);
  //PetscPrintf (PETSC_COMM_WORLD, "%d-QMR\n", sQMR);
  PetscPrintf (PETSC_COMM_WORLD, "%d-GMRES\n", sGMRES);
  PetscPrintf (PETSC_COMM_WORLD, "%d-CR\n", sCR);
  PetscPrintf (PETSC_COMM_WORLD, "%d-CG\n", sCG);
  PetscPrintf (PETSC_COMM_WORLD, "%d-CGS\n", sCGS);
  PetscPrintf (PETSC_COMM_WORLD, "%d-BiCG\n", sBICG);
  PetscPrintf (PETSC_COMM_WORLD, "%d-BiCGStab\n", sBICGS);
  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "d) Pre-conditioners list\n");
  PetscPrintf (PETSC_COMM_WORLD, "%d-None\n", pNONE);
  PetscPrintf (PETSC_COMM_WORLD, "%d-Jacobi\n", pJACOBI);
  PetscPrintf (PETSC_COMM_WORLD, "%d-SOR\n", pSOR);
  PetscPrintf (PETSC_COMM_WORLD, "%d-ILU\n", pILU);
  PetscPrintf (PETSC_COMM_WORLD, "%d-ASM\n", pASM);
  PetscPrintf (PETSC_COMM_WORLD, "\n");

  PetscPrintf (PETSC_COMM_WORLD, "Done.\n");

  return LOGICAL_TRUE;

}
