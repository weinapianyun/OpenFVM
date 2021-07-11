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
#include <malloc.h>
#include <math.h>
#include <time.h>

#ifdef WIN32
#include <windows.h>
#endif

#include "ttime.h"

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"

#include "variables.h"
#include "globals.h"
#include "post.h"
//#include "decomp.h"
#include "reorder.h"
#include "restart.h"
#include "setup.h"
#include "solve.h"
#include "fill.h"

int
Simulation (char *path)
{

  int i, n;

  int d;

  char var[6];

  double fres[6];
  int fiter[6];

  int iter;

  double starttime, endtime, curtime, dt;
  double wtime, wdt;

  int irestart;

  double maxCp;

  double Vc, Vt;
  double pf;

  char *file;

  FILE *fpresults;		// Output to gmsh post-processing file (results
  FILE *fpprobe;		// Output to gnuplot file (probe)
  FILE *fpresiduals;		// Output to gnuplot file (residuals)
  FILE *fprestart;		// Input/Output for restart

  var[iu] = 'u';
  var[iv] = 'v';
  var[iw] = 'w';
  var[ip] = 'p';
  var[iT] = 'T';
  var[is] = 's';

  // Allocate memory
  printf ("\n");
  printf ("Allocating memory...\n");

  AllocateMemory ();

  printf ("Memory allocated.\n");

  // Set initial conditions
  SetInitialConditions ();

  // Set initial flux
  SetInitialFlux ();

  // Set boundary velocity and pressure
  SetBoundary ();

  // Set material properties 
  SetMaterialProperties ();

  // Set time intervals
  starttime = parameter.t0;
  endtime = parameter.t1;
  curtime = parameter.t0;
  dt = parameter.dt;

  // Set restart intervals
  irestart = parameter.restart;

  // Allocate memory for file string
  file = calloc (strlen (path) + 9, sizeof (char));

  // Read restart file
  sprintf (file, "%s.ini", path);

  fprestart = fopen (file, "rb");

  if (fprestart != NULL)
    {

      ReadRestart (fprestart, &curtime);

      endtime += curtime;

      fclose (fprestart);

    }

  if (parameter.steady == LOGICAL_TRUE)
    {

      sprintf (file, "%s.res", path);

      // Open output files for residuals
      fpresiduals = fopen (file, "w");
    }

  sprintf (file, "%s.pos", path);

  // Open output file for results
  fpresults = fopen (file, "w");

  d = sizeof (double);

  if (parameter.wbinary == LOGICAL_TRUE)
    {
      fprintf (fpresults, "$PostFormat\n");
      fprintf (fpresults, "%g %d %d\n", 1.3, 1, d);
      fprintf (fpresults, "$EndPostFormat\n");
    }
  else
    {
      fprintf (fpresults, "$PostFormat\n");
      fprintf (fpresults, "%g %d %d\n", 1.2, 0, d);
      fprintf (fpresults, "$EndPostFormat\n");
    }

  WriteResults (fpresults, LOGICAL_TRUE, LOGICAL_TRUE, curtime);

  // Set write time intervals
  if (parameter.nsav != 0)
    wdt = (endtime - starttime) / parameter.nsav;
  else
    wdt = 2 * (endtime - starttime);

  wtime = wdt;

  iter = 0;

  fiter[iu] = 0;
  fiter[iv] = 0;
  fiter[iw] = 0;
  fiter[ip] = 0;
  fiter[iT] = 0;
  fiter[is] = 0;

  if (parameter.fill == LOGICAL_TRUE)
    {
      Vt = VolumeTotal ();
    }

  do
    {

      curtime += dt;

      iter++;

      printf ("\nTime = %.3E\n", curtime);

      Solve (var, fiter, dt, &maxCp, verbose, pchecks);

      if (parameter.fill == LOGICAL_TRUE)
	{

	  Vc = VolumeFilled ();

	  pf = Vc / Vt * 100;

	  printf ("\nPercentage filled: %.2f%%\n", pf);

	  if (pf > parameter.pf)
	    break;

	}

      if (parameter.adjdt == LOGICAL_TRUE)
	{

	  if (maxCp > parameter.maxCp)
	    dt *= 0.85 * parameter.maxCp / maxCp;

	  if (1.25 * maxCp < parameter.maxCp)
	    dt *= 1.05;

	}

      if (iter >= irestart)
	{

	  // Write restart file
	  sprintf (file, "%s.ini", path);

	  fprestart = fopen (file, "wb");

	  if (fprestart != NULL)
	    {
	      WriteRestart (fprestart, curtime);
	    }

	  fclose (fprestart);

	  irestart += parameter.restart;

	}

      if (curtime > wtime + dt)
	{

	  WriteResults (fpresults, LOGICAL_TRUE, LOGICAL_TRUE, curtime);

	  //fflush (fpresults);

	  // Open output file for probes          
	  sprintf (file, "%s.prb", path);

	  fpprobe = fopen (file, "w");

	  WriteProbeViews (fpprobe, var, curtime);

	  fclose (fpprobe);

	  wtime += wdt;

	}

      if (parameter.steady == LOGICAL_TRUE)
	{

	  // Get residual
	  if (parameter.calc[iu] == LOGICAL_TRUE)
	    fres[iu] = l2Norm_V (Sub_VV (&xu, &xu0));
	  else
	    fres[iu] = 0.0;

	  if (parameter.calc[iv] == LOGICAL_TRUE)
	    fres[iv] = l2Norm_V (Sub_VV (&xv, &xv0));
	  else
	    fres[iv] = 0.0;

	  if (parameter.calc[iw] == LOGICAL_TRUE)
	    fres[iw] = l2Norm_V (Sub_VV (&xw, &xw0));
	  else
	    fres[iw] = 0.0;

	  if (parameter.calc[ip] == LOGICAL_TRUE)
	    fres[ip] = l2Norm_V (Sub_VV (&xp, &xp0));
	  else
	    fres[ip] = 0.0;

	  if (parameter.calc[iT] == LOGICAL_TRUE)
	    fres[iT] = l2Norm_V (Sub_VV (&xT, &xT0));
	  else
	    fres[iT] = 0.0;

	  if (parameter.calc[is] == LOGICAL_TRUE)
	    fres[is] = l2Norm_V (Sub_VV (&xs, &xs0));
	  else
	    fres[is] = 0.0;

	  n = 0;

	  for (i = 0; i < nphi; i++)
	    {

	      if (verbose == LOGICAL_TRUE)
		printf ("\nVariable: %c Iteration: %d Final residual: %+E\n",
			var[i], fiter[i], fres[i]);

	      if (fres[i] > parameter.ftol[i])
		n++;

	    }

	  if (n == 0)
	    {
	      printf ("\nSteady state reached.\n");
	      break;
	    }

	  WriteResidual (fpresiduals, iter, fres);

	  fflush (fpresiduals);

	}
      else
	{

	  if (curtime + 0.5 * dt > endtime)
	    break;

	}

    }
  while (dt > 0.0);

  WriteResults (fpresults, LOGICAL_TRUE, LOGICAL_TRUE, curtime);

  // Close output files
  fclose (fpresults);

  if (parameter.steady == LOGICAL_TRUE)
    {

      WriteResidual (fpresiduals, iter, fres);

      fflush (fpresiduals);

      fclose (fpresiduals);
    }

  // Open output file for probes
  sprintf (file, "%s.prb", path);

  fpprobe = fopen (file, "w");

  WriteProbeViews (fpprobe, var, curtime);

  fclose (fpprobe);

  // Release memory 
  free (file);

  DeallocateMemory ();

  return LOGICAL_TRUE;

}

void
usage ()
{

  printf ("\n");
  printf ("Usage:     ../OpenFVM [problem] [options] [n]\n");
  printf ("Example 1: ../OpenFVM lid f 1\n");
  printf ("Example 2: ../OpenFVM lid f 1 > lid.log\n");
//  printf ("Example 3: ../OpenFVM lid d 4\n");

  printf ("\n");
  printf ("Main options:\n");
  printf ("  r - Reorder mesh\n");
//  printf ("  d - Decompose the mesh into n regions\n");
  printf ("  f - Start simulation\n");
  printf ("\n");

  printf ("Additional options:\n");
  printf ("  v - Verbose mode\n");
  printf ("  c - Perform checks\n");
  printf ("\n");

}

int
main (int argc, char **argv)
{

  char *path;
  char *file;

  double start, end;
  double elapsed;

  printf ("\n");
  printf ("*****************************************\n");
  printf ("*                                       *\n");
  printf ("*    OpenFVM-Flow v1.0 - Serial         *\n");
  printf ("*                                       *\n");
  printf ("*****************************************\n");
  printf ("\n");

#ifdef WIN32
  SetThreadPriority (GetCurrentThread (), THREAD_PRIORITY_BELOW_NORMAL);
#endif

  if (argc != 4)
    {
      usage ();
      printf ("\nError: Wrong number of arguments, expected %d found %d.\n\n",
	      4 - 1, argc - 1);
      return LOGICAL_ERROR;
    }

  if (strchr (argv[2], 'c') != NULL)
    pchecks = LOGICAL_TRUE;
  else
    pchecks = LOGICAL_FALSE;

  if (strchr (argv[2], 'v') != NULL)
    verbose = LOGICAL_TRUE;
  else
    verbose = LOGICAL_FALSE;

  path = calloc (strlen (argv[1]) + 1, sizeof (char));
  file = calloc (strlen (argv[1]) + 12, sizeof (char));

  // Simulate
  if (strchr (argv[2], 'f') != NULL)
    {

      strcpy (path, argv[1]);

      // Read parameter file
      sprintf (file, "%s.par", path);
      ParImportPAR (file);

      // Read mesh file
      sprintf (file, "%s.msh", path);
      MshImportMSH (file);

      // Read boundary conditions file
      sprintf (file, "%s.bcd", path);
      BcdImportBCD (file);

      // Read material file
      sprintf (file, "%s.mtl", path);
      MtlImportMTL (file);

      // Start clock to measure simulation time
      start = ttime ();

      // Start simulation
      Simulation (path);

      end = ttime ();
      elapsed = end - start;
      printf ("\nElapsed time: %.2f seconds.\n", elapsed);

    }

  // Decompose
  /*
  if (strchr (argv[2], 'd') != NULL)
    {

      strcpy (path, argv[1]);

      // Read mesh file
      sprintf (file, "%s.msh", path);

      if (MshImportMSH (file) == LOGICAL_ERROR)
	{
	  printf ("\nError: Mesh file not found!\n");
	  printf ("%s\n\n", file);
	  return LOGICAL_ERROR;
	}

      // Allocate memory
      printf ("\n");
      printf ("Allocating memory...\n");

      AllocateMemory ();

      printf ("Memory allocated.\n");

      SetBoundary ();

      DeallocateMemory ();

      DecomposeMesh (path, atoi (argv[3]));

    }
  */

  // Reorder
  if (strchr (argv[2], 'r') != NULL)
    {

      strcpy (path, argv[1]);

      // Read mesh file
      sprintf (file, "%s.msh", path);

      if (MshImportMSH (file) == LOGICAL_ERROR)
	{
	  printf ("\nError: Mesh file not found!\n");
	  printf ("%s\n\n", file);
	  return LOGICAL_ERROR;
	}

      ReorderMesh (path);

    }

  free (path);
  free (file);

  // Free memory          
  MshFreeMemory ();

  printf ("Done.\n\n");

  return 0;

}
