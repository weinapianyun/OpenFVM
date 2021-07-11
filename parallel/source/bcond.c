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
#include "parser.h"
#include "bcond.h"

void
GetEntrySurface (FILE * fp, int j)
{

  int ival;

  gs = calloc (MAXL, sizeof (char));

  fscanf (fp, "%d", &ival);
  bcsurfaces[nbbcsurfaces + j].physreg = ival;

  fscanf (fp, "%s", gs);
  bcsurfaces[nbbcsurfaces + j].fu = calloc (strlen (gs) + 1, sizeof (char));
  strcpy (bcsurfaces[nbbcsurfaces + j].fu, gs);

  fscanf (fp, "%s", gs);
  bcsurfaces[nbbcsurfaces + j].fv = calloc (strlen (gs) + 1, sizeof (char));
  strcpy (bcsurfaces[nbbcsurfaces + j].fv, gs);

  fscanf (fp, "%s", gs);
  bcsurfaces[nbbcsurfaces + j].fw = calloc (strlen (gs) + 1, sizeof (char));
  strcpy (bcsurfaces[nbbcsurfaces + j].fw, gs);

  fscanf (fp, "%s", gs);
  bcsurfaces[nbbcsurfaces + j].fp = calloc (strlen (gs) + 1, sizeof (char));
  strcpy (bcsurfaces[nbbcsurfaces + j].fp, gs);

  fscanf (fp, "%s", gs);
  bcsurfaces[nbbcsurfaces + j].fT = calloc (strlen (gs) + 1, sizeof (char));
  strcpy (bcsurfaces[nbbcsurfaces + j].fT, gs);

  fscanf (fp, "%s", gs);
  bcsurfaces[nbbcsurfaces + j].fs = calloc (strlen (gs) + 1, sizeof (char));
  strcpy (bcsurfaces[nbbcsurfaces + j].fs, gs);

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "Surface: %d, fu: %s\n",
	       bcsurfaces[nbbcsurfaces + j].physreg,
	       bcsurfaces[nbbcsurfaces + j].fu);
  PetscPrintf (PETSC_COMM_WORLD, "Surface: %d, fv: %s\n",
	       bcsurfaces[nbbcsurfaces + j].physreg,
	       bcsurfaces[nbbcsurfaces + j].fv);
  PetscPrintf (PETSC_COMM_WORLD, "Surface: %d, fw: %s\n",
	       bcsurfaces[nbbcsurfaces + j].physreg,
	       bcsurfaces[nbbcsurfaces + j].fw);
  PetscPrintf (PETSC_COMM_WORLD, "Surface: %d, fp: %s\n",
	       bcsurfaces[nbbcsurfaces + j].physreg,
	       bcsurfaces[nbbcsurfaces + j].fp);
  PetscPrintf (PETSC_COMM_WORLD, "Surface: %d, fT: %s\n",
	       bcsurfaces[nbbcsurfaces + j].physreg,
	       bcsurfaces[nbbcsurfaces + j].fT);
  PetscPrintf (PETSC_COMM_WORLD, "Surface: %d, fs: %s\n",
	       bcsurfaces[nbbcsurfaces + j].physreg,
	       bcsurfaces[nbbcsurfaces + j].fs);

  free (gs);

}

void
GetEntryVolume (FILE * fp, int j)
{

  int ival;

  gs = calloc (MAXL, sizeof (char));

  fscanf (fp, "%d", &ival);
  bcvolumes[nbbcvolumes + j].physreg = ival;

  fscanf (fp, "%s", gs);
  bcvolumes[nbbcvolumes + j].fu = calloc (strlen (gs) + 1, sizeof (char));
  strcpy (bcvolumes[nbbcvolumes + j].fu, gs);

  fscanf (fp, "%s", gs);
  bcvolumes[nbbcvolumes + j].fv = calloc (strlen (gs) + 1, sizeof (char));
  strcpy (bcvolumes[nbbcvolumes + j].fv, gs);

  fscanf (fp, "%s", gs);
  bcvolumes[nbbcvolumes + j].fw = calloc (strlen (gs) + 1, sizeof (char));
  strcpy (bcvolumes[nbbcvolumes + j].fw, gs);

  fscanf (fp, "%s", gs);
  bcvolumes[nbbcvolumes + j].fp = calloc (strlen (gs) + 1, sizeof (char));
  strcpy (bcvolumes[nbbcvolumes + j].fp, gs);

  fscanf (fp, "%s", gs);
  bcvolumes[nbbcvolumes + j].fT = calloc (strlen (gs) + 1, sizeof (char));
  strcpy (bcvolumes[nbbcvolumes + j].fT, gs);

  fscanf (fp, "%s", gs);
  bcvolumes[nbbcvolumes + j].fs = calloc (strlen (gs) + 1, sizeof (char));
  strcpy (bcvolumes[nbbcvolumes + j].fs, gs);

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "Volume: %d, fu: %s\n",
	       bcvolumes[nbbcvolumes + j].physreg,
	       bcvolumes[nbbcvolumes + j].fu);
  PetscPrintf (PETSC_COMM_WORLD, "Volume: %d, fv: %s\n",
	       bcvolumes[nbbcvolumes + j].physreg,
	       bcvolumes[nbbcvolumes + j].fv);
  PetscPrintf (PETSC_COMM_WORLD, "Volume: %d, fw: %s\n",
	       bcvolumes[nbbcvolumes + j].physreg,
	       bcvolumes[nbbcvolumes + j].fw);
  PetscPrintf (PETSC_COMM_WORLD, "Volume: %d, fp: %s\n",
	       bcvolumes[nbbcvolumes + j].physreg,
	       bcvolumes[nbbcvolumes + j].fp);
  PetscPrintf (PETSC_COMM_WORLD, "Volume: %d, fT: %s\n",
	       bcvolumes[nbbcvolumes + j].physreg,
	       bcvolumes[nbbcvolumes + j].fT);
  PetscPrintf (PETSC_COMM_WORLD, "Volume: %d, fs: %s\n",
	       bcvolumes[nbbcvolumes + j].physreg,
	       bcvolumes[nbbcvolumes + j].fs);

  free (gs);

}

int
BcdImportBCD (char *file)
{

  int i, j, n;

  int inull;

  int tcode;
  int nbbcd;

  FILE *fp;
  char descr[512];

  fp = fopen (file, "r");

  if (fp == NULL)
    {
      PetscPrintf (PETSC_COMM_WORLD,
		   "\nError: Boundary conditions file not found!\n");
      PetscPrintf (PETSC_COMM_WORLD, "%s\n\n", file);
      exit (LOGICAL_ERROR);
    }

  nbbcsurfaces = 0;
  nbbcvolumes = 0;

  PetscPrintf (PETSC_COMM_WORLD,
	       "\nReading boundary conditions file: %s ...\n", file);

  do
    {

      do
	{
	  fscanf (fp, "%s", descr);

	  if (strcmp (descr, "$Boundary") == 0)
	    break;

	  if (strcmp (descr, "$EndFile") == 0)
	    break;

	}
      while (!feof (fp));

      if (strcmp (descr, "$EndFile") == 0)
	break;

      if (strcmp (descr, "$Boundary") == 0)
	{

	  fscanf (fp, "%d %d ", &inull, &nbbcd);
	  GetLine (fp);

	  for (i = 0; i < nbbcd; i++)
	    {

	      fscanf (fp, "%d %d", &tcode, &n);
	      GetLine (fp);

	      switch (tcode)
		{
		case 10000:

		  // Empty - physical surfaces

		  bcsurfaces =
		    realloc (bcsurfaces,
			     (nbbcsurfaces + n) * sizeof (bcd_surface));

		  for (j = 0; j < n; j++)
		    {
		      bcsurfaces[nbbcsurfaces + j].bc = EMPTY;

		      GetEntrySurface (fp, j);

		    }

		  nbbcsurfaces += n;

		  break;

		case 10020:

		  // Cyclic - physical surfaces

		  bcsurfaces =
		    realloc (bcsurfaces,
			     (nbbcsurfaces + n) * sizeof (bcd_surface));

		  for (j = 0; j < n; j++)
		    {
		      bcsurfaces[nbbcsurfaces + j].bc = CYCLIC;

		      GetEntrySurface (fp, j);

		    }

		  nbbcsurfaces += n;

		  break;

		case 10040:

		  // Permeable - physical surfaces

		  bcsurfaces =
		    realloc (bcsurfaces,
			     (nbbcsurfaces + n) * sizeof (bcd_surface));

		  for (j = 0; j < n; j++)
		    {
		      bcsurfaces[nbbcsurfaces + j].bc = PERMEABLE;

		      GetEntrySurface (fp, j);

		    }

		  nbbcsurfaces += n;

		  break;

		case 10050:

		  // Open - physical surfaces

		  bcsurfaces =
		    realloc (bcsurfaces,
			     (nbbcsurfaces + n) * sizeof (bcd_surface));

		  for (j = 0; j < n; j++)
		    {
		      bcsurfaces[nbbcsurfaces + j].bc = OPEN;

		      GetEntrySurface (fp, j);

		    }

		  nbbcsurfaces += n;

		  break;

		case 10100:

		  // Inlet - physical surfaces

		  bcsurfaces =
		    realloc (bcsurfaces,
			     (nbbcsurfaces + n) * sizeof (bcd_surface));

		  for (j = 0; j < n; j++)
		    {

		      bcsurfaces[nbbcsurfaces + j].bc = INLET;

		      GetEntrySurface (fp, j);

		    }

		  nbbcsurfaces += n;

		  break;

		case 10110:

		  // Pressure inlet - physical surfaces

		  bcsurfaces =
		    realloc (bcsurfaces,
			     (nbbcsurfaces + n) * sizeof (bcd_surface));

		  for (j = 0; j < n; j++)
		    {

		      bcsurfaces[nbbcsurfaces + j].bc = PRESSUREINLET;

		      GetEntrySurface (fp, j);

		    }

		  nbbcsurfaces += n;

		  break;

		case 10150:

		  //  Outlet - physical surfaces

		  bcsurfaces =
		    realloc (bcsurfaces,
			     (nbbcsurfaces + n) * sizeof (bcd_surface));

		  for (j = 0; j < n; j++)
		    {
		      bcsurfaces[nbbcsurfaces + j].bc = OUTLET;

		      GetEntrySurface (fp, j);

		    }

		  nbbcsurfaces += n;

		  break;

		case 10160:

		  //  Slip - physical surfaces

		  bcsurfaces =
		    realloc (bcsurfaces,
			     (nbbcsurfaces + n) * sizeof (bcd_surface));

		  for (j = 0; j < n; j++)
		    {

		      bcsurfaces[nbbcsurfaces + j].bc = SLIP;

		      GetEntrySurface (fp, j);

		    }

		  nbbcsurfaces += n;

		  break;

		case 10170:

		  //  Wall - physical surfaces

		  bcsurfaces =
		    realloc (bcsurfaces,
			     (nbbcsurfaces + n) * sizeof (bcd_surface));

		  for (j = 0; j < n; j++)
		    {

		      bcsurfaces[nbbcsurfaces + j].bc = WALL;

		      GetEntrySurface (fp, j);

		    }

		  nbbcsurfaces += n;

		  break;

		case 10180:

		  //  Adiabtic wall - physical surfaces

		  bcsurfaces =
		    realloc (bcsurfaces,
			     (nbbcsurfaces + n) * sizeof (bcd_surface));

		  for (j = 0; j < n; j++)
		    {
		      bcsurfaces[nbbcsurfaces + j].bc = ADIABATICWALL;

		      GetEntrySurface (fp, j);

		    }

		  nbbcsurfaces += n;

		  break;

		case 10190:

		  //  Moving wall - physical surfaces

		  bcsurfaces =
		    realloc (bcsurfaces,
			     (nbbcsurfaces + n) * sizeof (bcd_surface));

		  for (j = 0; j < n; j++)
		    {

		      bcsurfaces[nbbcsurfaces + j].bc = MOVINGWALL;

		      GetEntrySurface (fp, j);

		    }

		  nbbcsurfaces += n;

		  break;

		case 10200:

		  // Boundary conditions - physical surfaces

		  bcsurfaces =
		    realloc (bcsurfaces,
			     (nbbcsurfaces + n) * sizeof (bcd_surface));

		  for (j = 0; j < n; j++)
		    {

		      bcsurfaces[nbbcsurfaces + j].bc = SURFACE;

		      GetEntrySurface (fp, j);

		    }

		  nbbcsurfaces += n;

		  break;

		case 10250:

		  // Initial conditions - physical volumes

		  bcvolumes =
		    realloc (bcvolumes,
			     (nbbcvolumes + n) * sizeof (bcd_volume));

		  for (j = 0; j < n; j++)
		    {
		      bcvolumes[nbbcvolumes + j].bc = VOLUME;

		      GetEntryVolume (fp, j);
		    }

		  nbbcvolumes += n;

		  break;

		default:

		  PetscPrintf (PETSC_COMM_WORLD,
			       "\nError: Unknown boundary code.\n");
		  exit (LOGICAL_ERROR);
		  break;

		}

	    }

	}

    }
  while (!feof (fp));

  fclose (fp);

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "Number of boundary conditions: \t\t\t%d\n",
	       nbbcsurfaces);
  PetscPrintf (PETSC_COMM_WORLD, "Number of initial conditions: \t\t\t%d\n",
	       nbbcvolumes);

  PetscPrintf (PETSC_COMM_WORLD, "Done.\n");

  return LOGICAL_TRUE;

}
