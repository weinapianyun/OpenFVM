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

#include <time.h>

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

void
WriteVector (Vec * x)
{

  PetscViewer view;

  PetscViewerASCIIOpen (PETSC_COMM_WORLD, "vector.m", &view);

  PetscViewerSetFormat (view, PETSC_VIEWER_ASCII_MATLAB);

  VecView (*x, view);

  PetscViewerDestroy (view);

}

void
WriteMatrix (Mat * A)
{

  PetscViewer view;

  PetscViewerASCIIOpen (PETSC_COMM_WORLD, "matrix.m", &view);

  PetscViewerSetFormat (view, PETSC_VIEWER_ASCII_MATLAB);

  MatView (*A, view);

  PetscViewerDestroy (view);

}

int
CheckIfDiagonalMatrix (Mat * A)
{

  int i, j;

  int nr, nc;

  int cond1, cond2;

  double sum_ap_parcial, sum_an_parcial;
  double sum_ap_total, sum_an_total;

  MatGetSize (*A, &nr, &nc);

  sum_ap_total = 0.0;
  sum_an_total = 0.0;

  cond1 = LOGICAL_FALSE;
  cond2 = LOGICAL_FALSE;

  sum_ap_parcial = 0.0;
  
  for (i = 0; i < nr; i++)
    {

      if (nc > 0)
	sum_ap_parcial = LABS (Q_GetVal (A, i, i));

      sum_an_parcial = 0.0;

      for (j = 0; j < nc; j++)
	{

	  if (i == j)
	    continue;

	  sum_an_parcial += LABS (Q_GetVal (A, i, j));

	}

      if (sum_ap_parcial >= sum_an_parcial)
	cond1 = LOGICAL_TRUE;

      sum_ap_total += sum_ap_parcial;

      sum_an_total += sum_an_parcial;

    }

  if (sum_ap_total >= sum_an_total)
    cond2 = LOGICAL_TRUE;

  if (!cond1 && !cond2)
    return LOGICAL_FALSE;
  else
    return LOGICAL_TRUE;

}

void
SolveMatrix (Mat * A, Vec * x, Vec * b, int *iter, double *res, double *mtime, int msolver, int mprecond, int miter, float mtol)
{

  double start, end;

  PCType PrecondProc;
  
  start = MPI_Wtime ();

  // (0-None, 1-Jacobi, 2-SOR, 3-ILU, 4-ASM)

  PrecondProc = 0;

  switch (mprecond)
    {
    case pNONE:

      PrecondProc = PCNONE;
      break;

    case pJACOBI:

      PrecondProc = PCJACOBI;
      break;

    case pSOR:

      PrecondProc = PCSOR;
      break;

    case pILU:

      PrecondProc = PCILU;
      break;

    case pASM:

      PrecondProc = PCASM;
      break;

    }

  // (3-GMRES, 4-CR, 5-CG, 6-CGS, 7-BiCG, 8-BiCGStab)
  switch (msolver)
    {
    case sGMRES:

      // GMRES
      GMRESIter (A, x, b, miter, iter, mtol, res, PrecondProc);
      break;

    case sCR:

      // CR
      CRIter (A, x, b, miter, iter, mtol, res, PrecondProc);
      break;

    case sCG:

      // CG
      CGIter (A, x, b, miter, iter, mtol, res, PrecondProc);
      break;

    case sCGS:

      // CGS
      CGSIter (A, x, b, miter, iter, mtol, res, PrecondProc);
      break;

    case sBICG:

      // BiCG
      BiCGIter (A, x, b, miter, iter, mtol, res, PrecondProc);
      break;

    case sBICGS:

      // BiCGSTAB
      BiCGSTABIter (A, x, b, miter, iter, mtol, res, PrecondProc);
      break;

    default:

      PetscPrintf (PETSC_COMM_WORLD,
		   "\nError: Solver method (%d) not implemented\n", msolver);
      exit (LOGICAL_ERROR);
      break;

    }

  end = MPI_Wtime ();

  *mtime = end - start;

  VecGhostUpdateBegin (*x, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd (*x, INSERT_VALUES, SCATTER_FORWARD);

}
