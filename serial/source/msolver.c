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

#include "ttime.h"

#include "globals.h"

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"

#include "variables.h"
#include "gradient.h"
#include "geocalc.h"
#include "setup.h"
#include "ttime.h"

#include "msolver.h"

void
PrintVector (Vector * x)
{

  int i;

  int nr;

  nr = V_GetDim (x);

  for (i = 0; i < nr; i++)
    {
      printf ("%.16e\n", V_GetCmp (x, i + 1));
    }

}

void
WriteVector (Vector * x)
{

  FILE *fp;

  int i;

  int nr;

  nr = V_GetDim (x);

  fp = fopen ("vector.m", "w");

  fprintf (fp, "Vec_s = [\n");

  for (i = 0; i < nr; i++)
    {
      fprintf (fp, "%.16e\n", V_GetCmp (x, i + 1));
    }

  fprintf (fp, "];\n");

  fclose (fp);

}

void
WriteMatrix (QMatrix * A, int symmetric)
{

  FILE *fp;

  int i, j;

  int nr;
  int nj;

  int nonzeros;

  nr = Q_GetDim (A);

  fp = fopen ("matrix.m", "w");

  nonzeros = 0;

  for (i = 0; i < nr; i++)
    {

      nj = Q_GetLen (A, i + 1);

      nonzeros += nj;
    }


  fprintf (fp, "%% Size = %d %d\n", nr, nr);
  fprintf (fp, "%% Nonzeros = %d\n", nonzeros);
  fprintf (fp, "zzz = zeros(%d,3);\n", nonzeros);
  fprintf (fp, "zzz = [\n");

  for (i = 0; i < nr; i++)
    {

      nj = Q_GetLen (A, i + 1);

      for (j = 0; j < nj; j++)
	{
	  fprintf (fp, "%d %d  %.16e\n", i + 1, Q_GetPos (A, i + 1, j),
		   Q_GetVal (A, i + 1, j));
	}
    }

  fprintf (fp, "];\n");
  fprintf (fp, "Mat_s = spconvert(zzz);\n");

  if (symmetric == LOGICAL_TRUE)
    fprintf (fp, "Mat_s = triu(Mat_s, 1)+tril(Mat_s', 0);\n");

}

int
CheckIfDiagonalMatrix (QMatrix * A)
{

  int i, j;

  int nr;

  int nj;

  int cond1, cond2;

  double sum_ap_parcial, sum_an_parcial;
  double sum_ap_total, sum_an_total;

  sum_ap_total = 0.0;
  sum_an_total = 0.0;

  cond1 = LOGICAL_FALSE;
  cond2 = LOGICAL_FALSE;
  nr = Q_GetDim (A);
  
  sum_ap_parcial = 0.0;
  
  for (i = 0; i < nr; i++)
    {

      nj = Q_GetLen (A, i + 1);

      if (nj > 0)
	sum_ap_parcial = LABS (Q_GetVal (A, i + 1, 0));

      sum_an_parcial = 0.0;

      for (j = 1; j < nj; j++)
	{
	  sum_an_parcial += LABS (Q_GetVal (A, i + 1, j));
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
SolveMatrix (QMatrix * A, Vector * x, Vector * b, int *iter, double *res,
	     double *time, int msolver, int mprecond, int miter)
{

  PrecondProcType PrecondProc;

  double start, end;

  start = ttime ();

  // (0-Null, 1-Jacobi, 2-SOR, 3-ILU)

  PrecondProc = 0;
 
  switch (mprecond)
    {
    case pNONE:

      PrecondProc = NULL;
      break;

    case pJACOBI:

      PrecondProc = JacobiPrecond;
      break;

    case pSOR:

      PrecondProc = SSORPrecond;
      break;

    case pILU:

      PrecondProc = ILUPrecond;
      break;

    case pASM:

      PrecondProc = ILUPrecond;
      break;

    }

  // (0-Jacobi, 1-SOR, 2-QMR, 3-GMRES, 4-CG, 5-CGN, 6-CGS, 7-BiCG, 8-BiCGStab) 

  switch (msolver)
    {
    case sJACOBI:

      // Jacobi
      JacobiIter (A, x, b, miter, PrecondProc, 1.0);

      *iter = GetLastNoIter ();
      *res = GetLastAccuracy ();
      break;

    case sSOR:

      // SOR
      SSORIter (A, x, b, miter, PrecondProc, 1.0);

      *iter = GetLastNoIter ();
      *res = GetLastAccuracy ();
      break;

    case sQMR:

      // QMR
      QMRIter (A, x, b, miter, PrecondProc, 1.0);

      *iter = GetLastNoIter ();
      *res = GetLastAccuracy ();
      break;

    case sGMRES:

      // GMRES
      GMRESIter (A, x, b, miter, PrecondProc, 1.0);

      *iter = GetLastNoIter ();
      *res = GetLastAccuracy ();
      break;

    case sCG:

      // CG
      CGIter (A, x, b, miter, PrecondProc, 1.0);

      *iter = GetLastNoIter ();
      *res = GetLastAccuracy ();
      break;

    case sCGN:

      // CGN
      CGNIter (A, x, b, miter, PrecondProc, 1.0);

      *iter = GetLastNoIter ();
      *res = GetLastAccuracy ();
      break;

    case sCGS:

      // CGS
      CGSIter (A, x, b, miter, PrecondProc, 1.0);

      *iter = GetLastNoIter ();
      *res = GetLastAccuracy ();
      break;

    case sBICG:

      // BICG
      BiCGIter (A, x, b, miter, PrecondProc, 1.0);

      *iter = GetLastNoIter ();
      *res = GetLastAccuracy ();
      break;

    case sBICGS:

      // BICGS
      BiCGSTABIter (A, x, b, miter, PrecondProc, 1.0);

      *iter = GetLastNoIter ();
      *res = GetLastAccuracy ();
      break;

    default:

      printf ("\nError: Solver method (%d) not implemented\n", msolver);
      exit (LOGICAL_ERROR);
      break;

    }

  end = ttime ();

  *time = end - start;

}
