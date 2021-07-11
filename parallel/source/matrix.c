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

// PETSc 
#include "petscmat.h"

#include "globals.h"
#include "mesh.h"

#include "matrix.h"

void
Q_Constr (Mat *A, int n, int symmetric)
{

  if (symmetric == LOGICAL_TRUE)
    {
      MatCreateMPISBAIJ (PETSC_COMM_WORLD, 1, n, n, PETSC_DETERMINE,
			 PETSC_DETERMINE, MAXFACES + 1, PETSC_NULL,
			 MAXFACES + 1, PETSC_NULL, A);
      MatSetOption (*A, MAT_SYMMETRIC);
      MatSetOption (*A, MAT_IGNORE_LOWER_TRIANGULAR);
    }
  else
    {
      MatCreateMPIBAIJ (PETSC_COMM_WORLD, 1, n, n, PETSC_DETERMINE,
			PETSC_DETERMINE, MAXFACES + 1, PETSC_NULL,
			MAXFACES + 1, PETSC_NULL, A);
    }

  //MatCreateMPIAIJ (PETSC_COMM_WORLD, n, n, PETSC_DETERMINE, PETSC_DETERMINE, MAXFACES, PETSC_NULL, MAXFACES, PETSC_NULL, A); 

}

void
Q_Destr (Mat *A)
{

  MatDestroy (*A);

}

double
Q_GetVal (Mat *A, int row, int col)
{

  double value;

  MatGetValues (*A, 1, &row, 1, &col, &value);

  return value;


}

void
Q_SetEntry (Mat *A, int row, int col, double value)
{

  MatSetValue (*A, row, col, value, INSERT_VALUES);

}

void
Q_SetEntries (Mat *A, int nrows, int *row, int ncols, int *col, double *values)
{

  MatSetValues (*A, nrows, row, ncols, col, values, INSERT_VALUES);

}

