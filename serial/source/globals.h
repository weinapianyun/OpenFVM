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

#define iu       0		// Velocity - u
#define iv       1		// Velocity - v
#define iw       2		// Velocity - w
#define ip       3		// Pressure
#define iT       4		// Temperature
#define is       5		// Indicator function

#define nphi     6		// Number of design variables

/*
#define idens    0		// Density
#define ivisc    1		// Viscosity
#define ishear   2		// Thermal conductivity
#define ispheat  3		// Specific heat
#define ithcond  4		// Thermal conductivity
#define ipsi     5		// Compressibility

#define nmtl     6		// Number of material properties
*/

#define LOGICAL_TRUE   1
#define LOGICAL_FALSE  0
#define LOGICAL_ERROR -1

#define PI 3.14159265358979323846264

#define LMAX(A,B) ((A)>(B) ? (A):(B))
#define LMIN(A,B) ((A)<(B) ? (A):(B))

#define LABS(X)   ((X) < 0 ? -(X):(X))

#define LSGN(X)   ((X) < 0 ? -(1):(1))

#define SMALL  1E-15
#define VSMALL 1E-300

#define GREAT  1E+15
#define VGREAT 1E+300

#define MAXFACES 6

int verbose;
int pchecks;
