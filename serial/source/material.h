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

typedef struct
{

  float cpsi;

} mat_psi;

typedef struct
{

  float cdens;

} mat_dens;

typedef struct
{

  float cvisc;

} mat_visc;

typedef struct
{

  float cspheat;
  float cthcond;

} mat_therm;

typedef struct
{

  float celastmod;
  float cpoisson;

} mat_mech;

typedef struct
{

  mat_psi psi[2];
  mat_dens dens[2];
  mat_visc visc[2];
  mat_therm therm[2];
  mat_mech mech[2];

  float tens;
  float bthcond;

} mat_material;

mat_material material;

int MtlImportMTL (char *file);
