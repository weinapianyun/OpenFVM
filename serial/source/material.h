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

    float cpsi; // 绝对压力

} mat_psi;

typedef struct
{
    float cdens;// 密度

} mat_dens;

typedef struct
{

    float cvisc; // 粘度

} mat_visc;

typedef struct
{
    float cspheat; // 比热容
    float cthcond; // 热导率

} mat_therm;

typedef struct
{
    float celastmod; // 弹性模量
    float cpoisson; // 泊松比

} mat_mech;

typedef struct
{
  // 便于设置初始时已充填的流场, 两相流动
  mat_psi psi[2];
  mat_dens dens[2];
  mat_visc visc[2];
  mat_therm therm[2];
  mat_mech mech[2];

  float tens;
  float bthcond; // 边界的热导率

} mat_material;

mat_material material;

int MtlImportMTL (char *file);
