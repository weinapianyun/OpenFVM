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

#include "globals.h"

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"

#include "variables.h"
#include "gradient.h"
#include "geocalc.h"
#include "setup.h"

#include "gamma.h"
#include "velocity.h"
#include "pressure.h"
#include "temperature.h"

#include "solve.h"

void
AllocateMemory () // 申请获取内存
{

  V_Constr (&Co, "Courant number", nbelements, Normal, True);
  V_Constr (&uf, "Face flux velocity", nbfaces, Normal, True);

  V_Constr (&dens, "Density", nbelements, Normal, True);
  V_Constr (&visc, "Dynamic viscosity", nbelements, Normal, True);
  V_Constr (&thcond, "Thermal conductivity", nbelements, Normal, True);
  V_Constr (&spheat, "Specific heat", nbelements, Normal, True);

  V_Constr (&xu0, "Velocity x-component at cell center (previous time step)", nbelements, Normal, True);
  V_Constr (&xv0, "Velocity y-component at cell center (previous time step)", nbelements, Normal, True);
  V_Constr (&xw0, "Velocity z-component at cell center (previous time step)", nbelements, Normal, True);
  V_Constr (&xp0, "Pressure at cell center (previous time step)", nbelements, Normal, True);
  V_Constr (&xT0, "Temperature at cell center (previous time step)", nbelements, Normal, True);
  V_Constr (&xs0, "Gamma at cell center (previous time step)", nbelements, Normal, True);

  V_Constr (&xu, "Velocity x-component at cell center", nbelements, Normal, True);
  V_Constr (&xv, "Velocity y-component at cell center", nbelements, Normal, True);
  V_Constr (&xw, "Velocity z-component at cell center", nbelements, Normal, True);
  V_Constr (&xp, "Pressure at cell center", nbelements, Normal, True);
  V_Constr (&xT, "Temperature at cell center", nbelements, Normal, True);
  V_Constr (&xs, "Gamma at cell center", nbelements, Normal, True);

  V_Constr (&xuf, "Velocity x-component at face center", nbfaces, Normal, True);
  V_Constr (&xvf, "Velocity y-component at face center", nbfaces, Normal, True);
  V_Constr (&xwf, "Velocity z-component at face center", nbfaces, Normal, True);
  V_Constr (&xpf, "Pressure at face center", nbfaces, Normal, True);
  V_Constr (&xTf, "Temperature at face center", nbfaces, Normal, True);
  V_Constr (&xsf, "Gamma at face center", nbfaces, Normal, True);

  V_Constr (&ap, "Momentum matrix diagonal", nbelements, Normal, True);
  V_Constr (&hu, "Momentum matrix source x-component without pressure", nbelements, Normal, True);
  V_Constr (&hv, "Momentum matrix source y-component without pressure", nbelements, Normal, True);
  V_Constr (&hw, "Momentum matrix source z-component without pressure", nbelements, Normal, True);

}

void
DeallocateMemory () // 销毁申请内存
{

  V_Destr (&Co);
  V_Destr (&uf);

  V_Destr (&dens);
  V_Destr (&visc);
  V_Destr (&thcond);
  V_Destr (&spheat);

  V_Destr (&xu0);
  V_Destr (&xv0);
  V_Destr (&xw0);
  V_Destr (&xp0);
  V_Destr (&xT0);
  V_Destr (&xs0);

  V_Destr (&xu);
  V_Destr (&xv);
  V_Destr (&xw);
  V_Destr (&xp);
  V_Destr (&xT);
  V_Destr (&xs);

  V_Destr (&xuf);
  V_Destr (&xvf);
  V_Destr (&xwf);
  V_Destr (&xpf);
  V_Destr (&xTf);
  V_Destr (&xsf);

  V_Destr (&ap);
  V_Destr (&hu);
  V_Destr (&hv);
  V_Destr (&hw);

}

void
CheckMassConservationError (double dt) // 检查质量守恒的误差
{

  int i, j;

  int face; // 单元面的编号，
 // int pair; // 相邻界面的状态

  int element; // 网格单元编号

  double mcp; // 单元速度通量 Uf*A

  double sum; // 统计 Mass conservation error

  sum = 0.0; // 初始化

  for (i = 0; i < nbelements; i++) // 遍历网格单元
  {

      element = i;

      mcp = 0.0; // 单元内初始化

      for (j = 0; j < elements[element].nbfaces; j++) // 遍历单元的界面
	 {
          face = elements[element].face[j];

	      // pair = faces[face].pair;

	      mcp += V_GetCmp (&uf, face + 1) * faces[face].Aj; // mcp = Uf*A

	 }

      sum += LABS (mcp); // 各单元误差累计

  }

  printf ("\nMass conservation error: %+E kg\n", sum);

}

void
Solve (char *var, int *fiter, double dt, double *maxCp, int verbose, int pchecks)
{  /*  SIMPLE算法的一个时间步迭代求解过程  */

  // int i;

  CalculateGamma (var, fiter, dt, maxCp, verbose, pchecks); // 计算单元相函数

  // Set material properties 
  SetMaterialProperties (); // 设置材料参数

  // 系数变量初始化 ap、H
  V_SetAllCmp (&ap, 1.0);

  V_SetAllCmp (&hu, 0.0);
  V_SetAllCmp (&hv, 0.0);
  V_SetAllCmp (&hw, 0.0);

  CalculateVelocity (var, fiter, dt, *maxCp, verbose, pchecks); // 求解速度场

  CalculatePressure (var, fiter, dt, *maxCp, verbose, pchecks); // 求解压力场

  CorrectVelocity (var, fiter, dt, *maxCp, verbose, pchecks); // 修正速度场

  if (pchecks == LOGICAL_TRUE)
    {
      // Check mass conservation
      CheckMassConservationError (dt);
    }

  CalculateTemperature (var, fiter, dt, *maxCp, verbose, pchecks); // 求解温度场

}
