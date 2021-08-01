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

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"

#include "variables.h"
#include "globals.h"
#include "fill.h"

double
VolumeTotal () // 计算流场的总体积
{
    int i;
    int element;
    double vol; // 流场的总体积

    vol = 0.0; // 初始化
    for (i = 0; i < nbelements; i++) // 遍历所有网格单元
    {
        element = i;
        vol += elements[element].Vp;
    }
    return vol;
}

double
VolumeFilled () // 计算流场已经充填的体积
{

    int i;
    int element;
    double vol; // 流场已经充填的总体积

    vol = 0.0; // 初始化
    for (i = 0; i < nbelements; i++) // 遍历所有网格单元
    {
        element = i;
        // vol = F * Vp , F 为单元的相函数
        vol += V_GetCmp (&xs, element + 1) * elements[element].Vp;
    }
    return vol;
}

double
VolumeEntered (double dt) // 计算一个时间步中, 浇口注入的总体积
{
    int i;
    int face;
    double vol; // 一个时间步dt中, 从浇口中注入流体的总体积

    vol = 0.0;
    for (i = 0; i < nbfaces; i++) // 遍历所有界面
    {
        face = i;
        if (faces[face].bc == INLET) // 寻找浇口界面
        {
            // vol = uf * Aj *dt
            vol += -V_GetCmp (&uf, face + 1) * faces[face].Aj * dt;
        }
    }
    return vol;
}
