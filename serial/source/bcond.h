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

enum {
    NONE = 0,
    EMPTY, // 空的
    CYCLIC, // 循环
    PROCESSOR, // 进程边界
    OPEN, // 开放
    INLET, // 入口
    PRESSUREINLET, // 恒压入口
    OUTLET, // 出口
    ADIABATICWALL, // 绝热墙
    MOVINGWALL, // 移动墙
    WALL, // 无滑移
    SLIP, // 滑移
    PERMEABLE, // 可渗透
    CONSTRAINTU, // 速度约束
    CONSTRAINTV,
    CONSTRAINTW,
    CONSTRAINT, // 温度约束
    PRESSURE, // 压力
    SURFACE, // 表面
    VOLUME // 体
} bcd_type;

typedef struct {
    float x, y, z;

} bcd_vector;

typedef struct {

    int physreg;

    int bc;

    char *fu, *fv, *fw;  // 从文件中读取的 字符串类型的 物理场数值
    char *fp;
    char *fT;
    char *fs;

} bcd_surface;

typedef struct {

    int physreg; // 物理编号

    int bc; // 边界条件类型

    char *fu, *fv, *fw; // 从文件中读取的 字符串类型的 物理场数值
    char *fp;
    char *fT;
    char *fs;

} bcd_volume;

// Boundary conditions
int nbbcsurfaces; // 设有边界条件的单元界面的数目
bcd_surface *bcsurfaces; // 指向设有边界条件的网格单元界面

int nbbcvolumes; // 设有边界条件的单元体的数目
bcd_volume *bcvolumes; // 指向设有边界条件的网格单元

int BcdImportBCD(char *file); // 边界条件参数(.bcd文件)导入
