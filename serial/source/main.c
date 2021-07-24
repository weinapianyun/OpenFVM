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
#include <malloc.h>
#include <math.h>
#include <time.h>

#ifdef WIN32
#include <windows.h>
#endif

#include "ttime.h"

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"

#include "variables.h"
#include "globals.h"
#include "post.h"
//#include "decomp.h"
#include "reorder.h"
#include "restart.h"
#include "setup.h"
#include "solve.h"
#include "fill.h"

int
Simulation (char *path) // 模拟执行函数
{

  int i, n; // 循环参数

  int d; // 储存 double-字节数

  char var[6]; // 流场变量标志

  double fres[6]; // 存储各变量残差
  int fiter[6]; // 存储各(变量)方程 求解的次数

  int iter; // 外迭代次数

  double starttime, endtime, curtime, dt; // 开始、结束时间，当前时间，时间步
  double wtime, wdt; // 写文件的时刻，写文件的时间间隔

  int irestart; // 重启间隔

  double maxCp; // 最大柯朗数

  double Vc, Vt; // 充填体积、总体积
  double pf; // 充填百分比

  char *file; // 存储输入文件的字符名

  // 文件指针
  FILE *fpresults;		// Output to gmsh post-processing file (results)
  FILE *fpprobe;		// Output to gnuplot file (probe)
  FILE *fpresiduals;	// Output to gnuplot file (residuals)
  FILE *fprestart;		// Input/Output for restart

  var[iu] = 'u'; // 存储字符 'u'
  var[iv] = 'v';
  var[iw] = 'w';
  var[ip] = 'p';
  var[iT] = 'T';
  var[is] = 's';

  // Allocate memory
  printf ("\n");
  printf ("Allocating memory...\n");

  AllocateMemory (); // 分配内存空间

  printf ("Memory allocated.\n");

  // Set initial conditions  设置初始条件
  SetInitialConditions ();

  // Set initial flux  设置初始流量
  SetInitialFlux ();

  // Set boundary velocity and pressure  设置边界条件
  SetBoundary ();

  // Set material properties  设置材料性质参数
  SetMaterialProperties ();

  // Set time intervals  设置时间间隔
  starttime = parameter.t0;
  endtime = parameter.t1;
  curtime = parameter.t0;
  dt = parameter.dt;

  // Set restart intervals
  irestart = parameter.restart; // 设置重新计算的 迭代次数/时间步数

  // Allocate memory for file string
  file = calloc (strlen (path) + 9, sizeof (char)); // 存储 文件的全名

  // Read restart file  读取重启文件
  sprintf (file, "%s.ini", path); // 将字符串"path.ini"格式化输出到 file指向的字符串中

  fprestart = fopen (file, "rb"); // 新建 path.ini二进制文件，将内容清空，以只读形式打开

  if (fprestart != NULL) // 若.ini文件不为空
  {

      ReadRestart (fprestart, &curtime); // 读取重启时间步的流场数据

      endtime += curtime; // 更新终止计算时间

      fclose (fprestart); // 关闭流 fprestart, 刷新所有的缓冲区

  }

  if (parameter.steady == LOGICAL_TRUE) // 若流场流动已经稳定
  {

      sprintf (file, "%s.res", path); // 将字符串 path.res 放到 file 中

      // Open output files for residuals
      fpresiduals = fopen (file, "w"); // 新建 path.res文件，将内容清空，以只写形式打开
  }

  sprintf (file, "%s.pos", path); // 将字符串 path.pos 放到 file 中

  // Open output file for results
  fpresults = fopen (file, "w"); // 新建 path.pos文件，将内容清空，以只写形式打开

  d = sizeof (double); // 存储 double类型的 字节数

  if (parameter.wbinary == LOGICAL_TRUE) // 判断建立的 pos文件是否为二进制类型
  {
      // 将相应的字符串写入 path.pos文件中
      fprintf (fpresults, "$PostFormat\n");
      fprintf (fpresults, "%g %d %d\n", 1.3, 1, d);
      fprintf (fpresults, "$EndPostFormat\n");
  }
  else
  {
      fprintf (fpresults, "$PostFormat\n");
      fprintf (fpresults, "%g %d %d\n", 1.2, 0, d);
      fprintf (fpresults, "$EndPostFormat\n");
  }

  // 生成 gmsh的后处理文件
  WriteResults (fpresults, LOGICAL_TRUE, LOGICAL_TRUE, curtime); // 将初始时刻的结果写入path.pos文件中

  // Set write time intervals  设置写入文件的时间间隔
  if (parameter.nsav != 0) // 若设置的文件保存次数不为 0
      wdt = (endtime - starttime) / parameter.nsav; // 文件保存时间间隔 = (末-初)/保存次数
  else
      wdt = 2 * (endtime - starttime); // 不保存文件

  wtime = wdt; // 写文件的时间 初始化

  iter = 0; // 迭代次数 初始化

  fiter[iu] = 0; // 求解次数初始化
  fiter[iv] = 0;
  fiter[iw] = 0;
  fiter[ip] = 0;
  fiter[iT] = 0;
  fiter[is] = 0;

  if (parameter.fill == LOGICAL_TRUE) // 判断是否计算、输出流场充填百分比
    {
      Vt = VolumeTotal (); // 计算流场单元 总体积
    }

  do
  {

      curtime += dt; // 计算当前时间

      iter++; // 外迭代次数 +1

      printf ("\nTime = %.3E\n", curtime); // 打印当前时间

      Solve (var, fiter, dt, &maxCp, verbose, pchecks); // 迭代求解流场一个时间步的变化

      if (parameter.fill == LOGICAL_TRUE) // 判断是否计算、输出流场充填百分比
	{

	  Vc = VolumeFilled (); // 计算流场充填的体积

	  pf = Vc / Vt * 100; // 计算充填百分比

	  printf ("\nPercentage filled: %.2f%%\n", pf); // 打印 充型体积百分比

	  if (pf > parameter.pf) // 若设定的充填百分比已经达到
	    break; // 结束外迭代

	}

      if (parameter.adjdt == LOGICAL_TRUE) // 判断是否启动时间步长调节
	{

	  if (maxCp > parameter.maxCp) // 若流场单元的最大courant数 > 设定的最大courant数
	    dt *= 0.85 * parameter.maxCp / maxCp; // 修正-减小时间步长

	  if (1.25 * maxCp < parameter.maxCp) // 若流场单元的最大courant数*1.25 < 设定的最大courant数
	    dt *= 1.05; // 修正-增大时间步长

	}

      if (iter >= irestart) // 判断是否到达 设定的重启的迭代次数(时间步数)
	{

	  // Write restart file 写重启文件
	  sprintf (file, "%s.ini", path); // 将字符串 path.ini 放到 file 中

	  fprestart = fopen (file, "wb"); // 新建 path.ini二进制文件，将内容清空，以只写形式打开

	  if (fprestart != NULL) // 若.ini文件不为空
	    {
	      WriteRestart (fprestart, curtime); // 将此刻的计算场景保存至.ini文件中
	    }

	  fclose (fprestart); // 关闭文件

	  irestart += parameter.restart; // 设定新的重启迭代次数

	}

      if (curtime > wtime + dt) // 判断当前时间步 是否计算完毕
	{

	  WriteResults (fpresults, LOGICAL_TRUE, LOGICAL_TRUE, curtime); // 将当前时刻的结果写入path.pos文件中

	  //fflush (fpresults); // 清空缓冲区

	  // Open output file for probes
	  sprintf (file, "%s.prb", path); // 将字符串 path.prb 放到 file 中

	  fpprobe = fopen (file, "w"); // 新建 path.prb文件，将内容清空，以只写形式打开

	  WriteProbeViews (fpprobe, var, curtime); // 在当前时刻创建探测视图，写入path.prb文件中

	  fclose (fpprobe); // 关闭文件

	  wtime += wdt; // 写文件的时刻

	}

      if (parameter.steady == LOGICAL_TRUE) // 判断流场的流动过程是否已经稳定
	{

	  // Get residual 获取各变量的相对偏差
	  if (parameter.calc[iu] == LOGICAL_TRUE) // 若求解时计算 u变量场
	    fres[iu] = l2Norm_V (Sub_VV (&xu, &xu0)); // 计算向量V的 相对偏差范数
	  else
	    fres[iu] = 0.0;

	  if (parameter.calc[iv] == LOGICAL_TRUE)
	    fres[iv] = l2Norm_V (Sub_VV (&xv, &xv0));
	  else
	    fres[iv] = 0.0;

	  if (parameter.calc[iw] == LOGICAL_TRUE)
	    fres[iw] = l2Norm_V (Sub_VV (&xw, &xw0));
	  else
	    fres[iw] = 0.0;

	  if (parameter.calc[ip] == LOGICAL_TRUE)
	    fres[ip] = l2Norm_V (Sub_VV (&xp, &xp0));
	  else
	    fres[ip] = 0.0;

	  if (parameter.calc[iT] == LOGICAL_TRUE)
	    fres[iT] = l2Norm_V (Sub_VV (&xT, &xT0));
	  else
	    fres[iT] = 0.0;

	  if (parameter.calc[is] == LOGICAL_TRUE)
	    fres[is] = l2Norm_V (Sub_VV (&xs, &xs0));
	  else
	    fres[is] = 0.0;

	  n = 0;

	  for (i = 0; i < nphi; i++) // 遍历各变量的相对偏差
	  {

	      if (verbose == LOGICAL_TRUE)
	          printf ("\nVariable: %c Iteration: %d Final residual: %+E\n",
                      var[i], fiter[i], fres[i]); // 打印 求解信息

	      if (fres[i] > parameter.ftol[i]) // 若变量相对偏差大于设定的偏差允许值
	          n++;

	  }

	  if (n == 0)  // 若各变量均满足设定的偏差允许值
	  {
	      printf ("\nSteady state reached.\n");
	      break;
	  }

	  WriteResidual (fpresiduals, iter, fres); // 将当前时刻的相对偏差信息写入 path.res文件中

	  fflush (fpresiduals); // 刷新流 fpresiduals 的输出缓冲区

	}
      else
	{

	  if (curtime + 0.5 * dt > endtime) // 到达结束时间
	    break;

	}

  } while (dt > 0.0);

  WriteResults (fpresults, LOGICAL_TRUE, LOGICAL_TRUE, curtime); // 将结束时刻的结果写入path.pos文件(gmsh后处理)中

  // Close output files
  fclose (fpresults); // 关闭文件

  if (parameter.steady == LOGICAL_TRUE) // 若流场的流动过程已经稳定
  {

      WriteResidual (fpresiduals, iter, fres); // 迭代完毕，将最后的相对偏差写入文件

      fflush (fpresiduals); // 刷新流 fpresiduals 的输出缓冲区

      fclose (fpresiduals); // 关闭文件
  }

  // Open output file for probes
  sprintf (file, "%s.prb", path); // 将字符串 path.prb 放到 file 中

  fpprobe = fopen (file, "w"); // 以只写形式打开 path.prb文件

  WriteProbeViews (fpprobe, var, curtime); // 在结束时刻创建探测视图，写入path.prb文件中

  fclose (fpprobe); // 关闭文件

  // Release memory 
  free (file);

  DeallocateMemory (); // 释放内存空间

  return LOGICAL_TRUE;

}

void
usage () // 使用方法提示
{

  printf ("\n");
  printf ("Usage:     ../OpenFVM [problem] [options] [n]\n");
  printf ("Example 1: ../OpenFVM lid f 1\n");
  printf ("Example 2: ../OpenFVM lid f 1 > lid.log\n");
//  printf ("Example 3: ../OpenFVM lid d 4\n");

  printf ("\n");
  printf ("Main options:\n");
  printf ("  r - Reorder mesh\n");
//  printf ("  d - Decompose the mesh into n regions\n");
  printf ("  f - Start simulation\n");
  printf ("\n");

  printf ("Additional options:\n");
  printf ("  v - Verbose mode\n");
  printf ("  c - Perform checks\n");
  printf ("\n");

}

int
main (int argc, char **argv) // argc为输入文件的数目, argv为指令流/指针数组
{ // 其中，argv[0~3]对应的指针指向的字符串分别为：  path\to\flow case r 1

  char *path;
  char *file;

  double start, end; // 模拟开始、结束时间
  double elapsed; // 模拟求解的时间

  printf ("\n");
  printf ("*****************************************\n");
  printf ("*                                       *\n");
  printf ("*    OpenFVM-Flow v1.0 - Serial         *\n");
  printf ("*                                       *\n");
  printf ("*****************************************\n");
  printf ("\n");

#ifdef WIN32
  SetThreadPriority (GetCurrentThread (), THREAD_PRIORITY_BELOW_NORMAL);
#endif

  if (argc != 4) // 检查输入参数是否为 4 个字符串
    {
      usage ();
      printf ("\nError: Wrong number of arguments, expected %d found %d.\n\n",
	      4 - 1, argc - 1);
      return LOGICAL_ERROR;
    }

  /* 判断输入的执行命令的特征字符 —— c v f r d  */

  if (strchr (argv[2], 'c') != NULL)
    pchecks = LOGICAL_TRUE; // 将 pchecks 设置为逻辑真
  else
    pchecks = LOGICAL_FALSE;

  if (strchr (argv[2], 'v') != NULL)
    verbose = LOGICAL_TRUE;
  else
    verbose = LOGICAL_FALSE;

  // 存储 4 个输入文件的目录名  eg: djc
  path = calloc (strlen (argv[1]) + 1, sizeof (char));
  // 存储 4 个输入文件的 完整文件名  eg: djc.msh / djc.bcd / djc.par /djc.mtl
  file = calloc (strlen (argv[1]) + 12, sizeof (char));

  // Simulate
  if (strchr (argv[2], 'f') != NULL) // 模拟执行指令

    {

      strcpy (path, argv[1]);

      // Read parameter file
      sprintf (file, "%s.par", path);
      ParImportPAR (file); // 读取 file 中字符串代表的 par文件

      // Read mesh file
      sprintf (file, "%s.msh", path);
      MshImportMSH (file);

      // Read boundary conditions file
      sprintf (file, "%s.bcd", path);
      BcdImportBCD (file);

      // Read material file
      sprintf (file, "%s.mtl", path);
      MtlImportMTL (file);

      // Start clock to measure simulation time
      start = ttime ();

      // Start simulation
      Simulation (path);

      end = ttime ();
      elapsed = end - start;
      printf ("\nElapsed time: %.2f seconds.\n", elapsed);

    }

  // Decompose
  /*
  if (strchr (argv[2], 'd') != NULL) // 网格分区指令
    {

      strcpy (path, argv[1]);

      // Read mesh file
      sprintf (file, "%s.msh", path);

      if (MshImportMSH (file) == LOGICAL_ERROR)
	{
	  printf ("\nError: Mesh file not found!\n");
	  printf ("%s\n\n", file);
	  return LOGICAL_ERROR;
	}

      // Allocate memory
      printf ("\n");
      printf ("Allocating memory...\n");

      AllocateMemory ();

      printf ("Memory allocated.\n");

      SetBoundary ();

      DeallocateMemory ();

      DecomposeMesh (path, atoi (argv[3]));

    }
  */

  // Reorder
  if (strchr (argv[2], 'r') != NULL) // 网格单元重排序指令

    {

      strcpy (path, argv[1]);

      // Read mesh file
      sprintf (file, "%s.msh", path);

      if (MshImportMSH (file) == LOGICAL_ERROR)
	{
	  printf ("\nError: Mesh file not found!\n");
	  printf ("%s\n\n", file);
	  return LOGICAL_ERROR;
	}

      ReorderMesh (path);

    }

  // 释放内存
  free (path);
  free (file);

  // Free memory
  MshFreeMemory ();

  printf ("Done.\n\n");

  return 0;

}
