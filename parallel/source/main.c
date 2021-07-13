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

#include <stdio.h> // 默认IO库
#include <stdlib.h> //基本函数库
#include <string.h> //字符串处理
#include <malloc.h> //动态内存
#include <math.h> //普通数学计算
#include <time.h> //时间信息处理

#include "variables.h" //引用PETSc库的类型建立的变量
#include "vector.h" //计算的向量

#include "mesh.h" //网格
#include "material.h" //材质
#include "bcond.h" //边界条件
#include "param.h" //计算参数

#include "globals.h" //全局变量
#include "geocalc.h" //物理向量的计算
#include "parser.h" //？？
#include "post.h" //后处理
#include "decomp.h" //网格划分
#include "reorder.h" //节点重排序
#include "parallel.h" //并行计算
#include "restart.h" //重新计算
#include "setup.h" //初始计算条件设置
#include "solve.h" //控制内存、计算过程
#include "fill.h" //单元充型计算

int
Simulation (char *path) //仿真过程
{

  int i, n;

  int d;

  char var[6];

  double fres[6];

  int fiter[6];

  int iter; //迭代次数

  int irestart;

  double starttime, endtime, curtime, dt; //开始、结束时间，当前时间，时间步长
  double wtime, wdt;

  double maxCp;

  double Vc, Vt;
  double Vcp, Vtp;

  double pf;

  char *file;
  // 输出文件
  FILE *fpresults;		// Output to gmsh post-processing file (results)
  FILE *fpprobe;		// Output to gnuplot file (probe)
  FILE *fpresiduals;		// Output to gnuplot file (residuals)
  FILE *fprestart;		// Input/Output for restart 

  var[iu] = 'u'; //储存物理量
  var[iv] = 'v';
  var[iw] = 'w';
  var[ip] = 'p';
  var[iT] = 'T';
  var[is] = 's';

  // Set ghosts 影子单元
  SetGhosts ();

  // Allocate memory  分配内存
  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "Allocating memory...\n"); 

  AllocateMemory ();

  PetscPrintf (PETSC_COMM_WORLD, "Memory allocated.\n");

  // Set elements center 设置单元中心
  SetCenters ();

  VecGhostGetLocalForm (cex, &cexl);
  VecGhostGetLocalForm (cey, &ceyl);
  VecGhostGetLocalForm (cez, &cezl);

  // Set initial conditions  设置初始条件
  SetInitialConditions ();

  // Set initial flux  设置初始流量
  SetInitialFlux ();

  // Set boundary velocity and pressure  设置初始速度和压力
  SetBoundary ();

  // Set material properties    设置材料属性
  SetMaterialProperties ();

  // Set time intervals 设置时间间隔
  starttime = parameter.t0;
  endtime = parameter.t1;
  curtime = parameter.t0;
  dt = parameter.dt;

  // Set restart intervals 设置重启间隔
  irestart = parameter.restart;

  // Allocate memory for file string  为文件字符串分配内存
  file = calloc (strlen (path) + 9, sizeof (char));

  // Read restart file  读取重启文件
  sprintf (file, "%s.%03d.ini", path, processor);

  fprestart = fopen (file, "rb");

  if (fprestart != NULL)
    {

      ReadRestart (fprestart, &curtime);

      endtime += curtime;

      fclose (fprestart);

    }

  if (parameter.steady == LOGICAL_TRUE)
    {

      sprintf (file, "%s.res", path);

      // Open output files for residuals  打开输出文件查找残差
      PetscFOpen (PETSC_COMM_WORLD, file, "w", &fpresiduals);
    }

  sprintf (file, "%s.%03d.pos", path, processor);

  // Open output file for results  打开输出文件获取结果
  fpresults = fopen (file, "w");

  d = sizeof (double);

  if (parameter.wbinary == LOGICAL_TRUE)
    {
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

  WriteResults (fpresults, LOGICAL_TRUE, LOGICAL_TRUE, curtime);

  // Set write time intervals  设置写入时间间隔
  if (parameter.nsav != 0)
    wdt = (endtime - starttime) / parameter.nsav;
  else
    wdt = 2 * (endtime - starttime);

  wtime = wdt;

  iter = 0;

  fiter[iu] = 0;
  fiter[iv] = 0;
  fiter[iw] = 0;
  fiter[ip] = 0;
  fiter[iT] = 0;
  fiter[is] = 0;

  if (parameter.fill == LOGICAL_TRUE) //判断求解区域的性质
    {
      Vtp = VolumeTotal (); //计算单个单元总体积

      Vc = 0.0;

      Vt = Vtp;  //单元体积的累加

      if (processor == 0) //判断当前的进程标号是否为 0
	{
	  for (i = 1; i < nbprocessors; i++) //对除了0号之外的所有进程操作
	    {
	      // 0号进程(可能是主进程) 接收一个double类型的数据-单元的体积
	      MPI_Recv (&Vtp, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, 0);
	      //函数的参数依次为 0进程 接收缓冲区的起始地址、最多可接收的数据的个数、接收数据的数据类型
	      //发送数据的进程的进程标识号、消息标识、本进程和发送进程所在的通信域、返回状态

	      Vt += Vtp;
	    }
	}
      else
	{
          //发送一个double类型的数据 单元的体积
	  MPI_Send (&Vtp, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}

      MPI_Barrier (MPI_COMM_WORLD); //保持各进程的同步
    }

  do
    {

      curtime += dt;

      iter++;

      PetscPrintf (PETSC_COMM_WORLD, "\nTime = %.3E\n", curtime);

      Solve (var, fiter, dt, &maxCp, verbose, pchecks);

      if (parameter.fill == LOGICAL_TRUE)
	{

	  Vcp = VolumeFilled ();

	  Vc = Vcp;

	  if (processor == 0)
	    {

	      for (i = 1; i < nbprocessors; i++)
		{
		  MPI_Recv (&Vcp, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, 0);// 0号进程接收一个double类型的数据-单元p已充填的体积

		  Vc += Vcp;
		}

	    }
	  else
	    {
	      MPI_Send (&Vcp, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);//其他进程发送一个double类型的数据 单元p已充填的体积
	    }

	  MPI_Barrier (MPI_COMM_WORLD);

	  if (processor == 0)
	    {
	      pf = Vc / Vt * 100;

	      for (i = 1; i < nbprocessors; i++)
		{
		  MPI_Send (&pf, 1, MPI_DOUBLE, i, 4, MPI_COMM_WORLD); //进程0 向其他进程发送单元界面上的压力值
		}
	    }
	  else
	    {
	      MPI_Recv (&pf, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, 0);//其他进程接收单元界面上的压力值
	    }

	  MPI_Barrier (MPI_COMM_WORLD);

	  PetscPrintf (PETSC_COMM_WORLD, "\nPercentage filled: %.2f%%\n", pf);

	  if (pf > parameter.pf)
	    break;

	}

      if (parameter.adjdt == LOGICAL_TRUE)
	{

	  if (maxCp > parameter.maxCp)
	    dt *= 0.85 * parameter.maxCp / maxCp;

	  if (1.25 * maxCp < parameter.maxCp)
	    dt *= 1.05;

	}

      if (iter >= irestart)
	{

	  // Write restart file  写入重启文件
	  sprintf (file, "%s.%03d.ini", path, processor);

	  fprestart = fopen (file, "wb");

	  if (fprestart != NULL)
	    {

	      WriteRestart (fprestart, curtime);

	    }

	  fclose (fprestart);

	  irestart += parameter.restart;

	}

      if (curtime > wtime + dt)
	{

	  WriteResults (fpresults, LOGICAL_TRUE, LOGICAL_TRUE, curtime);

	  //fflush (fpresults);

	  // Open output file for probes 
	  sprintf (file, "%s.%03d.prb", path, processor);

	  fpprobe = fopen (file, "w");

	  WriteProbeViews (fpprobe, var, curtime);

	  fclose (fpprobe);

	  wtime += wdt;

	}

      if (parameter.steady == LOGICAL_TRUE)
	{

	  // Get residual  获取残差

	  if (parameter.calc[iu] == LOGICAL_TRUE)
	    {

	      VecWAXPY (temp1, -1.0, xu, xu0);
	      VecNorm (temp1, NORM_2, &fres[iu]);

	    }
	  else
	    fres[iu] = 0.0;

	  if (parameter.calc[iv] == LOGICAL_TRUE)
	    {

	      VecWAXPY (temp1, -1.0, xv, xv0);
	      VecNorm (temp1, NORM_2, &fres[iv]);

	    }
	  else
	    fres[iv] = 0.0;

	  if (parameter.calc[iw] == LOGICAL_TRUE)
	    {

	      VecWAXPY (temp1, -1.0, xw, xw0);
	      VecNorm (temp1, NORM_2, &fres[iw]);

	    }
	  else
	    fres[iw] = 0.0;

	  if (parameter.calc[ip] == LOGICAL_TRUE)
	    {

	      VecWAXPY (temp1, -1.0, xp, xp0);
	      VecNorm (temp1, NORM_2, &fres[ip]);

	    }
	  else
	    fres[ip] = 0.0;

	  if (parameter.calc[iT] == LOGICAL_TRUE)
	    {

	      VecWAXPY (temp1, -1.0, xT, xT0);
	      VecNorm (temp1, NORM_2, &fres[iT]);

	    }
	  else
	    fres[iT] = 0.0;

	  if (parameter.calc[is] == LOGICAL_TRUE)
	    {

	      VecWAXPY (temp1, -1.0, xs, xs0);
	      VecNorm (temp1, NORM_2, &fres[is]);

	    }
	  else
	    fres[is] = 0.0;

	  n = 0;

	  for (i = 0; i < nphi; i++)
	    {

	      if (verbose == LOGICAL_TRUE)
		PetscPrintf (PETSC_COMM_WORLD,
			     "\nVariable: %c Iteration: %d Final residual: %+E\n",
			     var[i], fiter[i], fres[i]);

	      if (fres[i] > parameter.ftol[i])
		n++;

	    }

	  if (n == 0)
	    {
	      PetscPrintf (PETSC_COMM_WORLD, "\nSteady state reached.\n");
	      break;
	    }

	  WriteResidual (fpresiduals, iter, fres);

	  fflush (fpresiduals);

	}
      else
	{

	  if (curtime + 0.5 * dt > endtime)
	    break;

	}

    }
  while (dt > 0.0);

  WriteResults (fpresults, LOGICAL_TRUE, LOGICAL_TRUE, curtime);

  // Close output files  关闭输出文件

  fclose (fpresults);

  if (parameter.steady == LOGICAL_TRUE)
    {

      WriteResidual (fpresiduals, iter, fres);

      fflush (fpresiduals);

      PetscFClose (PETSC_COMM_WORLD, fpresiduals);

    }

  // Open output file for probes 
  sprintf (file, "%s.%03d.prb", path, processor);

  fpprobe = fopen (file, "w");

  WriteProbeViews (fpprobe, var, curtime);

  fclose (fpprobe);

  VecGhostRestoreLocalForm (cex, &cexl);
  VecGhostRestoreLocalForm (cey, &ceyl);
  VecGhostRestoreLocalForm (cez, &cezl);

  // Release memory  释放内存
  free (file);

  DeallocateMemory ();

  return LOGICAL_TRUE;

}

void
usage () //操作命令行提示
{

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "Usage:     mpirun -np [n] ../OpenFVM [problem] [options] [n] [petsc-options]\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "Example 1: mpirun -np [n] ../OpenFVM lid f 1\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "Example 2: mpirun -np [n] ../OpenFVM lid f 1 > lid.log\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "Example 3: mpirun -np [n] ../OpenFVM lid d 4\n");

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "Main options:\n");
  PetscPrintf (PETSC_COMM_WORLD, "  r - Reorder the mesh\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "  d - Decompose the mesh into n partitions/processors\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "  f - Start simulation with n partitions/processors\n");
  PetscPrintf (PETSC_COMM_WORLD, "\n");

  PetscPrintf (PETSC_COMM_WORLD, "Additional options:\n");
  PetscPrintf (PETSC_COMM_WORLD, "  v - Verbose mode\n");
  PetscPrintf (PETSC_COMM_WORLD, "  c - Perform checks\n");
  PetscPrintf (PETSC_COMM_WORLD, "\n");

}

int
main (int argc, char **argv) //主函数入口
{

  char *path;
  char *file;

  double start, end;
  double elapsed;

  static char help[] =
    "Three-dimensional unstructured finite-volume implicit flow solver.\n";

  // Start MPI   启动MPI并行过程
  PetscInitialize (&argc, &argv, (char *) 0, help);

  MPI_Comm_size (PETSC_COMM_WORLD, &nbprocessors); //当前通信域包含的进程的个数
  MPI_Comm_rank (PETSC_COMM_WORLD, &processor); // 获取当前的进程标识

  //PetscPrintf (PETSC_COMM_WORLD, "Processor: %d\n", processor);

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "*****************************************\n");
  PetscPrintf (PETSC_COMM_WORLD, "*                                       *\n");
  PetscPrintf (PETSC_COMM_WORLD, "*    OpenFVM-Flow v1.0 - Parallel       *\n");
  PetscPrintf (PETSC_COMM_WORLD, "*                                       *\n");
  PetscPrintf (PETSC_COMM_WORLD, "*****************************************\n");
  PetscPrintf (PETSC_COMM_WORLD, "\n");

  if (argc < 4)
    {
      usage ();
      PetscPrintf (PETSC_COMM_WORLD,
		   "\nError: Wrong number of arguments, expected %d found %d.\n\n",
		   4 - 1, argc - 1);

      // end MPI  
      PetscFinalize ();

      return LOGICAL_ERROR;
    }

  if (strchr (argv[2], 'c') != NULL)
    pchecks = LOGICAL_TRUE;
  else
    pchecks = LOGICAL_FALSE;

  if (strchr (argv[2], 'v') != NULL)
    verbose = LOGICAL_TRUE;
  else
    verbose = LOGICAL_FALSE;

  path = calloc (strlen (argv[1]) + 1, sizeof (char));
  file = calloc (strlen (argv[1]) + 12, sizeof (char));

  // Simulate 数值模拟过程
  if (strchr (argv[2], 'f') != NULL)
    {

      if (nbprocessors < 1)
	{
	  PetscPrintf (PETSC_COMM_WORLD,
		       "\nError: Number of processors cannot be lower than 1.\n\n");

	  // end MPI
	  PetscFinalize ();
	  
	  return LOGICAL_ERROR;
	}

      if (atoi (argv[3]) < 1)
	{
	  PetscPrintf (PETSC_COMM_WORLD,
		       "\nError: Number of partitions cannot be lower than 1.\n\n");

	  // end MPI
	  PetscFinalize ();

	  return LOGICAL_ERROR;
	}

      if (nbprocessors != atoi (argv[3]))
	{
	  PetscPrintf (PETSC_COMM_WORLD,
		       "\nError: Number of processors and partitions are different.\n\n");

	  // end MPI
	  PetscFinalize ();

	  return LOGICAL_ERROR;
	}

      strcpy (path, argv[1]);

      // Read parameter file  读入计算参数文件
      sprintf (file, "%s.par", path);  //将(path).par文件输出到 file 指向的内存中
      ParImportPAR (file);

      // Read mesh file  读入网格文件
      strcpy (file, path);
      sprintf (file, "%s.%03d.msh", path, processor);
      MshImportMSH (file);

      // Read boundary conditions file  读入边界条件文件
      sprintf (file, "%s.bcd", path);
      BcdImportBCD (file);

      // Read material file  读入材料参数文件
      sprintf (file, "%s.mtl", path);
      MtlImportMTL (file);

      MPI_Barrier (MPI_COMM_WORLD);

      // Start clock to measure simulation time 启动时钟，测量仿真计算时间
      start = MPI_Wtime ();

      // Start simulation 开始仿真
      Simulation (path);

      end = MPI_Wtime ();
      elapsed = end - start;
      printf ("\nProcessor: %d, Elapsed time: %.2f seconds.\n", processor,
	      elapsed);

    }

  // Decompose 网格分解
  if (strchr (argv[2], 'd') != NULL)
    {

      strcpy (path, argv[1]);

      // Read mesh file  
      sprintf (file, "%s.msh", path);

      if (MshImportMSH (file) == LOGICAL_ERROR)
	{
	  PetscPrintf (PETSC_COMM_WORLD, "\nError: Mesh file not found!\n");
	  PetscPrintf (PETSC_COMM_WORLD, "%s\n\n", file);
	  return LOGICAL_ERROR;
	}

      // Read boundary conditions file 读取边界条件
      sprintf (file, "%s.bcd", path);
      BcdImportBCD (file);

      // Allocate memory 分配内存
      PetscPrintf (PETSC_COMM_WORLD, "\n");
      PetscPrintf (PETSC_COMM_WORLD, "Allocating memory...\n");

      AllocateMemory ();

      PetscPrintf (PETSC_COMM_WORLD, "Memory allocated.\n");

      SetBoundary ();

      DeallocateMemory ();

      DecomposeMesh (path, atoi (argv[3]));

    }

  // Reorder  重新排序
  if (strchr (argv[2], 'r') != NULL)
    {

      strcpy (path, argv[1]);

      // Read mesh file
      sprintf (file, "%s.msh", path);

      if (MshImportMSH (file) == LOGICAL_ERROR)
	{
	  PetscPrintf (PETSC_COMM_WORLD, "\nError: Mesh file not found!\n");
	  PetscPrintf (PETSC_COMM_WORLD, "%s\n\n", file);
	  return LOGICAL_ERROR;
	}

      ReorderMesh (path);

    }

  free (path);
  free (file);

  // Free memory          
  MshFreeMemory ();

  MPI_Barrier (MPI_COMM_WORLD); //等待进程同步

  PetscPrintf (PETSC_COMM_WORLD, "Done.\n\n");

  // end MPI
  PetscFinalize (); //关闭 MPI调用

  return 0;
}
