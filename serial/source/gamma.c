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
#include "msolver.h"

#include "gamma.h"

double
CalculateMaxCourantNumber (double dt, int interface) // 计算最大柯朗数
{
    unsigned int i, j;
    // 使用cpu内部寄存器的变量，加快变量的读写速度; 一般int把变量放在内存中
    register unsigned int face; // 单元界面编号
    register unsigned int element; // 网格单元编号

    double Cpp; // 单元柯朗数
    double Cj; // 中间参数：各单元界面的贡献
    double s; // 单元体积的充填比例(相函数)
    double cs; // 中间参数：？
    double maxCp; // 最大柯朗数

    maxCp = 0.0;

    for (i = 0; i < nbelements; i++)
    {
        element = i;

        if (interface == 1)
        {
            // 保证得到 0<s<1 的网格单元的充填饱和度
            s = LMIN(LMAX(V_GetCmp (&xs, element + 1), 0.0), 1.0);
            cs = (1.0 - s) * (1.0 - s) * s * s * 16.0;
        }
        else
        {
            cs = 1.0;
        }

        Cpp = 0.0;

        for (j = 0; j < elements[element].nbfaces; j++)
        {
            face = elements[element].face[j];
            Cj = LMAX (-V_GetCmp (&uf, face + 1) * faces[face].Aj * dt /
                       elements[element].Vp, 0.0); // 计算单元每个界面的贡献 Co = u*dt*S/V
            Cpp += cs * Cj; // 计算单元的柯朗数
        }

        maxCp = LMAX (maxCp, Cpp);
        V_SetCmp (&Co, element + 1, Cpp);
    }
    return maxCp;
}

void
CorrectFaceGamma () // 修正单元界面上的相函数
{
    int i;
    // 使用cpu内部寄存器的变量，加快变量的读写速度
    register unsigned int face, pair;
    register unsigned int element, neighbor;

    double betaj;

    for (i = 0; i < nbfaces; i++)
    {
        face = i;
        element = faces[face].element;
        pair = faces[face].pair;

        if (pair != -1)
        {
            neighbor = faces[pair].element;

            if (V_GetCmp (&uf, face + 1) > 0.0)
                betaj = V_GetCmp (&betaf, face + 1);
            else
                betaj = 1.0 - V_GetCmp (&betaf, face + 1);

            V_SetCmp (&xsf, face + 1,
                      LMAX (LMIN
                                    ((1.0 - betaj) * V_GetCmp (&xs,
                                                               element + 1) +
                                     betaj * V_GetCmp (&xs, neighbor + 1), 1.0), 0.0));
        }
        else
        {
            if (faces[face].bc == OUTLET)
            {
                // zero gradient
                V_SetCmp (&xsf, face + 1, V_GetCmp (&xs, element + 1));
            }
        }
    }
}

void
PredictBeta () // 预测 beta 值
{
    unsigned int i;

    register unsigned int face, pair;
    register unsigned int element, neighbor;
    register unsigned int donor, acceptor;

    double dot, l1, l2;
    double su, sdn;
    double sjnCBC, sjnUQ;
    double sjn;

    double qj, tetaj;
    double betaj;
    double Cod;
    double ang;

    msh_vector grads;

    for (i = 0; i < nbfaces; i++)
    {
        face = i;
        element = faces[face].element;
        pair = faces[face].pair;
        betaj = 0.0;

        if (parameter.ncicsamcor != 0)
        {
            if (pair != -1)
            {
                neighbor = faces[pair].element;

                if (V_GetCmp (&uf, face + 1) != 0.0)
                {
                    if (V_GetCmp (&uf, face + 1) > 0.0)
                    {
                        acceptor = neighbor;
                        donor = element;

                        grads = Gradient (&xs, &xsf, LOGICAL_FALSE, donor);

                        dot = grads.x * faces[face].d.x + grads.y * faces[face].d.y + grads.z * faces[face].d.z;
                        l1 = GeoMagVector (grads);
                        l2 = GeoMagVector (faces[face].d);
                    }
                    else
                    {
                        acceptor = element;
                        donor = neighbor;

                        grads = Gradient (&xs, &xsf, LOGICAL_FALSE, donor);

                        dot = grads.x * faces[pair].d.x + grads.y * faces[pair].d.y + grads.z * faces[pair].d.z;
                        l1 = GeoMagVector (grads);
                        l2 = GeoMagVector (faces[pair].d);
                    }

                    su = LMIN (LMAX (V_GetCmp (&xs, acceptor + 1) - 2 * dot, 0.0), 1.0);
                    Cod = LMIN(V_GetCmp (&Co, donor + 1), 1.0);

                    if (LABS (V_GetCmp (&xs, acceptor + 1) - su) > SMALL)
                    {
                        sdn = (V_GetCmp (&xs, donor + 1) - su) / (V_GetCmp (&xs, acceptor + 1) - su);

                        if (sdn >= 0.0 && sdn <= 1.0 && LABS (Cod) > SMALL)
                            sjnCBC = LMIN (1.0, sdn / Cod);
                        else
                            sjnCBC = sdn;

                        if (sdn >= 0.0 && sdn <= 1.0)
                            sjnUQ =  LMIN ((8.0 * Cod * sdn + (1.0 - Cod) * (6.0 * sdn + 3.0)) / 8.0, sjnCBC);
                        else
                            sjnUQ = sdn;

                        if (LABS (l1 * l2) > SMALL)
                            ang = LABS (dot / (l1 * l2));
                        else
                            ang = LABS (dot / SMALL);

                        if (ang > 1.0)
                            ang = 1.0;

                        tetaj = acos (ang);

                        qj = LMIN (parameter.kq * 0.5 * (cos (2 * tetaj) + 1.0), 1.0);
                        sjn = qj * sjnCBC + (1 - qj) * sjnUQ;

                        if (LABS (1.0 - sdn) > SMALL)
                        {
                            betaj = LMIN (LMAX ((sjn - sdn) / (1.0 - sdn), 0.0), 1.0);
                        }
                    }
                }
            }
        }
        V_SetCmp (&betaf, face + 1, betaj);
    }

}

void
CorrectBeta (double dt) // 修正 beta 值
{
    unsigned int i;

    unsigned int face, pair;
    register unsigned int element;
    unsigned int neighbor;
    unsigned int donor, acceptor;

    double Cj;
    double cbetaj, betaj;
    double ds, Ep, Em;

    for (i = 0; i < nbfaces; i++)
    {
        face = i;
        element = faces[face].element;
        pair = faces[face].pair;
        cbetaj = 0.0;

        betaj = V_GetCmp (&betaf, face + 1);

        if (betaj < 1E-2) continue;

        if (pair != -1)
        {
            neighbor = faces[pair].element;

            if (V_GetCmp (&uf, face + 1) != 0.0)
            {
                if (V_GetCmp (&uf, face + 1) > 0.0)
                {
                    acceptor = neighbor;
                    donor = element;
                }
                else
                {
                    acceptor = element;
                    donor = neighbor;
                }

                Cj =LMIN (LMAX(-V_GetCmp (&uf, face + 1) *
                               faces[face].Aj * dt /elements[element].Vp, 0.0), 1.0);
                ds =0.5 * (V_GetCmp (&xs0, acceptor + 1) +V_GetCmp (&xs, acceptor + 1)) -
                    0.5 * (V_GetCmp (&xs0, donor + 1) +V_GetCmp (&xs, donor + 1));
                if (V_GetCmp (&xs, donor + 1) < 0.0)
                {
                    Em = LMAX (-V_GetCmp (&xs, donor + 1), 0.0);
                    // Donor value < 0.0 Ex: sd = -0.1 -> Em = +0.1
                    if (Em > SMALL && Cj > SMALL)
                    {
                        if (ds > Em)
                        {
                            cbetaj =Em * (2 + Cj -2 * Cj * betaj) / (2 * Cj * (ds - Em));
                            cbetaj = LMIN (cbetaj, betaj);
                        }
                    }
                }
                if (V_GetCmp (&xs, donor + 1) > 1.0)
                {
                    Ep = LMAX (V_GetCmp (&xs, donor + 1) - 1.0, 0.0);
                    // Donor value > 1.0 Ex: sd = 1.1 -> Ep = +0.1
                    if (Ep > SMALL && Cj > SMALL)
                    {
                        if (ds < -Ep)
                        {
                            cbetaj =Ep * (2 + Cj -2 * Cj * betaj) / (2 * Cj * (-ds - Ep));
                            cbetaj = LMIN (cbetaj, betaj);
                        }
                    }
                }
            }
        }
        betaj -= cbetaj;
        betaj = LMAX (betaj, 0.0);
        V_SetCmp (&betaf, face + 1, betaj);
    }
}

void
BoundScalar (Vector * x, double min, double max)
{
    unsigned int i;
    register unsigned int element;

    for (i = 0; i < nbelements; i++)
    {
        element = i;

        V_SetCmp (x, element + 1,LMAX (LMIN (V_GetCmp (x, element + 1), max), min));
    }
}

void
SmoothScalar (Vector * xm, Vector * x, int n)
{
    unsigned int i, j, k;
    unsigned int face, pair;

    register unsigned int element;
    unsigned int neighbor;

    double sj;

    //double dNf, dPf;
    double lambda;
    double sum1, sum2;
    double *sa;
    double *sm;

    sa = calloc (nbelements, sizeof (double));
    sm = calloc (nbelements, sizeof (double));

    for (i = 0; i < nbelements; i++)
    {
        element = i;

        sa[i] = V_GetCmp (x, element + 1);
        sm[i] = 0.0;
    }

    for (k = 0; k < n; k++)
    {
        for (i = 0; i < nbelements; i++)
        {
            element = i;

            sum1 = 0.0;
            sum2 = 0.0;
            for (j = 0; j < elements[element].nbfaces; j++)
            {
                face = elements[element].face[j];
                pair = faces[face].pair;

                if (pair != -1)
                {
                    neighbor = faces[pair].element;
                    /*
                    dNf =
                      GeoMagVector (GeoSubVectorVector
                            (elements[neighbor].celement,
                             faces[pair].cface));
                    dPf =
                      GeoMagVector (GeoSubVectorVector
                            (elements[element].celement,
                             faces[face].cface));

                    lambda = dPf / (dPf + dNf);
                    */
                    lambda = 0.5;
                    sj = sa[neighbor] * lambda + sa[element] * (1.0 - lambda);
                }
                else
                {
                    sj = sa[element];
                }

                sum1 += sj * faces[face].Aj;
                sum2 += faces[face].Aj;
            }
            sm[i] = sum1 / sum2;
        }

        for (i = 0; i < nbelements; i++)
        {
            element = i;
            sa[i] = sm[i];
        }
    }

    for (i = 0; i < nbelements; i++)
    {
        element = i;

        V_SetCmp (xm, element + 1, sm[i]);
    }
    free (sa);
    free (sm);
}

void
CalculateMassFraction () // 计算流场充填流体的质量
{
    unsigned int i;
    register unsigned int element;

    double f[2];
    double vol[2];

    vol[0] = 0.0;
    vol[1] = 0.0;

    for (i = 0; i < nbelements; i++)
    {
        element = i;

        f[0] = (1.0 - V_GetCmp (&xs, element + 1));
        f[1] = V_GetCmp (&xs, element + 1);

        vol[0] += f[0] * V_GetCmp (&dens, element + 1) * elements[element].Vp;
        vol[1] += f[1] * V_GetCmp (&dens, element + 1) * elements[element].Vp;
    }
    printf ("\nMass of fluid %d: %+E kg\n", 0, vol[0]);
    printf ("\nMass of fluid %d: %+E kg\n", 1, vol[1]);
}

void
BuildVolumeOfFluidMatrix (double dt) // 创建 VOF 矩阵
{
    unsigned int i, j, n;
    unsigned int face, pair;
    register unsigned int element;
    unsigned int neighbor;

    double aip;
    double ain[MAXFACES];
    unsigned int ani[MAXFACES];
    double bip;

    double betaj;

    // Equation: ds/dt + div(s*U) = 0

    for (i = 0; i < nbelements; i++)
    {
        element = i;
        aip = 0;
        bip = 0;
        n = 0;

        for (j = 0; j < elements[element].nbfaces; j++)
        {
            face = elements[element].face[j];
            pair = faces[face].pair;

            if (parameter.scheme[is] == UDS)
            {
                // UDS
                if (V_GetCmp (&uf, face + 1) > 0.0)
                    betaj = 0.0;
                else
                    betaj = 1.0;
            }
            else
            {
                // CDS
                if (V_GetCmp (&uf, face + 1) > 0.0)
                    betaj = V_GetCmp (&betaf, face + 1);
                else
                    betaj = 1.0 - V_GetCmp (&betaf, face + 1);
            }

            if (pair != -1)
            {
                neighbor = faces[pair].element;

                aip += 0.5 * (1.0 - betaj) * V_GetCmp (&uf, face + 1) * faces[face].Aj;
                ain[n] = 0.5 * betaj * V_GetCmp (&uf, face + 1) * faces[face].Aj;
                ani[n] = neighbor;
                n++;

                bip += -0.5 * (1.0 - betaj) * V_GetCmp (&uf, face + 1) *
                        faces[face].Aj * V_GetCmp (&xs, element + 1);
                bip += -0.5 * betaj * V_GetCmp (&uf, face + 1) *
                        faces[face].Aj * V_GetCmp (&xs, neighbor + 1);
            }
            else
            {
                bip += -1.0 * V_GetCmp (&uf, face + 1) * faces[face].Aj * V_GetCmp (&xsf, face + 1);
            }
        }
        if (dt > 0.0)
        {
            aip += elements[element].Vp / dt;
            bip += elements[element].Vp / dt * V_GetCmp (&xs0, element + 1);
        }

        if (aip == 0.0 || aip != aip)
        {
            printf ("\nError: Problem setting up volume-of-fluid matrix\n");
            exit (LOGICAL_ERROR);
        }

        Q_SetLen (&As, element + 1, n + 1);
        Q_SetEntry (&As, element + 1, 0, element + 1, aip);

        for (j = 0; j < n; j++)
        {
            Q_SetEntry (&As, element + 1, j + 1, ani[j] + 1, ain[j]);
        }
        V_SetCmp (&bs, element + 1, bip);
    }
}

void
SolveVolumeOfFluidExplicit (double dt) // 显式求解流场单元体积比
{
    unsigned int i, j, n;

    // 直接调用CPU寄存器中变量
    unsigned int face, pair;
    register unsigned int element;
    unsigned int neighbor;

    double aip;
    double ain[MAXFACES];
    unsigned int ani[MAXFACES];
    double bip;
    double betaj;
    double sums;

    // Equation: ds/dt + div(s*U) = 0

    for (i = 0; i < nbelements; i++)
    {
        element = i;
        aip = 0;
        bip = 0;
        n = 0;
        for (j = 0; j < elements[element].nbfaces; j++)
        {
            face = elements[element].face[j];
            pair = faces[face].pair;
            if (parameter.scheme[is] == UDS)
            {
                // UDS
                if (V_GetCmp (&uf, face + 1) > 0.0)
                    betaj = 0.0;
                else
                    betaj = 1.0;
            }
            else
            {
                // CDS
                if (V_GetCmp (&uf, face + 1) > 0.0)
                    betaj = V_GetCmp (&betaf, face + 1);
                else
                    betaj = 1.0 - V_GetCmp (&betaf, face + 1);
            }

            if (pair != -1)
            {
                neighbor = faces[pair].element;

                aip += 0.5 * (1.0 - betaj) * V_GetCmp (&uf, face + 1) * faces[face].Aj;
                ain[n] = 0.5 * betaj * V_GetCmp (&uf, face + 1) * faces[face].Aj;
                ani[n] = neighbor;
                n++;

                bip += -0.5 * (1.0 - betaj) * V_GetCmp (&uf, face + 1) *
                       faces[face].Aj * V_GetCmp (&xs, element + 1);
                bip += -0.5 * betaj * V_GetCmp (&uf, face + 1) *
                       faces[face].Aj * V_GetCmp (&xs, neighbor + 1);

            }
            else
            {
                bip += -1.0 * V_GetCmp (&uf, face + 1) * faces[face].Aj *
                       V_GetCmp (&xsf, face + 1);
            }
        }

        if (dt > 0.0)
        {
            aip += elements[element].Vp / dt;
            bip += elements[element].Vp / dt * V_GetCmp (&xs0, element + 1);
        }

        sums = 0.0;
        for (j = 0; j < n; j++)
        {
            sums += ain[j] * V_GetCmp (&xs0, ani[j] + 1);
        }

        V_SetCmp (&xs, element + 1, (bip - sums) / aip);
    }
}

void
CalculateGamma (char *var, int *fiter, double dt, double *maxCp, int verbose,
                int pchecks) // 计算流场的相函数
{
    unsigned int i, j;

    double mres; // 残差
    int miter; // 迭代
    double mtime; // 时间

    if (parameter.calc[is] == LOGICAL_FALSE)
        return;

    if (parameter.vofastemp == LOGICAL_TRUE)
    {
        // The value of T is assigned to s
        Asgn_VV (&xs, &xT);
        return;
    }

    *maxCp = CalculateMaxCourantNumber (dt, 1);
    parameter.ncicsamsteps = LMIN ((int) (*maxCp * 2) + 1, 100);

    if (verbose == LOGICAL_TRUE)
        printf ("\nMaximum Courant number at interface: %.3f\n", *maxCp);

    V_Constr (&betaf, "CICSAM interpolation factor", nbfaces, Normal, True);

    for (i = 0; i < parameter.ncicsamsteps; i++)
    {
        *maxCp = CalculateMaxCourantNumber (dt / parameter.ncicsamsteps, 0);
        // Store previous time step values
        Asgn_VV (&xs0, &xs);

        fiter[is]++;

        // Predict beta - CICSAM
        PredictBeta (&betaf, &xsf, &xs, &Co, &uf);

        for (j = 0; j <= parameter.ncicsamcor; j++)
        {
            if (parameter.timemethod[is] == IMPLICITEULER)
            {
                Q_Constr (&As, "Indicator function matrix", nbelements, False,
                          Rowws, Normal, True);
                V_Constr (&bs, "Gamma source", nbelements, Normal, True);
            }

            if (parameter.timemethod[is] == EXPLICITEULER)
            {
                // Matrix free VOF
                SolveVolumeOfFluidExplicit (dt / parameter.ncicsamsteps);
            }

            if (parameter.timemethod[is] == IMPLICITEULER)
            {
                // Build VOF matrix
                BuildVolumeOfFluidMatrix (dt / parameter.ncicsamsteps);

                if (pchecks == LOGICAL_TRUE)
                {
                    if (!CheckIfDiagonalMatrix (&As))
                    {
                        printf
                                ("\nWarning: Volume-of-fluid matrix is not diagonal dominant\n");
                        WriteMatrix (&As, LOGICAL_FALSE);
                        WriteVector (&bs);
                        //exit (LOGICAL_ERROR);
                    }
                }

                // Set matrix solution accuracy
                SetRTCAccuracy (parameter.mtol[is]);

                // Solve matrix to get indicator function s
                SolveMatrix (&As, &xs, &bs, &miter, &mres, &mtime,
                             parameter.msolver[is], parameter.mprecond[is],
                             parameter.miter[is]);

                if (verbose == LOGICAL_TRUE)
                    printf
                            ("\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
                             var[is], miter, mres, mtime);

                if ((mres > parameter.mtol[is] && miter == parameter.miter[is])
                    || LASResult () != LASOK)
                {
                    printf ("\nError: Problem solving matrix %c\n", var[is]);
                    exit (LOGICAL_ERROR);
                }
            }

            // Correct beta
            CorrectBeta (dt / parameter.ncicsamsteps);

            // Correct face and boundary
            CorrectFaceGamma ();

            if (parameter.timemethod[is] == IMPLICITEULER)
            {
                Q_Destr (&As);
                V_Destr (&bs);
            }
        }

        if (parameter.fill == LOGICAL_TRUE)
        {
            // Bound volume fraction
            BoundScalar (&xs, 0.0, 1.0);
        }
        // Smooth volume fraction
        //SmoothScalar (&xsm, &xs, 2);
    }

    if (pchecks == LOGICAL_TRUE)
    {
        // Calculate mass fractions
        CalculateMassFraction ();
    }

    V_Destr (&betaf);
}
