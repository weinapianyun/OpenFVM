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

//PETSc
#include "petscksp.h"

Vec cex, cey, cez;
Vec cexl, ceyl, cezl;

Vec Co;
Vec Col;

Vec uf;

Vec dens, visc, spheat, thcond;
Vec densl, viscl, spheatl, thcondl;

Vec xu0, xv0, xw0, xp0, xT0, xs0;
Vec xu0l, xv0l, xw0l, xp0l, xT0l, xs0l;

Vec xu, xv, xw, xp, xT, xs;
Vec xul, xvl, xwl, xpl, xTl, xsl;

Vec xuf, xvf, xwf, xpf, xTf, xsf;

Mat Am, Ac, Ae, As;

Vec bu, bv, bw, bp, bT, bs;

Vec hu, hv, hw;
Vec hul, hvl, hwl;

Vec ap;
Vec apl;

Vec xpp, xTp;

Vec betaf;

Vec temp1, temp2;
