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

#include <stdlib.h>
#include <time.h>

#ifdef WIN32
#include <windows.h>
#endif /*  */
#include "ttime.h"
double
ttime (void) // 以秒为单位返回系统时间
{

  double sec;

#ifdef WIN32

  sec = GetCurrentTime () * 0.001;

#else

  struct timeval now;

  // Get the system time.
  if (gettimeofday (&now, (struct timezone *) 0))
    {
      return (0);
    }

  //Return time in seconds.  
  sec = (double) now.tv_sec + (((double) now.tv_usec) / 1000000.0);

#endif

  return (sec);

}
