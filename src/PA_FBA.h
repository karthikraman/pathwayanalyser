/***************************************************************************
 *   Copyright (C) 2006 by Karthik Raman   *
 *   karthik@rishi.serc.iisc.ernet.in   *
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

#include <iostream>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <unistd.h>

#include <map>
#include <vector>

#include "sp_linprog.h"

#include "sbml/SBMLDocument.h"
#include "sbml/SBMLWriter.h"
#include "sbml/SBMLTypes.h"

using namespace std;
