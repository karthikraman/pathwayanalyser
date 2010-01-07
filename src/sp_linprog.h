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

#include <cstdio>
#include <cstdlib>
extern "C" {
#include "glpk.h"
}
using namespace std;

/*************************************************************************
 * MINIMISATION ROUTINE (not *MAX* !!!!!!!)                              *
 *                                                                       *
 * This routine is aimed to solve linear programming problems of the     *   
 * $\min \vec{f}^T \vec{x} $ s.t. $\mat{A} \vec{x}=\vec{b}$              *
 * and $x_{lb_i} \le x \le x_{ub_i}, i = 1..m$                           *
 * $\mat{A}$ is an $m \times n$ matrix                                   *
 * $\vec{x}$ is $n \times 1$                                             *
 * $\vec{b}$ is obviously $m \times 1$                                   *
 * $\vec{f}$ is $n \times 1$                                             *
 *                                                                       *
 * This routine attempts to partially mimic the syntax of the linprog    *
 * routine the MATLAB Optimization Toolbox                               *
 *                                                                       *
 *                                                                       *
 * Important WARNING: All arrays passed will be indexed from one onwards *
 * to size. Therefore, declare all arrays with one element extra in each *
 * dimension. x[0], A[0][i], A[i][0] will all never be used!!!!!!!!!!!!! * 
 *************************************************************************/

int sp_linprog (double f[], int nz, int iRow[], int jCol[], double dA[], double b[], double x_lb[], double x_ub[], double x[], double *fval, int m, int n); 
