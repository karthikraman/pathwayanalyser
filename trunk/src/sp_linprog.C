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
3B
 * dimension. x[0], A[0][i], A[i][0] will all never be used!!!!!!!!!!!!! * 
 *************************************************************************/

#define no_simplex
#define interior

int sp_linprog (double f[], int nz, int iRow[], int jCol[], double dA[], double b[], double x_lb[], double x_ub[], double x[], double *fval, int m, int n) 
{
	LPX *lp;
	
	lp = lpx_create_prob();
	lpx_set_prob_name(lp, "linprog");
	lpx_set_obj_dir(lp, LPX_MIN);
	
	lpx_add_rows(lp, m);
	for (int i=1; i<=m; i++)
		lpx_set_row_bnds(lp, i, LPX_FX, x[i], x[i]);

	lpx_add_cols(lp, n);
	for (int i=1; i<=n; i++)
	{
		lpx_set_col_bnds(lp, i, LPX_DB, x_lb[i], x_ub[i]);
		lpx_set_obj_coef(lp, i, f[i]);
	}
	
	#ifdef DEBUG	
	for (int i=1; i<=nz; i++)
		printf("ia[%d]=%d;\tja[%d]=%d;\tar[%d]=%f\n",i,iRow[i],i,jCol[i],i,dA[i]);
	#endif
	
	lpx_load_matrix(lp, nz, iRow, jCol, dA);
	
	#ifdef simplex
	int exitflag = lpx_simplex(lp);
	*fval = lpx_get_obj_val(lp);
	for (int i=1; i<=n; i++)
		x[i]=lpx_get_col_prim(lp,i);
	#endif
	
	#ifdef interior
	int exitflag = lpx_interior(lp);
	*fval = lpx_ipt_obj_val(lp);
	for (int i=1; i<=n; i++)
		x[i]=lpx_ipt_col_prim(lp,i);
	#endif
	
	lpx_delete_prob(lp);
	free(lp);
	switch(exitflag)
	{
		case LPX_E_OK:
		{
			printf("The LP problem has been successfully solved (to optimality).\n");
			return 0;
		}
		case LPX_E_FAULT:
		{
			printf("The solver can't start the search becuse either the problem has no\n");
			printf("rows and/or no columns, or some row has non-zero objective coefficient.\n");
			return 1;
		}
		case LPX_E_NOFEAS:
		{
			printf("The LP problem has no feasible (primal or dual) solution.\n");
			return 2;
		}
		case LPX_E_NOCONV:
		{
			printf("The search was prematurely terminated due to very slow convergence or divergence.\n");
			return 3;
		}
		case LPX_E_ITLIM:
		{
			printf("The search was prematurely terminated becase the simplex iterations limit has been exceeded.\n");
			return 4;
		}
		case LPX_E_INSTAB:
		{
			printf("The search was prematurely terminated due to numerical instability on solving Newtonian system.\n");
			return 5;
		}
	}

	return exitflag;
}
