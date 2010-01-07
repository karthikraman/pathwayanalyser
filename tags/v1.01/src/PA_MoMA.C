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

#include "PA_MoMA.h"
const int maxFilenameLength=1024;

int main(int argc, char *argv[])
{
	int ch;
	char *modelName=NULL;
	
	extern char *optarg;
	extern int optopt;
	
	const char *options="f:h";

	int nErrors=0;
	

	/* Parse Command-Line Arguments */
	
	if ((ch = getopt(argc, argv, options)) == -1)
	{
		printf("Usage: %s -f <modelName (without extension)>\n", argv[0]);
		return 1;
	}
		
	do
	{
		switch(ch)
		{
			case 'f':
				modelName = optarg;
				printf("Model Name is %s\n", modelName);
				break;

			case 'h':
				printf("Usage: %s -f <modelName (without extension)> [-t]\n", argv[0]);
				printf("\t-f\tName of the SBML file containing the model, without extension\n");
				printf("\t-h\tDisplays this help\n");
				return 0;

			case '?':
				return 1;
				break;
			
			case ':':
				printf("unknown arg %c\n", optopt);
				return 1;
				break;
			
			default:
				printf("%s requires at least one argument!\n", argv[0]);
				return 1;
				break;
		}
	} while ((ch = getopt(argc, argv, options)) != -1);
	
	if (modelName==NULL)
	{
		fprintf(stderr, "Error! Model name not specified!\n");
		return 1;
	}	

	char objFnFileName[maxFilenameLength];
	sprintf(objFnFileName,"%s.objfn",modelName);
	ifstream objFnFile(objFnFileName,ios::in);

	char lbFileName[maxFilenameLength];
	sprintf(lbFileName,"%s.lb",modelName);
	ifstream lbFile(lbFileName, ios::in);

	char ubFileName[maxFilenameLength];
	sprintf(ubFileName,"%s.ub",modelName);
	ifstream ubFile(ubFileName, ios::in);


	if (!objFnFile)
	{
		fprintf(stderr, "Could not open objective function file %s!!!\n",objFnFileName);
		return 1;
	}
	if (!lbFile)
	{
		fprintf(stderr, "Could not open lower bounds file %s!!!\n",lbFileName);
		return 1;
	}
	if (!ubFile)
	{
		fprintf(stderr, "Could not open upper bounds file %s!!!\n",ubFileName);
		return 1;
	}
	/* Argument parsing complete */

	SBMLDocument* PA_doc;
	Model* PA_model;

	char modelFileName[maxFilenameLength];
	sprintf(modelFileName,"%s.sbml",modelName);

	cout << "Trying to open " << modelFileName << "..." << flush;
	PA_doc=readSBML(modelFileName);
  	int errors = PA_doc->getNumErrors();

  	if (errors > 0)
  	{
		cerr << "Error(s): " << errors       << endl;
		cerr << endl;
    		PA_doc->printErrors(cerr);
		return 1;
  	}
	else
		cout << "successful." << endl;

	PA_model=PA_doc->getModel();	
	
	int nReactions=PA_model->getNumReactions();
	int nSpecies=PA_model->getNumSpecies();
	string modelDesc=PA_model->getName();
	int nzStoich=0;
	
	int nModifiers=0;
	int speciesCounter=0;
	int nGenes=0;

	for (int i=0; i<nReactions; i++)
	{
		Reaction* currentReaction=PA_model->getReaction(i);
		int nReactants=currentReaction->getNumReactants();
		int nProducts=currentReaction->getNumProducts();
		nModifiers+=currentReaction->getNumModifiers();	//this is approximate, but is certainly an upper bound
		nzStoich+=nReactants+nProducts;
	}

	typedef map<string, int> StringIntMap;	
	StringIntMap speciesId;
	StringIntMap geneId;
	map<int, string> geneName;
	
	int 	iRow_S[nzStoich+1];
	double 	data_S[nzStoich+1];
	int 	jCol_S[nzStoich+1];
	
	for (int i=0; i<=nzStoich; i++)
	{
		iRow_S[i]=0;
		jCol_S[i]=0;
		data_S[i]=0;
	}

	int 	iRow_geneRn[nModifiers+1];
	int 	jCol_geneRn[nModifiers+1];

	for (int i=0; i<=nModifiers; i++)
	{
		iRow_geneRn[i]=0;
		jCol_geneRn[i]=0;
	}
	
	bool ext[nReactions+1];
	double rev[nReactions+1];
	
	for (int j=0; j<=nReactions; j++)
		ext[j]=true;

	int metabolite_num=0;
	int gene_num=0;
	int nzStoich_counter=0;
	int nzgeneRn_counter=0;
	StringIntMap::iterator pos;
	
	for (int i=0; i<nReactions; i++)
	{
		Reaction* currentReaction=PA_model->getReaction(i);
		int nReactants=Reaction_getNumReactants(currentReaction);
		rev[i+1]=Reaction_getReversible(currentReaction);

		for (int j=0; j<nReactants; j++)
		{
			SpeciesReference *currentSpeciesReference = currentReaction->getReactant(j);
			const string currentSpecies=currentSpeciesReference->getSpecies();
			pos=speciesId.find(currentSpecies);
			if (pos==speciesId.end())	
			{	
				speciesCounter++;				//new metabolite
				speciesId[currentSpecies]=speciesCounter; 	//set the metabolite's id in the map
			}	
			metabolite_num=speciesId[currentSpecies];
			nzStoich_counter++;
			iRow_S[nzStoich_counter]=metabolite_num;
			jCol_S[nzStoich_counter]=i+1;
			data_S[nzStoich_counter]=data_S[nzStoich_counter]-currentSpeciesReference->getStoichiometry();
		}
		int nProducts=Reaction_getNumProducts(currentReaction);
		for (int j=0; j<nProducts; j++)
		{
			ext[i+1]=false;
			/* 	AN important assumption that is made is that all external reactions are 
				indicated as material exiting the system, rather than coming into the
				system, i.e., all reactions that have no products but only reactants
				are external reactions.
			*/
			SpeciesReference *currentSpeciesReference = currentReaction->getProduct(j);
			const string currentSpecies=currentSpeciesReference->getSpecies();
			pos=speciesId.find(currentSpecies);
			if (pos==speciesId.end())
			{	
				speciesCounter++;				//new metabolite
				speciesId[currentSpecies]=speciesCounter; 	//set the metabolite's id in the map
			}	
			metabolite_num=speciesId[currentSpecies];
			nzStoich_counter++;
			iRow_S[nzStoich_counter]=metabolite_num;
			jCol_S[nzStoich_counter]=i+1;
			data_S[nzStoich_counter]=data_S[nzStoich_counter]+currentSpeciesReference->getStoichiometry();
		}
		
		int nModifiers=currentReaction->getNumModifiers();
		for (int j=0; j<nModifiers; j++)
		{
			ModifierSpeciesReference *currentModifierSpeciesReference = currentReaction->getModifier(j);
			string currentModifierSpecies=currentModifierSpeciesReference->getSpecies();
			pos=geneId.find(currentModifierSpecies);
			if (pos==geneId.end())
			{
				nGenes++;
				geneId[currentModifierSpecies]=nGenes;
				geneName[nGenes]=PA_model->getSpecies(currentModifierSpecies)->getName();
			}
			gene_num=geneId[currentModifierSpecies];
			nzgeneRn_counter++;
			iRow_geneRn[nzgeneRn_counter]=gene_num;
			jCol_geneRn[nzgeneRn_counter]=i+1;
		}	
	}

	cout << "nSpecies      : " << speciesCounter << endl;
	cout << "nReactions    : " << nReactions << endl;
	
	/* Begin FBA here! */

	int m=nSpecies;
	int n=nReactions;
	double f[n+1];

	double b[m+1], x[n+1];
	double x_lb[n+1], x_ub[n+1], LB[n+1], UB[n+1];
	double fval;

	for (int i=0; i<=m; i++)
		b[i]=0;

	f[0]=0;
	x_lb[0]=0;
	x_ub[0]=0;
	UB[n+1]=0;
	LB[0]=0;
	x[0]=0;
	for (int i=1; i<=n; i++)
	{
		f[i]=0;
		lbFile >> LB[i];
		x_lb[i]=LB[i];
		ubFile >> UB[i];
		x_ub[i]=UB[i];
		x[i]=0;
	}
	
	//Objective function file does exist
	int k=0;
	double coeff;
	while (!objFnFile.eof())
	{
		objFnFile >> k >> coeff;
		f[k]=coeff;
	}
	objFnFile.close();
	
	int exitflag;
	
	/*
	Comprehensive Gene deletion is to be done here:
	*/
	
	char reportFileName[maxFilenameLength];
	sprintf(reportFileName,"%s.MoMA.delReport",modelName);
	ofstream reportFile(reportFileName, ios::out);

	char logFileName[maxFilenameLength];
	sprintf(logFileName,"%s.MoMA.log",modelName);
	ofstream logFile(logFileName, ios::out);
	
	char fluxFileName[maxFilenameLength];
	sprintf(fluxFileName,"%s.MoMA.flux",modelName);
	ofstream fluxFile(fluxFileName, ios::out);

	reportFile << "Gene Deletion Report for " << modelFileName << " (" << modelDesc << ")" << endl;
	reportFile << endl;
	reportFile << "Mutant        \t" << "Function Value\tEssentiality" << endl;
	reportFile << endl;
	

	exitflag=sp_linprog(f, nzStoich, iRow_S, jCol_S, data_S, b, x_lb, x_ub, x, &fval, m, n);
	if (exitflag!=0)
	{
		nErrors++;
		logFile << "Error while performing FBA for wild type! " << geneName[k] << "; linprog exited with code " << exitflag << endl;
	 	reportFile << setiosflags(ios::left) << setw(15) << "Wild" << "\t" << setiosflags(ios::fixed) << setprecision(6) << fval << " \t - \t ERR!" << endl;
	}
	else
	{
	 	reportFile << setiosflags(ios::left) << setw(15) << "Wild" << "\t" << setiosflags(ios::fixed) << setprecision(6) << fval << " \t -" << endl;
	}
	
	fluxFile << "Wild:" << endl;
	fluxFile.precision(6);
	for (int j=1; j<=nReactions; j++)
	{
		fluxFile << x[j] << endl;
	}

	double wild_fval=fval;

	/*	No more of GLPK. SO make sure all arrays and matrices are 
		indexed 0->n-1!!!!!
	*/
	
	const int nx = nReactions;
		
	char ix_lb[nx];	
	char ix_ub[nx];

	for (int i=0; i<nx; i++)
	{
		ix_lb[i]=1;
		ix_ub[i]=1;
	}	

	const int nnzQ = nx;
	int	irowQ[nx];
	int	jcolQ[nx];
	double dQ[nx];
	double c[nx];

	for (int i=0; i<nx; i++)
	{
		irowQ[i]=i;
		jcolQ[i]=i;
		dQ[i]=1.0;	//Q is identity matrix!
		c[i]=-x[i+1]; 	//c is now wild-type function value
	}

	int my = m;
	int nnzA = nzStoich;
	int irowA[nzStoich];
	int jcolA[nzStoich];
	double dA[nzStoich];

	for (int i=0; i<nnzA; i++)
	{
		irowA[i]=iRow_S[i+1]-1;	//zero indexing vs one indexing
		jcolA[i]=jCol_S[i+1]-1;	//zero indexing vs one indexing
		dA[i]=data_S[i+1];
	}
	
	//doubleLexSort(irowA,nnzA,jcolA,dA); -- Where is it in OOQP!!!

	/* Sorting in row major order...
	   The following code uses bubble sort O(nz^2)
	   need to change it to heapsort O(nz lg nz) */

	int itmp, jtmp;
	double dtmp;

	for (int i=0; i<nnzA; i++)
	{
		for (int j=0; j<nnzA; j++)
		{
			if (irowA[j]>irowA[j+1] || (irowA[j]==irowA[j+1] && jcolA[j]>jcolA[j+1]))
			{
				itmp=irowA[j];
				irowA[j]=irowA[j+1];
				irowA[j+1]=itmp;
				jtmp=jcolA[j];
				jcolA[j]=jcolA[j+1];
				jcolA[j+1]=jtmp;
				dtmp=dA[j];
				dA[j]=dA[j+1];
				dA[j+1]=dtmp;
			}
		}
	}


	const int mz = 0;
	double *clow = 0;
	char  *iclow = 0;

	double *cupp = 0;
	char  *icupp = 0;

	const int nnzC = 0;
	int   * irowC  = 0;
	int   * jcolC  = 0;
	double  * dC   = 0;

	cout << endl;

	for (int k=1; k<=nGenes; k++)
	{

		cout << endl;
		cout << "Performing MoMA for " << geneName[k] << endl;

		for (int i=0; i<=n; i++)
		{
			x[i]=0;		//reset x value
			x_ub[i]=UB[i+1];//also shift
			x_lb[i]=LB[i+1];
		}

		for (int j=1; j<=nzgeneRn_counter; j++)
		{
			if(iRow_geneRn[j]==k)
			{
				x_ub[jCol_geneRn[j]-1]=0;
				x_lb[jCol_geneRn[j]-1]=0;
			}
		}

		QpGenSparseMa27 *qp = new QpGenSparseMa27(nx, my, mz, nnzQ, nnzA, nnzC);
		QpGenData *prob = (QpGenData*) qp->copyDataFromSparseTriple(
				c,      irowQ,  nnzQ,   jcolQ,  dQ,
				x_lb,   ix_lb,  x_ub,   ix_ub,
				irowA,  nnzA,   jcolA,  dA,  	b,
				irowC,  nnzC,   jcolC,  dC,
				clow,   iclow,  cupp,   icupp );
		QpGenVars *vars = (QpGenVars*) qp->makeVariables(prob);
		QpGenResiduals *resid = (QpGenResiduals*) qp->makeResiduals(prob);
		GondzioSolver  *s = new GondzioSolver(qp, prob);

		bool quiet=1; 
		if(!quiet)
			s->monitorSelf();

		int ierr = s->solve(prob,vars, resid);
		  
		if(ierr==0)
		{
			fluxFile << geneName[k] << ":" << endl;
			fluxFile.precision(6);
			vars->x->writefToStream(fluxFile, "%{value}");
		}
		else
		{
			fluxFile << geneName[k] << ": ERROR!" << endl;
			fluxFile.precision(6);
			vars->x->writefToStream(fluxFile, "%{value}");
			cout << "Could not solve the problem.\n";
		}

		exitflag=ierr;

		if (exitflag!=0)
		{
			nErrors++;
			logFile << "Error while performing deletion of Gene " << geneName[k] << "; linprog exited with code " << exitflag << endl;
		}

		//Calculate fval here
	
		fval=0;
		double q[n];
		vars->x->copyIntoArray(q);
		for (int i=0; i<n; i++)
		{
			fval+=f[i+1]*q[i];
//			cout << i << ": f["<<i<<"] = "<<f[i+1]<<"; x["<<i<<"] = " << q[i] << "; " << x_lb[i] << ", "<<x_ub[i] << endl;
		}
	 	reportFile << setiosflags(ios::left) << setw(15) << geneName[k] << "\t" << setiosflags(ios::fixed) << setprecision(6) << fval << " \t";

		//We use the 5% viability criterion here

		if (fabs(fval)>=0.05*fabs(wild_fval))
			reportFile << "NON-ESSENTIAL" << endl;
		else
			reportFile << "ESSENTIAL" << endl;

		cout << "done." << endl;
	}
	fluxFile.close();
	logFile.close();
	reportFile.close();
	
	cout << endl << "Deletion report written to " << reportFileName << endl;

	if (nErrors)
	{
		cout << "There were " << nErrors << " errors during the execution of the program." << endl;
		cout << "Please see log file for details" << endl;
	}
	else
	{
		cout << "The in silico gene deletion studies were successfully completed for all the genes." << endl;
	}
}
