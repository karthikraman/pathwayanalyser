/**
 * \file    getStoichiometricMatrix.cpp
 * \brief   Extracts stoichiometric matrix from an SBML file
 * \author  Karthik Raman
 *
 * $Id: getStoichiometricMatrix.cpp$
 * $Source: $ 
 */
/* This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
 * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
 * documentation provided hereunder is on an "as is" basis, and the
 * California Institute of Technology and Japan Science and Technology
 * Corporation have no obligations to provide maintenance, support,
 * updates, enhancements or modifications.  In no event shall the
 * California Institute of Technology or the Japan Science and Technology
 * Corporation be liable to any party for direct, indirect, special,
 * incidental or consequential damages, including lost profits, arising
 * out of the use of this software and its documentation, even if the
 * California Institute of Technology and/or Japan Science and Technology
 * Corporation have been advised of the possibility of such damage.  See
 * the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 *
 * The original code contained here was initially developed by:
 *
 * 	Karthik Raman
 *	Supercomputer Education & Research 
 *		Centre/Bioinformatice Centre
 *	Indian Institute of Science
 *	Bangalore - 560 012
 *	INDIA
 *	
 *	mailto: k.raman AT bioc.uzh.ch
 */

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <map>
#include <vector>

#include "sbml/SBMLReader.h"
#include "sbml/SBMLTypes.h"

using namespace std;

/*   Description:
	This program extracts the stoichiometry matrix for a give SBML file
	and writes it to <filename>.stoich.

	The format of the output file is as follows:

	Line 1: "METABOLITES:"
	Line 2: <Metabolites listed one after the other, tab delimited>
	Line 3: BLANK
	Line 4: "REACTIONS:"
	Line 5: <Reaction IDs listed one after the other, tab delimited>
	Line 6: BLANK
	Line 7: "SPARSE STOICHIOMETRIC MATRIX:"
	Line 8: <Number of non-zero elements>
	Line 9 onwards: 
		Tab separated list of elements, in co-ordinate format
		i	j	A[i][j]

*/
int main(int argc, char *argv[])
{
	typedef map<string, int> StringIntMap ;

	if (argc!=2)
	{
		cerr << "Usage: " << argv[0] << " <filename>" << endl;
		return 1;
	}

	SBMLDocument* d;
	Model* m;

	cout << "Trying to open " << argv[1] << "..." << flush;
	d=readSBML(argv[1]);
  	int errors = d->getNumErrors();

  	if (errors > 0)
  	{
		cout << "Error(s): " << errors       << endl;
		cout << endl;
    		d->printErrors(cout);
		return 1;
  	}
	else
		cout << "successful." << endl; 

	int level=d->getLevel();

	m=d->getModel();	

	int nReactions=m->getNumReactions();
	int nSpecies=m->getNumSpecies();
	int nzStoich=0;
	
	char stoichFileName[1024];
	sprintf(stoichFileName,"%s.stoich",argv[1]);
	ofstream stoichFile(stoichFileName,ios::out);
	
	for (int i=0; i<nReactions; i++)
	{
		Reaction* currentReaction=m->getReaction(i);
		int nReactants=currentReaction->getNumReactants();
		int nProducts=currentReaction->getNumProducts();
		nzStoich+=nReactants+nProducts;
	}

	int nMetab=0;
	
	StringIntMap speciesId;
	StringIntMap geneId;
	
	int 	iRow_S[nzStoich+1];
	double 	data_S[nzStoich+1];
	int 	jCol_S[nzStoich+1];
	
	for (int i=0; i<=nzStoich; i++)
	{
		iRow_S[i]=0;
		jCol_S[i]=0;
		data_S[i]=0;
	}

	int metabolite_num=0;

	StringIntMap::iterator pos;
	int nzStoich_counter=0;
	

	stoichFile << "%METABOLITES: " << endl << "%";
	for (int i=0; i<nReactions; i++)
	{
		Reaction* currentReaction=m->getReaction(i);

		int nReactants=currentReaction->getNumReactants();
		for (int j=0; j<nReactants; j++)
		{
			SpeciesReference *currentSpeciesReference = currentReaction->getReactant(j);
			const string currentSpecies=currentSpeciesReference->getSpecies();
			
			pos=speciesId.find(currentSpecies);	
			if(pos==speciesId.end())	
			{	
				nMetab++;							//new metabolite
				speciesId[currentSpecies]=nMetab; 	//set the metabolite's id in the map
				stoichFile << currentSpecies << "\t";
			}
			
			metabolite_num=speciesId[currentSpecies];
			
			nzStoich_counter++;
			iRow_S[nzStoich_counter]=metabolite_num;
			jCol_S[nzStoich_counter]=i+1;
			data_S[nzStoich_counter]=data_S[nzStoich_counter]-currentSpeciesReference->getStoichiometry();
		}

		int nProducts=currentReaction->getNumProducts();
		for (int j=0; j<nProducts; j++)
		{
			SpeciesReference *currentSpeciesReference = currentReaction->getProduct(j);
			const string currentSpecies=currentSpeciesReference->getSpecies();
			
			pos=speciesId.find(currentSpecies);	
			if(pos==speciesId.end())	
			{	
				nMetab++;							//new metabolite
				speciesId[currentSpecies]=nMetab; 	//set the metabolite's id in the map
				stoichFile << currentSpecies << "\t";
			}
			metabolite_num=speciesId[currentSpecies];
			
			nzStoich_counter++;
			iRow_S[nzStoich_counter]=metabolite_num;
			jCol_S[nzStoich_counter]=i+1;
			data_S[nzStoich_counter]=data_S[nzStoich_counter]+currentSpeciesReference->getStoichiometry();
		}
	}
	cout << endl;
	
	cout << "nSpecies      : " << nSpecies	 << endl;
	cout << "nMetabolites  : " << nMetab << endl;
	cout << "nReactions    : " << nReactions << endl;

	stoichFile << endl;
	stoichFile << "%REACTIONS: " << endl << "%";
	
	for (int i=0; i<nReactions; i++)
	{
		Reaction* currentReaction=m->getReaction(i);
		string rid;
		if (level==1)
			rid = currentReaction->getName();
		else
			rid = currentReaction->getId();

		stoichFile << rid << "\t";
	}
	
	stoichFile << endl;
	
	/* Sorting in row major order...
	   The following code uses bubble sort O(nz^2)
	   need to change it to heapsort O(nz lg nz) */

	int itmp, jtmp;
	double dtmp;

	for (int i=1; i<=nzStoich; i++)
	{
		for (int j=1; j<nzStoich; j++)
		{
			if (iRow_S[j]>iRow_S[j+1] || (iRow_S[j]==iRow_S[j+1] && jCol_S[j]>jCol_S[j+1]))
			{
				itmp=iRow_S[j];
				iRow_S[j]=iRow_S[j+1];
				iRow_S[j+1]=itmp;
				jtmp=jCol_S[j];
				jCol_S[j]=jCol_S[j+1];
				jCol_S[j+1]=jtmp;
				dtmp=data_S[j];
				data_S[j]=data_S[j+1];
				data_S[j+1]=dtmp;
			}
		}
	}

	/*The following code is necessary only if the stoichiometry property has not been used 
	  properly, and instead, we have `n' identical reactants/products in a reaction*/

	/* Begin Patch Code Segment */
	int NZ=nzStoich;
	for (int i=1; i<nzStoich; i++)
	{
		int dup=0;	//duplication of reactants/products, without using stoichiometry aptly
		for (int j=i; iRow_S[j]==iRow_S[j+1] && jCol_S[j]==jCol_S[j+1]; j++)
			dup++;
		for (int j=i; j<=i+dup; j++)
			data_S[j]*=(dup+1);

		i+=dup;
		NZ-=dup;
	}
	/* End Patch Code Segment*/

	cout << "nZ            : " << NZ <<  " (" << 100.0*(double)NZ/((double)(nMetab)*(double)nReactions)<< "%)" << endl;
	stoichFile << "%SPARSE STOICHIOMETRIC MATRIX: " << endl;
	stoichFile << "%" << NZ << endl;
	for (int i=1; i<=nzStoich; i++)
	{
		if (iRow_S[i]==iRow_S[i+1] && jCol_S[i]==jCol_S[i+1])
			continue;
		stoichFile << iRow_S[i] << "\t" << jCol_S[i] << "\t" << data_S[i] << endl;
	}
}
