#!/bin/bash
echo "This script will help you set up PathwayAnalyser on your system..."
echo
echo -n "Directory where PA should be installed [$HOME/PathwayAnalyser-v1]? "
read installdir
[[ $installdir == "" ]] && installdir=$HOME/PathwayAnalyser-v1
echo -n "SBML Libraries Directory (e.g. /usr/lib/sbml)? "
read sbml_lib_dir
echo -n "SBML Includes Directory (e.g. /usr/include/sbml)? "
read sbml_incl_dir
echo -n "Xerces-C Libraries Directory (e.g. /usr/local/xercesc/lib)? "
read xercesc_lib_dir
echo -n "Xerces-C Includes Directory (e.g. /usr/local/xercesc/include)? "
read xercesc_incl_dir
echo -n "GLPK Libraries Directory (e.g. /usr/lib/glpk)? "
read glpk_lib_dir
echo -n "GLPK Includes Directory (e.g. /usr/include/glpk)? "
read glpk_incl_dir
echo -n "OOQP Libraries Directory (e.g. /usr/local/OOQP/lib)? "
read ooqp_lib_dir
echo -n "OOQP Includes Directory (e.g. /usr/local/OOQP/include)? "
read ooqp_incl_dir
echo -n "OOQP MA27 libMA27.a location (e.g. /usr/local/OOQP/extras/MA27/libMA27.a)? "
read libMA27

echo "Installing in '$installdir'..."
mkdir -p $installdir #Actually this should be the last step.. after the make
cd src
g++ sp_linprog.C -c
g++ PA_FBA.C sp_linprog.o -O2 -I $sbml_inc_dir -I $glpk_inc_dir -I $xercesc_inc_dir -L $xercesc_lib_dir -L $sbml_lib_dir -L $glpk_lib_dir -lsbml -lglpk -lxerces-c -o $install_dir/PA_FBA
g++ PA_MoMA.C sp_linprog.o -O2 -I $sbml_inc_dir -I $glpk_inc_dir -I $xercesc_inc_dir -I $ooqp_inc_dir -L $xercesc_lib_dir -L $sbml_lib_dir -L $glpk_lib_dir -L $ooqp_lib_dir -lsbml -lxerces-c -lglpk -looqpgensparse -looqpsparse  -looqpgondzio -looqpbase -lblas $libMA27 -lg2c -lgfortran -o $install_dir/PA_MoMA

