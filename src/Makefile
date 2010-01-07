CC = g++
CFLAGS = -c 
OFLAGS = -o
OPTFLAGS = -O3 -Wall -Wno-deprecated -g
#GLPK_LIB = -L /usr/local/lib
#GLPK_INC = -I /usr/local/include
SBML_LIB = -L /usr/local/lib
SBML_INC = -I /usr/local/include/sbml
XERCES_LIB = -L /usr/include/xercesc/lib
XERCES_INC = -I /usr/include/xercesc/include
OOQP_LIB = -L /home/karthik/OOQP/OOQP/lib
OOQP_INC = -I /home/karthik/OOQP/OOQP/include

LIB_OPTIONS = ${GLPK_LIB} ${SBML_LIB} ${XERCES_LIB} ${OOQP_LIB}
INC_OPTIONS = ${GLPK_INC} ${SBML_INC} ${XERCES_INC} ${OOQP_INC}

all:	PA_MoMA PA_FBA getStoich

PA_FBA:	sp_linprog.o PA_FBA.h PA_FBA.C
	${CC} ${OPTFLAGS} PA_FBA.C ${INC_OPTIONS} ${LIB_OPTIONS} sp_linprog.o -lglpk /usr/local/lib/libsbml.so -lxerces-c ${OFLAGS} ../bin/PA_FBA

PA_MoMA: sp_linprog.o PA_MoMA.h PA_MoMA.C
	${CC} ${OPTFLAGS} PA_MoMA.C ${INC_OPTIONS} ${LIB_OPTIONS} sp_linprog.o -lglpk -o ../bin/PA_MoMA -lsbml -looqpgensparse -looqpsparse  -looqpgondzio -looqpbase -lblas /home/karthik/OOQP/OOQP/extras/MA27/libMA27.a -lg2c -lgfortran

sp_linprog.o:	sp_linprog.C sp_linprog.h 
	${CC} ${OPTFLAGS} ${CFLAGS} ${INC_OPTIONS} ${LIB_OPTIONS} sp_linprog.C

getStoich:	getStoichiometricMatrix.cpp
	${CC} ${OPTFLAGS} getStoichiometricMatrix.cpp -lsbml -lxerces-c -O3 ${INC_OPTIONS} ${LIB_OPTIONS} -o ../bin/getStoich

clean: 
	rm -f *.o 
	rm -f ../bin/PA_FBA ../bin/PA_MoMA 
	rm -f ../bin/getStoich
