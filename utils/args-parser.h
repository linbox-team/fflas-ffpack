/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/tests/test-common.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified 2011 Brice Boyer (more types,...)
 *
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 27, 2002.
 *
 * Added parametrization to the VectorCategory tags to make them fit the
 * Rootbeer meeting design of VectorCategories being parametrized by
 * VectorTraits.
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __FFLAFLAS_test_common_H
#define __FFLAFLAS_test_common_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "fflas-ffpack/utils/print-utils.h"


using namespace std;

enum ArgumentType {
	TYPE_NONE, TYPE_INT, TYPE_INTEGER, TYPE_DOUBLE, TYPE_INTLIST, TYPE_STR
};
#define TYPE_BOOL TYPE_NONE

#define END_OF_ARGUMENTS \
     { '\0', "\0", "\0", TYPE_NONE, NULL }

struct Argument
{
	char             c;
	const char      *example    ;
	const char      *helpString ;
	ArgumentType     type       ;
	void            *data       ;
};
// example may be passed as null and will be generated intelligently
// eg "-b {YN+-}" for bools, "-v v" for all else


void parseArguments (int argc, char **argv, Argument *args, bool printDefaults = true);

/** writes the values of all arguments, preceded by the programName */
std::ostream& writeCommandString (std::ostream& os, Argument *args, char* programName);


/* Give an approximation of the value of the incomplete gamma function at a, x,
 * to within the tolerance tol */



void printHelpMessage (const char *program, Argument *args, bool printDefaults = false)
{
	int i, l;

	// Skip past libtool prefix in program name
	if (!strncmp (program, "lt-", strlen ("lt-")))
		program += strlen ("lt-");

	std::cout << "Usage: " << program << " [options] [<report file>]" << std::endl;
	std::cout << std::endl;
	std::cout << "Where [options] are the following:" << std::endl;

	for (i = 0; args[i].c != '\0'; ++i) {
		if (args[i].example != 0) {
			std::cout << "  " << args[i].example;
			l = 10 - strlen (args[i].example);
			do std::cout << ' '; while (--l > 0);
		}
		else if (args[i].type == TYPE_NONE)
			std::cout << "  -" << args[i].c << " {YN+-} ";
		else
			std::cout << "  -" << args[i].c << ' ' << args[i].c << "      ";

		std::cout << args[i].helpString;
		if (printDefaults) {
			l = 54 - strlen (args[i].helpString);
			do std::cout << ' '; while (--l > 0);
			std::cout << " (default ";
			switch (args[i].type) {
			case TYPE_NONE:
				cout << ((*(bool *)args[i].data)?"ON":"OFF");
				break;
			case TYPE_INT:
				cout << *(int *) args[i].data;
				break;
			case TYPE_INTEGER:
				cout << *(long int *) args[i].data;
				break;
			case TYPE_DOUBLE:
				cout << *(double *) args[i].data;
				break;
			case TYPE_INTLIST:
				cout << *(std::list<int> *) args[i].data ;
			case TYPE_STR:
				cout << *(std::string *) args[i].data ;
			}
			std::cout << ")";
		}
		std::cout << std::endl;
	}

	std::cout << "  -h or -?  Display this message" << std::endl;
	std::cout << "For boolean switches, the argument may be omitted, meaning the switch should be ON" << std::endl;
	std::cout << std::endl;
	std::cout << "If <report file> is '-' the report is written to std output.  If <report file> is" << std::endl;
	std::cout << "not given, then no detailed reporting is done. This is suitable if you wish only" << std::endl;
	std::cout << "to determine whether the tests succeeded." << std::endl;
	std::cout << std::endl;
	std::cout << "[1] N.B. This program does not verify the primality of Q, and does not use a" << std::endl;
	std::cout << "    field extension in the event that Q=p^n, n > 1" << std::endl;
	std::cout << std::endl;
}

/* Find an argument in the argument list for a character */

Argument *findArgument (Argument *args, char c)
{
	int i;

	for (i = 0; args[i].c != '\0' && args[i].c != c; ++i) ;

	if (args[i].c != '\0')
		return &(args[i]);
	else
		return (Argument *) 0;
}

/* Parse command line arguments */

/*! @internal
 *  @brief transforms a string list of ints to a list of int
 *  string "12,13,15"  is turned into list of ints {12,13,15}
 *  @param outlist list once converted
 *  @param instring list to be converted
 *  @return status message.
 */
int getListArgs(std::list<int> & outlist, std::string & instring)
{
	int start = 0 ;
	int count = 0 ;
	size_t i = 0 ;
	for( ; i < instring.size() ; ++i) {
		if (isdigit(instring[i])) {
			++count;
			continue ;
		}
		if (ispunct(instring[i])) {
			if (!count) {
				std::cout  << std::endl << "ill formed list " << instring << std::endl;
				for (size_t sp = 0 ; sp < 16+i ; ++sp)
					std::cout << '-' ;
				std::cout << '^' << std::endl;
				return(1);
			}
			int j = atoi(instring.substr(start,count).c_str());
			outlist.push_front(j);
			count =  0 ;
			start = i+1 ;
		}
		else {
			std::cout << std::endl << "ill formed list " << instring << std::endl;
			for (size_t sp = 0 ; sp < 16+i ; ++sp)
				std::cout << '-' ;
			std::cout << '^' << std::endl;
			return(1);
		}

	}
	std::cout << std::endl;
	if (!count) {
		std::cout  << std::endl << "ill formed list " << instring << std::endl;
		for (size_t sp = 0 ; sp < 15+i ; ++sp)
			std::cout << '-' ;
		std::cout << '^' << std::endl;
		return(1);
	}

	int j = atoi(instring.substr(start,count).c_str());
	outlist.push_front(j);

	return 0 ;
}


void parseArguments (int argc, char **argv, Argument *args, bool printDefaults)
{
	int i;
	Argument *current;

	for (i = 1; i < argc; ++i) {
		if (argv[i][0] == '-') {
			if (argv[i][1] == 0) {
				// LinBox::commentator.setBriefReportStream (cout);
				// LinBox::commentator.setReportStream (cout);
				std::cout << "Writing report data to cout (intermingled with brief report)" << std::endl << std::endl;
				std::cout.flush ();
			}
			else if (argv[i][1] == 'h' || argv[i][1] == '?') {
				printHelpMessage (argv[0], args, printDefaults);
				exit (1);
			}
			else if ((current = findArgument (args, argv[i][1])) != (Argument *) 0) {
				switch (current->type) {
				case TYPE_NONE:
					{
					if (argc == i+1 || (argv[i+1][0] == '-' && argv[i+1][1] != '\0')) {
						// if at last argument, or next argument is a switch, set to true
						*(bool *) current->data = true;
						break;
					}
					*(bool *) current->data =
					(argv[i+1][0] == '+'
					 || argv[i+1][0] == 'Y'
					 || argv[i+1][0] == 'y'
					 || argv[i+1][0] == 'T'
					 || argv[i+1][0] == 't') ;
					++i;
				}
					break;

				case TYPE_INT:
					{
					*(int *) current->data = atoi (argv[i+1]);
					++i;
					}
					break;

				case TYPE_INTEGER:
					{
						long int tmp = atoi(argv[i+1]);
						*(long int *) current->data = tmp;
					}
					++i;
					break;

				case TYPE_DOUBLE:
					{
						*(double *) current->data = atof (argv[i+1]);
						++i;
					}
					break;

				case TYPE_INTLIST:
					{
						std::string lst = argv[i+1] ;
						std::list<int> LST ;
						getListArgs(LST,lst);
						*(std::list<int> *) current->data = LST ;
						++i;
					}
					break;

				case TYPE_STR:
					*(std::string *) current->data = argv[i+1] ;
					break;

				}
			} else {
				std::cerr << "ERROR: Bad argument " << argv[i] << std::endl;
				break;
			}
		} else {
			// LinBox::commentator.setBriefReportStream(cout);
			// LinBox::commentator.setDefaultReportFile (argv[i]);
			std::cout << "Writing report data to " << argv[i] << std::endl << std::endl;
			std::cout.flush ();
		}
	}
}

std::ostream& writeCommandString (std::ostream& os, Argument *args, char* programName)
{
	os << programName;
	for (int i = 0; args[i].c != '\0'; ++i) {
		cout << " -" << args[i].c;
		switch (args[i].type) {
		case TYPE_NONE:
			if (! (*(bool *)args[i].data)) os << " N";
			break;
		case TYPE_INT:
			os << ' ' << *(int *) args[i].data;
			break;
		case TYPE_INTEGER:
			os << ' ' << *(long int *) args[i].data;
			break;
		case TYPE_DOUBLE:
			os << ' ' << *(double *) args[i].data;
			break;
		case TYPE_INTLIST:
			os << ' ' << *(std::list<int> *) args[i].data;
			break;
		case TYPE_STR:
			os << ' ' << *(std::string *) args[i].data;
			break;
		}
	}
	return os << std::endl;
}


#endif // __FFLAFLAS_test_common_H
