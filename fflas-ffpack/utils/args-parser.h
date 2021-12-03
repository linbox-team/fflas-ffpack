/* utils/args-parser.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 27, 2002.
 * Modified 2011 Brice Boyer (more types,...)
 *
 * Added parametrization to the VectorCategory tags to make them fit the
 * Rootbeer meeting design of VectorCategories being parametrized by
 * VectorTraits.
 *
 * ------------------------------------
 * ========LICENCE========
 * This file is part of the library FFLAS-FFPACK.
 *
 * FFLAS-FFPACK is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __FFLASFFPACK_args_parser_H
#define __FFLASFFPACK_args_parser_H

#include <fflas-ffpack/fflas-ffpack-config.h>
#include <givaro/givinteger.h>
#include <givaro/givprint.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <list>
#include <stdlib.h>

enum ArgumentType {
    TYPE_NONE, TYPE_INT, TYPE_UINT64, TYPE_LONGLONG, TYPE_INTEGER, TYPE_DOUBLE, TYPE_INTLIST, TYPE_STR
};
#define TYPE_BOOL TYPE_NONE

#define END_OF_ARGUMENTS \
{ '\0', "\0", "\0", TYPE_NONE, NULL }

#ifdef _GIVARO_CONFIG_H
#define type_integer Givaro::Integer
#else
#define type_integer long int
#endif

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

namespace FFLAS {
    void parseArguments (int argc, char **argv, Argument *args, bool printDefaults = true);
}

void printHelpMessage (const char *program, Argument *args, bool printDefaults = false)
{
    int i, l;

    // Skip past libtool prefix in program name
    if (!strncmp (program, "lt-", strlen ("lt-")))
        program += strlen ("lt-");

    std::cout << "Usage: " << program << " [options] [<report file>]" << std::endl;
    std::cout << std::endl;
    std::cout << "Where [options] are the following:" << std::endl;

    bool messageboolean(false),messageprimality(false);

    for (i = 0; args[i].c != '\0'; ++i) {
        if (args[i].example != 0) {
            std::cout << "  " << args[i].example;
            l = 10 - (int)strlen (args[i].example);
            do std::cout << ' '; while (--l > 0);
        }
        else if (args[i].type == TYPE_NONE) {
            std::cout << "  -" << args[i].c << " {YN+-} ";
            messageboolean = true;
        }
        else
            std::cout << "  -" << args[i].c << ' ' << args[i].c << "      ";

        std::cout << args[i].helpString;
        if (strncmp(args[i].helpString,"Operate over the \"field\"",24) == 0)
            messageprimality = true;
        if (printDefaults) {
            l = 54 - (int)strlen (args[i].helpString);
            do std::cout << ' '; while (--l > 0);
            std::cout << " (default ";
            switch (args[i].type) {
            case TYPE_NONE:
                std::cout << ((*(bool *)args[i].data)?"ON":"OFF");
                break;
            case TYPE_INT:
                std::cout << *(int *) args[i].data;
                break;
            case TYPE_UINT64:
                std::cout << *(uint64_t *) args[i].data;
                break;
            case TYPE_LONGLONG:
                std::cout << *(long long *) args[i].data;
                break;
            case TYPE_INTEGER:
                std::cout << *(type_integer *) args[i].data;
                break;
            case TYPE_DOUBLE:
                std::cout << *(double *) args[i].data;
                break;
            case TYPE_INTLIST:
                std::cout << *(std::list<int> *) args[i].data ;
                break;
            case TYPE_STR:
                std::cout << "\"" << *(std::string *) args[i].data << "\"" ;
                break;
            }
            std::cout << ")";
        }
        std::cout << std::endl;
    }

    std::cout << "  -h or -?  Display this message" << std::endl;
    if (messageboolean)
        std::cout << "For boolean switches, the argument may be omitted, meaning the switch should be ON" << std::endl;
    std::cout << std::endl;
    std::cout << "If <report file> is '-' the report is written to std output.  If <report file> is" << std::endl;
    std::cout << "not given, then no detailed reporting is done. This is suitable if you wish only" << std::endl;
    std::cout << "to determine whether the tests succeeded." << std::endl;
    std::cout << std::endl;
    if (messageprimality)
        std::cout << "[1] N.B. This program does not verify the primality of Q, and does not use a" << std::endl
        << "    field extension in the event that Q=p^n, n > 1" << std::endl;
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
            int j = atoi(instring.substr((size_t)start,(size_t)count).c_str());
            outlist.push_front(j);
            count =  0 ;
            start = int(i+1) ;
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

    int j = atoi(instring.substr((size_t)start,(size_t)count).c_str());
    outlist.push_front(j);

    return 0 ;
}


namespace FFLAS {
    /**
     * @brief Get the value of an argument 
     *        and avoid core dump when no value was given after an argument
     * 
     * @param argv argument value list
     * @param i    argument index
     * @return char* argument value
     */
    char* getArgumentValue(int argc, char **argv,int i){
        if (i+1 < argc) {
            return argv[i+1];
        } else {
            std::cerr << "ArgumentParser error: Expected a value after argument " << argv[i] << std::endl;
            exit(-1);
        }
    }

    void parseArguments (int argc, char **argv, Argument *args, bool printDefaults)
    {
        int i;
        Argument *current;

        for (i = 1; i < argc; ++i) {
            // std::cout << "i=" << i << std::endl;
            if (argv[i][0] == '-') {
                if (argv[i][1] == 0) {
                    std::cout << "Writing report data to cout (intermingled with brief report)" << std::endl << std::endl;
                    std::cout.flush ();
                }
                else if (argv[i][1] == 'h' || argv[i][1] == '?' || argv[i][1] == '-') {
                    printHelpMessage (argv[0], args, printDefaults);
                    exit (1);
                }
                else if ((current = findArgument (args, argv[i][1])) != (Argument *) 0) {
                    switch (current->type) {
                    case TYPE_NONE:
                        {
                            if (argc == i+1 || (argv[i+1][0] == '-' && argv[i+1][1] != '\0')) {
                                // if at last argument, or next argument is a switch, set to true
                                *((bool *) current->data) = true;
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
                            *(int *) current->data = atoi (getArgumentValue(argc, argv, i));
                            ++i;
                        }
                        break;

                    case TYPE_UINT64:
                        {
                            *(uint64_t *) current->data = atol (getArgumentValue(argc, argv, i));
                            ++i;
                        }
                        break;

                    case TYPE_LONGLONG:
                        {
                            *(long long *) current->data = atol (getArgumentValue(argc, argv, i));
                            ++i;
                        }
                        break;

                    case TYPE_INTEGER:
                        {
#ifdef _GIVARO_CONFIG_H
                            type_integer tmp(getArgumentValue(argc, argv, i));
#else
                            type_integer tmp = atol(getArgumentValue(argc, argv, i));
#endif
                            *(type_integer *) current->data = tmp;
                        }
                        ++i;
                        break;

                    case TYPE_DOUBLE:
                        {
                            *(double *) current->data = atof (getArgumentValue(argc, argv, i));
                            ++i;
                        }
                        break;

                    case TYPE_INTLIST:
                        {
                            std::string lst = getArgumentValue(argc, argv, i) ;
                            std::list<int> LST ;
                            getListArgs(LST,lst);
                            *(std::list<int> *) current->data = LST ;
                            ++i;
                        }
                        break;

                    case TYPE_STR:
                        {
                            *(std::string *) current->data = getArgumentValue(argc, argv, i);
                            ++i;
                        }
                        break;

                    }
                } else {
                    std::cerr << "ERROR: Bad argument " << argv[i] << std::endl;
                    break;
                }
            } else {
                std::cout << "Writing report data to " << argv[i] << std::endl << std::endl;
                std::cout.flush ();
            }
        }
    }

    /** writes the values of all arguments, preceded by the programName */
    std::ostream& writeCommandString (std::ostream& os, Argument *args, const char* programName = nullptr)
    {
        if (programName != nullptr)
            os << programName;

        for (int i = 0; args[i].c != '\0'; ++i) {
            os << " -" << args[i].c;
            switch (args[i].type) {
            case TYPE_NONE:
                if ((*(bool *)args[i].data)) os << " Y";
                else os << " N";
                break;
            case TYPE_INT:
                os << ' ' << *(int *) args[i].data;
                break;
            case TYPE_UINT64:
                os << ' ' << *(uint64_t *) args[i].data;
                break;
            case TYPE_LONGLONG:
                os << ' ' << *(long long *) args[i].data;
                break;
            case TYPE_INTEGER:
                os << ' ' << *(Givaro::Integer *) args[i].data;
                break;
            case TYPE_DOUBLE:
                os << ' ' << *(double *) args[i].data;
                break;
            case TYPE_INTLIST:
                os << ' ' << *(std::list<int> *) args[i].data;
                break;
            case TYPE_STR:
                os << " \"" << *(std::string *) args[i].data << "\"";
                break;
            }
        }

        return os;
    }
}

#undef type_integer

#endif // __FFLASFFPACK_args_parser_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
