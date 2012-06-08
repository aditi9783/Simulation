//
//  main.cpp
//  cpp_code
//
//  Created by Aditi Gupta on 5/29/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include "boost/regex.hpp"
#include "snp.h"
#include <stdio.h>
using namespace std;

static int usage() {
    cout << "USAGE:\n"
        << "-r\tHapread file\n"
        << "-s\tSNP positions as obtained by first error pass (redErr.pl)\n"
        << "-h\tHaplotype file\n";
    return 0;
}

void hapread_open(const char *);
void snp_open(const char *);
void hapfile_open(const char *);

int main(int argc, const char * argv[])
{
    // usage();
    cout << "argc is" << argc;
    
    // check if enough parameters have been passed
    if (argc < 3)
        usage();
    else { 
        for(int i = 1; i < argc; i++) {
            if ( strcmp(argv[i], "-s") == 0 ) {
                snp_open(argv[i+1]);
            } else if ( strcmp(argv[i], "-h") == 0 ) {
                hapfile_open(argv[i+1]);
            } else if ( strcmp(argv[i], "r") == 0 ) {
                hapread_open(argv[i+1]);
            }
        }        
    }
    return 0;
}

void snp_open(const char * snpfile) {
    FILE * file;
    int i = 1;
    float pos, nreads, cov, frac;
    char base;
    file = fopen (snpfile , "r");
    if (file == NULL) perror ("Error opening snpfile");
    else {
        while ( ! feof (file) ) {
            snp s;
            fscanf(file, "%d\t%c\t%d\t%d\n", &pos, &base, &nreads, &cov);
            frac = nreads/cov;
            s.define_snp(pos, base, frac);
            cout << "line " << i << ": " << s.position() << ", " << s.base() << ", " << s.cov() << endl;
            i++;
        }
        fclose (file);
    }
}
