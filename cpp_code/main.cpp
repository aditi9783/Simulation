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

/*
class stringProcessor:
{
private:
    vector<string*> container;
public:
    stringProcessor();
    ~stringProcessor()
    {
        for(int x = 0; x<container.size(); x++)
            if(container[x]) delete container[x];
    }
    void analyzeFile(string fileName)
    {
        ifstream fileIn;
        fileIn.open(fileName.c_str());
        while(!fileIn.eof())
        {
            if(fileIn>>
        }
    }
    
    
};
*/
static int usage() {
    cout << "USAGE:\n"
        << "-r\tHapread file\n"
        << "-s\tSNP positions as obtained by first error pass (redErr.pl)\n"
        << "-h\tHaplotype file\n";
    return 0;
}

//void hapread_open(const char *);
void snp_open(const char *);
//void hapfile_open(const char *);

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
    //        } else if ( strcmp(argv[i], "-h") == 0 ) {
        //        hapfile = argv[i+1];
        //        hapfile_open(hapfile);
        //        cout << "Haplotype file is " << hapfile << endl;
            } 
        }        
    }
    return 0;
}

void snp_open(const char * snpfile) {
    FILE * file;
    int pos, nreads, cov, i = 1;
    char base;
    file = fopen (snpfile , "r");
    if (file == NULL) perror ("Error opening snpfile");
    else {
        while ( ! feof (file) ) {
            fscanf(file, "%d\t%c\t%d\t%d\n", &pos, &base, &nreads, &cov);
            cout << "line " << i << ": " << pos << ", " << base << ", " << nreads << ", " << cov << endl;
            i++;
        }
        fclose (file);
    }
/*    ifstream file (snpfile);
    string line;
    int pos, nreads, cov;
    char base;
    
    if (file.is_open()) {
        boost::regex pattern("(\d+)\t([A-Za-z])\t(\d+)\t(\d+)");
        boost::smatch result;
        while ( file.good() ) {
            getline (file, line);
            stringstream linestream(line);
            bool isMatchFound = boost::regex_match(line, result, pattern); 
            if (isMatchFound) {
//                int value(result[1], result[3], result[4]);
                cout<<result[1]<<endl;
//                base = result[2];
//                nreads = result[3];
//                cov = result[4];
//                for (unsigned int i=1; i < what.size(); i++) 
//                { 
//                    cout << "WHAT " << i << " " << what[i] << endl; 
//                } 
            } 
        }
    } else {
        cout << "cannot open file\n";
    }
*/    
}
