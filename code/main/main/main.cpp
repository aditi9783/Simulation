//
//  main.cpp
//  main
//
//  Created by Aditi Gupta on 5/24/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main()
{
    
    string line;
    ofstream data;
    ifstream seqfile;
    
    seqfile.open("rangenome10k.fa");
    data.open("duplicate.fa");
        
    if (seqfile.is_open()) {
        while (seqfile.good()) {
            getline(seqfile, line);
            data << line << endl;
        }
    } else {
        cout << "cannot open file\n";
    }
    
    seqfile.close();
    data.close();
    
    return 0;
 }
