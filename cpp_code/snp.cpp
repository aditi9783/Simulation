//
//  snp.cpp
//  cpp_code
//
//  Created by Aditi Gupta on 6/6/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//

#include <iostream>
#include "main.h"
#include "snp.h"

void snp::define_snp (int a, char b, float c) {     // set values for snp
    pos = a;
    nt = b;
    coverage = c;
}
    

void snp::setrank (int a) {
    rank = a;
}

void snp::hap_assoc (int a) {
    
}
    int* gethaps();
    


