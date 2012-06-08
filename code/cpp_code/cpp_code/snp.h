//
//  snp.h
//  cpp_code
//
//  Created by Aditi Gupta on 6/6/12.
//  Copyright (c) 2012 Michigan State University. All rights reserved.
//

#ifndef cpp_code_snp_h
#define cpp_code_snp_h


class snp{
    int pos;
    char nt;
    float frac_cov;
    int rank;
    int hap;
    
public:
    void define_snp(int, char, float); 
    int position() { return pos; }
    char base() { return nt; }
    float cov() { return frac_cov; }
    void setrank(int);
    int getrank() { return rank; }
    void hap_assoc(int);
    int gethaps();
    
};


#endif
