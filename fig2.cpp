// reference:
//  http://www.bioinformaticsonline.org/pr/pr.html?c=03
// data from:
//  https://www.ncbi.nlm.nih.gov/genbank/

#include<iostream>
#include "align.h"

main() {
    int i,j,n,m,e(2),o(2);
    std::string s,t,u,v;

    read_fasta(s, "CAA35120.1.fasta");
    read_fasta(t, "P17678.1.fasta");

    Mat<int> score = PAM250();
    n = GlobalAlignment(u,v,s,t,score,e,o);
    PrintAlignment(u,v,48);
    m = CalcScore(u,v,score,e,o);
    std::cout << n << ' ' << m << '\n';
}