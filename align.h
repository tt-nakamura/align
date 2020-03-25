#include<string>
#include "Mat.h"

Mat<int> eye();
Mat<int> PAM250();

int GlobalAlignment(std::string& u, std::string& v,
                    const std::string& s, const std::string& t,
                    const Mat<int>& score = eye(),
                    int e = 0,// gap extension penalty
                    int o = 0);// gap open penalty

int LocalAlignment(std::string& u, std::string& v,
                   int& u_begin, int& v_begin,// index of beginning
                   const std::string& s, const std::string& t,
                   const Mat<int>& score = eye(),
                   int e = 0,// gap extension penalty
                   int o = 0);// gap open penalty

void read_fasta(std::string& s, const char *filename);

void PrintAlignment(const std::string& u, const std::string& v,
                    int width = 64);

void PlotPath(const char *filename,
              const std::string& u, const std::string& v,
              int u_begin=0, int v_begin=0);

int CalcScore(const std::string& u, const std::string& v,
              const Mat<int>& score= eye(),
              int e = 0,// gap extension penalty
              int o = 0);// gap open penalty