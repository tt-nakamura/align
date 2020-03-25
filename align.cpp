#include<string>
#include<iostream>
#include<fstream>
#include "Mat.h"
using namespace std;

#define GAP_CHAR   '-'
#define ALIGN_CHAR ':'

int GlobalAlignment(string& u, string& v,
                    const string& s, const string& t,
                    const Mat<int>& score,
                    int e, int o)
// input:
//   s,t = biological sequence (protein or DNA)
//   score = matching score between i-th and j-th character
//           (symmetric matrix)
//   e = gap extension penalty
//   o = gap open penalty
// output:
//   u,v = globally aligned sequence with max score
//         (including gaps) by Needleman-Wunsch method
// return:
//   total score of alignment u,v
// reference:
//   N. C. Jones and P. A. Pevzner
//    "An Introduction to Bioinformatics Algorithms" section 6.6
{
    int i,j,k,l,a,b;
    Mat<int> A(0, s.length()+1, t.length()+1), S(A), T(A), B(A);
    T[1][0] = A[1][0] = S[0][1] = A[0][1] = -o-e;
    for(i=2; i<=s.length(); i++) T[i][0] = A[i][0] = A[i-1][0] - e;
    for(j=2; j<=t.length(); j++) S[0][j] = A[0][j] = A[0][j-1] - e;
    for(i=1; i<=s.length(); i++) {
        for(j=1; j<=t.length(); j++) {
            a = S[i-1][j] - e;
            b = A[i-1][j] - e - o;
            if(a>b) { S[i][j] = a; k=4; }
            else { S[i][j] = b; k=0; }
            a = T[i][j-1] - e;
            b = A[i][j-1] - e - o;
            if(a>b) { T[i][j] = a; k|=8; }
            else T[i][j] = b;
            A[i][j] = A[i-1][j-1] + score[s[i-1]][t[j-1]];
            if(S[i][j] > A[i][j]) { A[i][j] = S[i][j]; B[i][j] = 1; }
            if(T[i][j] > A[i][j]) { A[i][j] = T[i][j]; B[i][j] = 2; }
            B[i][j] |= k;
        }
    }
    u.clear();
    v.clear();
    i = s.length(); k=0;
    j = t.length(); l=0;
    while(i>0 || j>0) {
        if     (i==0) { k=1; l=0; --j; }
        else if(j==0) { k=0; l=1; --i; }
        else if(l==1 || B[i][j]&1) { l=(B[i][j]&4 ? 1:2); k=0; --i; }
        else if(k==1 || B[i][j]&2) { k=(B[i][j]&8 ? 1:2); l=0; --j; }
        else { k=l=0; --i; --j; }
        u.insert(0, 1, k ? GAP_CHAR : s[i]);
        v.insert(0, 1, l ? GAP_CHAR : t[j]);
    }
    return A[s.length()][t.length()];
}

int LocalAlignment(string& u, string& v,
                   int& u_begin, int& v_begin,
                   const string& s, const string& t,
                   const Mat<int>& score,
                   int e, int o)
// input:
//   s,t = biological sequence (protein or DNA)
//   score = matching score between i-th and j-th character
//           (symmetric matrix)
//   e = gap extension penalty
//   o = gap open penalty
// output:
//   u_begin,v_begin = index of where alignment starts
//   u,v = locally aligned sequence with max score
//         (including gaps) by Smith-Waterman method
// return:
//   total score of alignment u,v
// reference:
//   N. C. Jones and P. A. Pevzner
//    "An Introduction to Bioinformatics Algorithms" section 6.8
{
    int i,j,k,l,a,b;
    Mat<int> A(0, s.length()+1, t.length()+1), S(A), T(A), B(A);
    for(i=1; i<=s.length(); i++) {
        for(j=1; j<=t.length(); j++) {
            a = S[i-1][j] - e;
            b = A[i-1][j] - e - o;
            if(a>b) { S[i][j] = a; k=4; }
            else { S[i][j] = b; k=0; }
            a = T[i][j-1] - e;
            b = A[i][j-1] - e - o;
            if(a>b) { T[i][j] = a; k|=8; }
            else T[i][j] = b;
            a = A[i-1][j-1] + score[s[i-1]][t[j-1]];
            if(S[i][j] > a) { a = S[i][j]; B[i][j] = 1; }
            if(T[i][j] > a) { a = T[i][j]; B[i][j] = 2; }
            A[i][j] = (a<0 ? 0:a);
            B[i][j] |= k;
        }
    }
    a=0;
    for(k=s.length(); k>0; k--)
        for(l=t.length(); l>0; l--)
            if(A[k][l] > a) a = A[i=k][j=l];
    u.clear();
    v.clear();
    while(A[i][j]) {
        if     (l==1 || B[i][j]&1) { l=(B[i][j]&4 ? 1:2); k=0; --i; }
        else if(k==1 || B[i][j]&2) { k=(B[i][j]&8 ? 1:2); l=0; --j; }
        else { k=l=0; --i; --j; }
        u.insert(0, 1, k ? GAP_CHAR : s[i]);
        v.insert(0, 1, l ? GAP_CHAR : t[j]);
    }
    u_begin = i;
    v_begin = j;
    return a;
}

void read_fasta(string& s, const char *filename) {
    ifstream ifs(filename);
    string t;
    s.clear();
    while(!ifs.eof()) {
        getline(ifs, t);
        if(t[0]!='>') s += t;
    }
}

void PrintAlignment(const string& u, const string& v, int width)
// u,v = output of GlobalAlignment or LocalAlignment
// width = number of characters in one line
{
    int i,j,k;
    for(k=0; k<u.length(); k+=width) {
        j = (k+width < u.length() ? k+width : u.length());
        cout << u.substr(k,j-k) << '\n';
        for(i=k; i<j; i++)
            cout << (u[i]==v[i] ? ALIGN_CHAR : ' ');
        cout << '\n';
        cout << v.substr(k,j-k) << '\n';
        if(j<u.length()) cout << '\n';
    }
}

void PlotPath(const char *filename,
              const string& u, const string& v,
              int u_begin, int v_begin)
// write txt file to be loaded by gnuplot
// u,v = output of GlobalAlignment or LocalAlignment
// u_begin, v_begin = output of LocalAlignment
{
    int i(u_begin), j(v_begin), k;
    ofstream ofs(filename);
    for(k=0; k<u.length(); k++) {
        if(u[k]!=GAP_CHAR) i++;
        if(v[k]!=GAP_CHAR) j++;
        ofs << i << ' ' << j << '\n';
    }
}

int CalcScore(const string& u, const string& v,
              const Mat<int>& score,
              int e, int o)
// total score of alignment u,v
{
    int i,k;
    for(i=k=0; i<u.length(); i++) {
        if(u[i]==GAP_CHAR) {
            if(i && u[i-1]==GAP_CHAR) k -= e;
            else k -= o+e;
        }
        else if(v[i]==GAP_CHAR) {
            if(i && v[i-1]==GAP_CHAR) k -= e;
            else k -= o+e;
        }
        else k += score[u[i]][v[i]];
    }
    return k;
}