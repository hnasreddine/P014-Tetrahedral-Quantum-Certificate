#include <gmpxx.h>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
int main(int argc,char** argv){ifstream in(argv[1]);int k;in>>k;bool all=true;
for(int z=0;z<k;z++){
string name;int n;in>>name>>n;
vector<vector<mpq_class>> A(n, vector<mpq_class>(n));
for(int i=0;i<n;i++)for(int j=0;j<n;j++){string s;in>>s;A[i][j]=mpq_class(s);}
// symmetric LDL, fraction-accurate Schur complement update
bool ok=true;double minp=1e12;
for(int p=0;p<n;p++){
 mpq_class pivot=A[p][p]; if(pivot<=0){ok=false;break;} if(pivot.get_d()<minp)minp=pivot.get_d();
 for(int i=p+1;i<n;i++)if(A[i][p]!=0){mpq_class factor=A[i][p]/pivot; for(int j=i;j<n;j++) A[j][i] -= factor * A[j][p]; }
 // mirror lower triangular to row p? Correction: updates using A[j][p] prior untouched
}
cout<<name<<" ";cout<<(ok?"EXACT_INTEGER_LDL_PASS":"FAIL")<<" lowerPivotBoundApprox "<<minp<<"\n";all &=ok;
}
return all?0:1;}
