#include <gmpxx.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
using namespace std;
using Q=mpq_class;using Z=mpz_class;
static Q qr(const Z &a,const Z &b) { Q v(a,b);v.canonicalize();return v; }
static bool readZ(istream& is,Z& a){string s; if(!(is>>s))return false;a=Z(s);return true;}
int main(int argc,char** argv){
 if(argc<2){cerr<<"usage certify_inverse_rational data.txt\n";return 2;}
 ifstream f(argv[1]);int n;Z sa,sr,sb,sy;if(!(f>>n)||!readZ(f,sa)||!readZ(f,sr)||!readZ(f,sb)||!readZ(f,sy))return 2;
 Z eAn,eAd,eBn,eBd;readZ(f,eAn);readZ(f,eAd);readZ(f,eBn);readZ(f,eBd);
 vector<vector<pair<int,Z>>> A(n);int entries=0;
 for(int k=0;k<n;k++){int cnt;f>>cnt;entries+=cnt;for(int z=0;z<cnt;z++){int j;Z val;f>>j;readZ(f,val);A[k].push_back({j,val});}}
 vector<vector<Z>> R(n,vector<Z>(n));for(int i=0;i<n;i++)for(int k=0;k<n;k++)readZ(f,R[i][k]);
 vector<Z> b(n),x(n);for(int j=0;j<n;j++)readZ(f,b[j]);for(int j=0;j<n;j++)readZ(f,x[j]);
 cout<<"Exact rational verifier: n="<<n<<", nz="<<entries<<"\n"<<flush;
 // rho0, ||R||inf
 Z maxRowsum(0),maxRsum(0);
 Z scaledidentity=sa*sr;
 for(int i=0;i<n;i++){
  vector<Z> row(n);
  Z rsum=0;
  for(int k=0;k<n;k++){
   const Z &factor=R[i][k];if(factor==0)continue;
   rsum+=abs(factor);
   for(auto &it:A[k]) mpz_addmul(row[it.first].get_mpz_t(),factor.get_mpz_t(),it.second.get_mpz_t());
  }
  maxRsum=max(maxRsum,rsum);
  Z sum=0;
  for(int j=0;j<n;j++) {Z err=(i==j ? scaledidentity:Z(0))-row[j];sum+=abs(err);}
  maxRowsum=max(maxRowsum,sum);
 }
 Q rho0=qr(maxRowsum,sa*sr),Rnorm=qr(maxRsum,sr);
 Q delA=qr(eAn,eAd),delB=qr(eBn,eBd);
 Q rho=rho0+Rnorm*delA;
 cout<<"Certified rho0 = "<<rho0.get_d()<<"; ||Rhat||inf="<<Rnorm.get_d()<<"; delta_A="<<delA.get_d()<<"; rho upper="<<rho.get_d()<<"\n"<<flush;
 // residual z=(b-A*x) and R*z
 if(sa!=sb){cerr<<"Implementation currently requires SA=SB\n";return 2;}
 vector<Z> z(n);
 for(int i=0;i<n;i++){
  z[i]=b[i]*sy;
  for(auto &it:A[i]) mpz_submul(z[i].get_mpz_t(),it.second.get_mpz_t(),x[it.first].get_mpz_t());
 }
 Z maxnu0=0;
 for(int i=0;i<n;i++){
  Z acc=0;
  for(int k=0;k<n;k++)mpz_addmul(acc.get_mpz_t(),R[i][k].get_mpz_t(),z[k].get_mpz_t());
  maxnu0=max(maxnu0,Z(abs(acc)));
 }
 Z xnorm=0;for(auto &el:x)xnorm=max(xnorm,Z(abs(el)));
 Q nu0=qr(maxnu0,sr*sa*sy),xNorm=qr(xnorm,sy);
 Q nu=nu0+Rnorm*(delB+delA*xNorm);
 if(rho>=1){cout<<"FAIL: rho>=1\n";return 1;}
 Q eps=nu/(Q(1)-rho);Q limit(1,1200000);
 cout<<"Certified nu0 = "<<nu0.get_d()<<"; delta_b="<<delB.get_d()<<"; ||xhat||inf="<<xNorm.get_d()<<"\n";
 cout<<"Certified nu upper="<<nu.get_d()<<"; pivot error bound="<<eps.get_d()<<"; limit="<<limit.get_d()<<"\n"<<flush;
 cout<<"rho < 1 exact: "<<(rho<1?"PASS":"FAIL")<<"\n";
 cout<<"max pivot error < positivity margin exact: "<<(eps<limit?"PASS":"FAIL")<<"\n";
 if(!(eps<limit))return 1;
 cout<<"PIVOT CERTIFICATION VERIFIED (rational inequalities).\n";
 return 0;
}
