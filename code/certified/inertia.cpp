// Deterministic symmetric fraction-free elimination. Input: n, then n*n integers.
// Congruence signs follow from successive Bareiss pivots. Zero diagonal pivots
// are resolved by a simultaneous row/column permutation. A zero diagonal with
// a nonzero off-diagonal requires a 2x2 pivot and is rejected, never guessed.
#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>
int main(){try {
 size_t n; if(!(std::cin>>n)||n==0||n>10000)throw std::runtime_error("bad dimension");
 std::vector<mpz_class>A(n*n);for(auto &a:A)if(!(std::cin>>a))throw std::runtime_error("short input");
 for(size_t i=0;i<n;i++)for(size_t j=0;j<i;j++)if(A[i*n+j]!=A[j*n+i])throw std::runtime_error("not symmetric");
 mpz_class old=1,piv,tmp,num; int neg=0,pos=0,zero=0,swaps=0; std::vector<int>signs;
 for(size_t k=0;k<n;k++){
  size_t p=k;while(p<n && A[p*n+p]==0)p++;
  if(p==n){bool allzero=true;for(size_t i=k;i<n;i++)for(size_t j=k;j<n;j++)if(A[i*n+j]!=0)allzero=false;
   if(!allzero)throw std::runtime_error("2x2 pivot required");zero=n-k;break;}
  if(p!=k){for(size_t j=0;j<n;j++)std::swap(A[k*n+j],A[p*n+j]);for(size_t i=0;i<n;i++)std::swap(A[i*n+k],A[i*n+p]);swaps++;}
  piv=A[k*n+k];int s=mpz_sgn(piv.get_mpz_t())*mpz_sgn(old.get_mpz_t());signs.push_back(s);if(s<0)neg++;else pos++;
  for(size_t i=k+1;i<n;i++)for(size_t j=i;j<n;j++){
   num=piv*A[i*n+j]-A[i*n+k]*A[k*n+j];
   if(!mpz_divisible_p(num.get_mpz_t(),old.get_mpz_t()))throw std::runtime_error("nonexact Bareiss division");
   mpz_divexact(tmp.get_mpz_t(),num.get_mpz_t(),old.get_mpz_t());A[i*n+j]=tmp;A[j*n+i]=tmp;
  }
  for(size_t i=k+1;i<n;i++)A[i*n+k]=A[k*n+i]=0;
  old=piv;
 }
 std::cout<<"{\"n\":"<<n<<",\"positive\":"<<pos<<",\"negative\":"<<neg<<",\"zero\":"<<zero<<",\"symmetric_swaps\":"<<swaps<<",\"exact_division_checked\":true,\"pivot_signs\":[";
 for(size_t i=0;i<signs.size();i++)std::cout<<(i?",":"")<<signs[i];std::cout<<"],\"status\":\"PASS\"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;} }
