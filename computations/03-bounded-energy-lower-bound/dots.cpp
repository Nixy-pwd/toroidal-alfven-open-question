#include <cstdint>
#include <vector>
extern "C" int rounded_residual(int rows,int cols,const int64_t* starts,
 const int64_t* indices,const int64_t* vals,const int64_t* witness,
 const int64_t* target,int tshift,int rshift,int64_t* result) {
  using wide=__int128_t;
  std::vector<wide> accum(cols);
  const wide tscale=wide(1)<<tshift, rscale=wide(1)<<rshift;
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++) accum[j]=wide(target[(int64_t)i*cols+j])*tscale;
    for(int64_t k=starts[i];k<starts[i+1];k++){
      const auto *w=witness+indices[k]*cols;
      const auto v=vals[k];
      for(int j=0;j<cols;j++) accum[j]-=wide(v)*w[j];
    }
    for(int j=0;j<cols;j++) {
      wide m=accum[j]/rscale;
      if(m>INT64_MAX || m<INT64_MIN) return 1;
      result[(int64_t)i*cols+j]=int64_t(m);
    }
  }
  return 0;
}
