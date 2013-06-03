#include <math.h> 
#include <iostream>
#include <mpi.h>
#include <unistd.h>
#include <string>
#include <time.h>
#include <algorithm> 
#include <functional> 
#include <locale>
#include <sys/utsname.h>

// trim from start
static inline std::string ltrim(std::string const &_s) {
        std::string s = _s;
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string rtrim(std::string const &_s) {
        std::string s = _s;
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string trim(std::string const &_s) {
        return ltrim(rtrim(_s));
}

 
int main(int argc, char *argv[]) 
{ 
  int n, rank, size, i, s; 
  const double PI25DT = 3.1415926535897932384626433832795028841971693993751; 
  double mypi, pi, h, sum, x; 
 
  MPI::Init(argc, argv); 
  size = MPI::COMM_WORLD.Get_size(); 
  rank = MPI::COMM_WORLD.Get_rank(); 
 
  n = 1;
  s = 0;
  if (rank == 0)
    for(i = 0; i < argc-1; ++i)
    {
      if(trim(argv[i]) == "--order") { n = atoi(argv[i+1]); }
      else if(trim(argv[i]) == "--sleep") { s = atoi(argv[i+1]); }
    }
  MPI::COMM_WORLD.Bcast(&n, 1, MPI::INT, 0); 
  MPI::COMM_WORLD.Bcast(&s, 1, MPI::INT, 0); 

  MPI::COMM_WORLD.Bcast(&n, 1, MPI::INT, 0); 
  h = 1.0 / (double) std::max(1, n); 
  sum = 0.0; 
  for (i = rank + 1; i <= n; i += size) { 
    x = h * ((double)i - 0.5); 
    sum += (4.0 / (1.0 + x*x)); 
    if(s != 0) ::sleep(s);
  } 
  mypi = h * sum; 
 
  MPI::COMM_WORLD.Reduce(&mypi, &pi, 1, MPI::DOUBLE, 
  		   MPI::SUM, 0); 
  if (rank == 0) 
  {
    std::cout << "pi to order " << n << " is approximately " << std::scientific << pi 
              << ", Error is " << fabs(pi - PI25DT) 
              << "  -- slept " << s << " seconds at each iteration -- "
              << " mpi world size is " << size
              << std::endl; 
    if(n == 666) 
    {
      // cray crapware does not check exit code when launched with more than
      // one process. hence, can't MPI::Finalize(); 
      return 2;
    }
    utsname sys;
    uname(&sys);
    std::cout << "sysname: " << sys.sysname << "\n";
    std::cout << "nodename: " << sys.nodename << "\n";
    std::cout << "release: " << sys.release << "\n";
    std::cout << "version: " << sys.version << "\n";
    std::cout << "machine: " << sys.machine << "\n";
  }
  if (n == 6666) return 2;
  MPI::Finalize(); 
  return n == 6666 ? 2: 0; 
} 

