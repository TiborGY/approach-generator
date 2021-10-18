#include <memory>
#include <random>
#include <iomanip>
#include <algorithm>

inline void pauser(){
    /// a portable way to pause a program
    std::string dummy;
    std::cout << "Press enter to continue...";
    std::getline(std::cin, dummy);
}

inline int64_t getFilesize(const std::string& fname){
    /// return with the size of a file in bytes, return -1 for error
    std::ifstream f;
    f.open(fname, std::ios::in | std::ios::binary); //open the file and check for error
    if (!f.is_open()){
       std::cout<< "error: open file for size check failed!" <<std::endl;
       std::cout<< "Cannot open file: " << fname << std::endl;
       return -1;
    }
    f.seekg (0, f.end); //getting the size of the file
    return f.tellg(); //no need to close f, the destructor handles it
}

inline bool FileExists(const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {return false;}
}

inline int CallBash(const std::string& cmd){
  const int ret = std::system(cmd.c_str());
  #pragma omp critical
  std::cout << "Return value of " << cmd << " : " << ret << std::endl;
  return ret;
}

inline int CallBash_quiet(const std::string& cmd){
	const int ret = std::system(cmd.c_str());
	#pragma omp critical
	if (ret!=0) std::cout << "Return value of " << cmd << " : " << ret << std::endl;
	return ret;
}

std::string GetCommandOutput(const std::string& cmdstr) {
    ///based on https://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c-using-posix
    const char * const cmd = cmdstr.c_str();
    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe){
        std::cout<<"popen() failed!"<<std::endl;
        assert(pipe);
    }
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
    return result;
}

template <typename T>
std::string to_string_prec(const T a_value, const int prec = 8){
    std::ostringstream out;
    out << std::setprecision(prec) << std::fixed << a_value;
    return out.str();
}

template <typename T>
std::string getSpacer(const T value, const uint32_t maxWidth){
	std::string spacer;
	for(uint32_t i=1; i<maxWidth; i++){
		if (value < std::pow(10, i)) spacer += ' ';
	}
	return spacer;
}

std::mt19937_64 initPRNG(){
  //prepare a uniform PRNG
  std::random_device rd;
  std::mt19937_64 PRNG(rd()); //seed PRNG using /dev/urandom or similar OS provided RNG
  std::uniform_int_distribution<> tmpdist{0, 255};
  //make sure the internal state of the PRNG is properly mixed by generating 10M random numbers
  //PRNGs often have unreliable distribution uniformity and other statistical properties before their internal state is sufficiently mixed
  int tmprngint;
  for (uint32_t i=0;i<10000000;i++) tmprngint=tmpdist(PRNG);
  return PRNG;
}

inline bool fpeq(const double a, const double b){
	if(std::abs(a-b) < 1.0E-14) return true;
	return false;
}

inline bool fpwhole(const double a){
	return fpeq(std::round(a), a);
}

inline std::string prettybool(const bool a){
	if (a){return "YES";}else{return "NO";}
}

std::vector<std::string> getCLIargs(const int argc, const char * const argv[]){
	std::vector<std::string> cliArgs;
	for (int i=1; i<argc; i++) cliArgs.push_back(argv[i]);
	return cliArgs;
}

void saveUint(const std::string fname, const uint32_t i){
	std::ofstream outfile;
	outfile.open(fname);
	if (!outfile.is_open()){
		std::cout<< "Error: could not open "<< fname<< " for writing!"<< std::endl;
		abort();
	}
	outfile<< i<< std::endl;
	if (!outfile.good()){
		std::cout<< "Error: failure while writing to "<< fname<< " !"<< std::endl;
		abort();
	}
}

void saveDouble(const std::string fname, const double d, const int prec = 8){
	std::ofstream outfile;
	outfile.open(fname);
	if (!outfile.is_open()){
		std::cout<< "Error: could not open "<< fname<< " for writing!"<< std::endl;
		abort();
	}
	outfile << std::setprecision(prec) << std::fixed;
	outfile<< d<< std::endl;
	if (!outfile.good()){
		std::cout<< "Error: failure while writing to "<< fname<< " !"<< std::endl;
		abort();
	}
}

uint32_t loadUint(const std::string fname){
	uint32_t tmp;
	std::ifstream infile;
	infile.open(fname);
	if (!infile.is_open()){
		std::cout<< "Error: could not open "<< fname<< " for reading!"<< std::endl;
		abort();
	}
	infile >> tmp;
	return tmp;
}

double loadDouble(const std::string fname){
	double tmp;
	std::ifstream infile;
	infile.open(fname);
	if (!infile.is_open()){
		std::cout<< "Error: could not open "<< fname<< " for reading!"<< std::endl;
		abort();
	}
	infile >> tmp;
	return tmp;
}

std::vector<double> loadDouble(const std::string fname, const uint32_t count){
	std::vector<double> tmp(count);
	std::ifstream infile;
	infile.open(fname);
	if (!infile.is_open()){
		std::cout<< "Error: could not open "<< fname<< " for reading!"<< std::endl;
		abort();
	}
	for(uint32_t i=0; i<count; i++) infile >> tmp.at(i);
	return tmp;
}

inline double square(const double a){
	return a*a;
}

inline double calcScalingFactor(const double energy, const double E_tol, const double maxPotE){
	/// This functtion expexts all three energies to be given in "absolute" energies.
	if (energy > maxPotE){
		return (E_tol + energy - maxPotE)/E_tol;
	}else{
		return 1.0;
	}
}

inline double calcDwtWeight(const double energy, const double Edwt0, const double Edwt1){
	/// This function expects all three energies to be given as energies relative the current minimum
	/// of the fitting set.
	return ( Edwt0/(energy+Edwt0) )*( Edwt1/(energy+Edwt1) );
}

inline double SUM(const std::vector<double>& input){
	return std::accumulate(input.cbegin(), input.cend(), 0.0);
}

inline double AVG(const std::vector<double>& input){
	return SUM(input)/input.size();
}

inline double MAX(const std::vector<double>& input){
	return *std::max_element(input.cbegin(), input.cend());
}

inline double RMS(const std::vector<double>& input){
	return std::sqrt(
		std::accumulate(
			input.cbegin(),
			input.cend(),
			0.0,
			[](double squared_sum, const double unsquared_in){
				return squared_sum+square(unsquared_in);
			}
		)
		/input.size()
	);
}

template<typename T>
void readStringOption(std::ifstream& infile, const std::string targetStr, T& output, uint32_t& linectr){
	///read a line of an already opened config file and attempt to extract the value of the option, specified by a given string
	///lines begginging with '#' are considered comments, and skipped
	///also increments the config file line counter, to make finding mistakes in the config file easier
	std::array<char,4096> tmp;
	tmp.fill(0);
	//validate intial state
	assert(targetStr.size() < 4096);
	assert(infile.is_open());
	assert(infile.good());
	//skip lines if they are just comments
	while(infile.peek()=='#'){
		std::string junk;
		std::getline(infile, junk);
		linectr++;
	}
	//extract a given number of characters, convert to string
	infile.getline(&tmp[0], targetStr.size()+1);
	const std::string line(&tmp[0]);
	//clear the failbit set by ifstream::getline(...), as it stops before the end of the line
	infile.clear(infile.rdstate() &~ std::ifstream::failbit); //bitwise AND NOT
	if (line!=targetStr){
		std::cout<<"Error: | line "<<linectr<<" of config file does not start with \""+targetStr+"\""<<std::endl;
		std::cout<<"       | Found this instead: "<<line<<std::endl;
		abort();
	}
	//at this point the remaining characters in the line should be the value of the option
	std::string s;
	std::getline(infile, s);
	linectr++;
	//validate and process the option based on type
	if constexpr (std::is_same<T,std::string>::value){
		output = s;
		return;
	}
	if constexpr (std::is_same<T,bool>::value){
		if ((s=="1")||(s=="true")||(s=="TRUE")||(s=="yes")||(s=="YES")||(s=="y")||(s=="Y")){
			output=true;
			return;
		}
		if ((s=="0")||(s=="false")||(s=="FALSE")||(s=="no")||(s=="NO")||(s=="n")||(s=="N")){
			output=false;
			return;
		}
		abort();
	}
	if constexpr (std::is_same<T,double>::value){
		output = std::stod(s);
		return;
	}
	if constexpr (std::is_same<T,uint32_t>::value){
		output = std::stoul(s);
		return;
	}
	//if this is reached, then an unimplemented type was requested
	abort();
}
