#include "system_types.hpp"

namespace braams_interface{
	bool braams_init_done = false; //"global" C++side flag to track if the Braams PES library is already initialized
}

extern "C"{
	void pes_init_x3y1z1u1(char* coeff_path);
	void pes_init_x4y1z1u1(char* coeff_path);
	void pes_init_x4y2z1(char* coeff_path);
	void pes_init_x5y2z1(char* coeff_path);
	void pes_init_x3y2z1(char* coeff_path);
	//energy evaluation function
	//coordinates are stored in a 1D array
	double pes_energy_x3y1z1u1(double* x);
	double pes_energy_bohr_x3y1z1u1(double* x);
	double pes_energy_x4y1z1u1(double* x);
	double pes_energy_bohr_x4y1z1u1(double* x);
	double pes_energy_x4y2z1(double* x);
	double pes_energy_bohr_x4y2z1(double* x);
	double pes_energy_x5y2z1(double* x);
	double pes_energy_bohr_x5y2z1(double* x);
	double pes_energy_x3y2z1(double* x);
	double pes_energy_bohr_x3y2z1(double* x);
}

void pes_init_str(std::string coeff_path){
	///C++side function for initializing the PES library with the path of fitting coefficients
	assert(coeff_path.size() < 500);
	static_assert(IsSystemBraamsSupported(), "System type is not supported by Braams library.");
	if(braams_interface::braams_init_done == true){
		std::cout<<"ERROR: You tried to call the initialization function of the Braams PES library more than once!"<<std::endl;
		std::cout<<"This is not possible to do, because the PES library has hidden state that cannot be reset from the C++ side."<<std::endl;
		std::cout<<"If you are seeing this message, it means you have run into a bug. Program will now abort."<<std::endl;
		abort();
	}
	braams_interface::braams_init_done = true;
	//this requires C++17
	if constexpr (PES_SYSTEM_TYPE == PES_X3Y1Z1U1){ pes_init_x3y1z1u1(coeff_path.data()); return;}
	if constexpr (PES_SYSTEM_TYPE == PES_X4Y1Z1U1){ pes_init_x4y1z1u1(coeff_path.data()); return;}
	if constexpr (PES_SYSTEM_TYPE == PES_X4Y2Z1){ pes_init_x4y2z1(coeff_path.data()); return;}
	if constexpr (PES_SYSTEM_TYPE == PES_X5Y2Z1){ pes_init_x5y2z1(coeff_path.data()); return;}
	if constexpr (PES_SYSTEM_TYPE == PES_X3Y2Z1){ pes_init_x3y2z1(coeff_path.data()); return;}
	std::cout<<"ERROR: You forgot to add your new system to pes_interface_cppside.hpp"<<std::endl;
	abort();
}

template <auto NUM_COORDS>
void ExtractLinearizedCoords(const xyzfile& input, std::array<double,NUM_COORDS>& output){
	///take an xyzfile object and extract the Cartesian coordinates into a 1D (linearized) array (as opposed to a 3 column 2D array)
	///this variant outputs to an std::array reference
	constexpr uint32_t atomcnt = NUM_COORDS/3;
	static_assert((atomcnt*3) == NUM_COORDS);
	for(uint32_t i=0; i<atomcnt; i++){
		output[i*3]     = input.atomvec[i].xcoord;
		output[(i*3)+1] = input.atomvec[i].ycoord;
		output[(i*3)+2] = input.atomvec[i].zcoord;
	}
}

template <uint32_t NUM_ATOMS>
std::array<double,NUM_ATOMS*3> ExtractLinearizedCoords(const xyzfile& input){
	///take an xyzfile object and extract the Cartesian coordinates into a 1D (linearized) array (as opposed to a 3 column 2D array)
	///this variant outputs by returning an std::array
	std::array<double,NUM_ATOMS*3> output;
	for(uint32_t i=0; i<NUM_ATOMS; i++){
		output[i*3]     = input.atomvec[i].xcoord;
		output[(i*3)+1] = input.atomvec[i].ycoord;
		output[(i*3)+2] = input.atomvec[i].zcoord;
	}
	return output;
}

template <auto NUM_COORDS>
void InsertLinearizedCoords(xyzfile& output, const std::array<double,NUM_COORDS>& input){
	///take a 1D (linearized) std::array of Cartesian coordinates, and insert them into an xyzfile object
	constexpr uint32_t atomcnt = NUM_COORDS/3;
	static_assert((atomcnt*3) == NUM_COORDS);
	for(uint32_t i=0; i<atomcnt; i++){
        output.atomvec[i].xcoord = input[i*3];
        output.atomvec[i].ycoord = input[(i*3)+1];
        output.atomvec[i].zcoord = input[(i*3)+2];
	}
}

template <uint32_t NUM_ATOMS>
double EvaluateEnergy(const xyzfile& coords){
	///return the potential enegy at a point specified by an xyzfile object
	///this function assumes the coordinates are in Angstroms
	std::array<double,NUM_ATOMS*3> x;
	ExtractLinearizedCoords(coords, x);
	static_assert(IsSystemBraamsSupported(), "System type is not supported by Braams library.");
	if constexpr (PES_SYSTEM_TYPE == PES_X3Y1Z1U1) return pes_energy_x3y1z1u1(x.data());
	if constexpr (PES_SYSTEM_TYPE == PES_X4Y1Z1U1) return pes_energy_x4y1z1u1(x.data());
	if constexpr (PES_SYSTEM_TYPE == PES_X4Y2Z1) return pes_energy_x4y2z1(x.data());
	if constexpr (PES_SYSTEM_TYPE == PES_X5Y2Z1) return pes_energy_x5y2z1(x.data());
	if constexpr (PES_SYSTEM_TYPE == PES_X3Y2Z1) return pes_energy_x3y2z1(x.data());
	std::cout<<"ERROR: You forgot to add your new system to pes_interface_cppside.hpp"<<std::endl;
	abort();
}

template <bool geomInBohrs, auto NUM_COORDS>
double EvaluateEnergy(std::array<double,NUM_COORDS>& coords){
	///return the potential enegy at a point specified by a 1D (linearized) std::array of Cartesian coordinates
	static_assert(IsSystemBraamsSupported(), "System type is not supported by Braams library.");
	if constexpr (geomInBohrs){
		if constexpr (PES_SYSTEM_TYPE == PES_X3Y1Z1U1) return pes_energy_bohr_x3y1z1u1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X4Y1Z1U1) return pes_energy_bohr_x4y1z1u1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X4Y2Z1) return pes_energy_bohr_x4y2z1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X5Y2Z1) return pes_energy_bohr_x5y2z1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X3Y2Z1) return pes_energy_bohr_x3y2z1(coords.data());
		std::cout<<"ERROR: You forgot to add your new system to pes_interface_cppside.hpp"<<std::endl;
		abort();
	}else{
		if constexpr (PES_SYSTEM_TYPE == PES_X3Y1Z1U1) return pes_energy_x3y1z1u1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X4Y1Z1U1) return pes_energy_x4y1z1u1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X4Y2Z1) return pes_energy_x4y2z1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X5Y2Z1) return pes_energy_x5y2z1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X3Y2Z1) return pes_energy_x3y2z1(coords.data());
		std::cout<<"ERROR: You forgot to add your new system to pes_interface_cppside.hpp"<<std::endl;
		abort();
	}
}
template <bool geomInBohrs, auto NUM_COORDS>
double EvaluateEnergy(std::array<double,NUM_COORDS>&& coords){
	///return the potential enegy at a point specified by a 1D (linearized) std::array of Cartesian coordinates
	///this version uses an rvalue reference to allow more flexible usage
	static_assert(IsSystemBraamsSupported(), "System type is not supported by Braams library.");
	if constexpr (geomInBohrs){
		if constexpr (PES_SYSTEM_TYPE == PES_X3Y1Z1U1) return pes_energy_bohr_x3y1z1u1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X4Y1Z1U1) return pes_energy_bohr_x4y1z1u1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X4Y2Z1) return pes_energy_bohr_x4y2z1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X5Y2Z1) return pes_energy_bohr_x5y2z1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X3Y2Z1) return pes_energy_bohr_x3y2z1(coords.data());
		std::cout<<"ERROR: You forgot to add your new system to pes_interface_cppside.hpp"<<std::endl;
		abort();
	}else{
		if constexpr (PES_SYSTEM_TYPE == PES_X3Y1Z1U1) return pes_energy_x3y1z1u1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X4Y1Z1U1) return pes_energy_x4y1z1u1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X4Y2Z1) return pes_energy_x4y2z1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X5Y2Z1) return pes_energy_x5y2z1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X3Y2Z1) return pes_energy_x3y2z1(coords.data());
		std::cout<<"ERROR: You forgot to add your new system to pes_interface_cppside.hpp"<<std::endl;
		abort();
	}
}

template <bool geomInBohrs, auto NUM_COORDS>
double EvaluateEnergy_copy(std::array<double,NUM_COORDS> coords){
	///return the potential enegy at a point specified by a 1D (linearized) std::array of Cartesian coordinates
	///this version gets the coordinates by value, instead of reference, to allow more flexible usage
	static_assert(IsSystemBraamsSupported(), "System type is not supported by Braams library.");
	if constexpr (geomInBohrs){
		if constexpr (PES_SYSTEM_TYPE == PES_X3Y1Z1U1) return pes_energy_bohr_x3y1z1u1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X4Y1Z1U1) return pes_energy_bohr_x4y1z1u1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X4Y2Z1) return pes_energy_bohr_x4y2z1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X5Y2Z1) return pes_energy_bohr_x5y2z1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X3Y2Z1) return pes_energy_bohr_x3y2z1(coords.data());
		std::cout<<"ERROR: You forgot to add your new system to pes_interface_cppside.hpp"<<std::endl;
		abort();
	}else{
		if constexpr (PES_SYSTEM_TYPE == PES_X3Y1Z1U1) return pes_energy_x3y1z1u1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X4Y1Z1U1) return pes_energy_x4y1z1u1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X4Y2Z1) return pes_energy_x4y2z1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X5Y2Z1) return pes_energy_x5y2z1(coords.data());
		if constexpr (PES_SYSTEM_TYPE == PES_X3Y2Z1) return pes_energy_x3y2z1(coords.data());
		std::cout<<"ERROR: You forgot to add your new system to pes_interface_cppside.hpp"<<std::endl;
		abort();
	}
}
