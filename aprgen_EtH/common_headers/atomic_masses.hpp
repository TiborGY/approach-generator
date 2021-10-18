constexpr double H_mass  = 1.0078250322;    // hydrogen-1        http://www.ciaaw.org/atomic-masses.htm
constexpr double C_mass  =         12.0;    // carbon-12
constexpr double N_mass  = 14.003074004;    // nitrogen-14       http://www.ciaaw.org/atomic-masses.htm
constexpr double O_mass  = 15.994914619;    // oxygen-16         http://www.ciaaw.org/atomic-masses.htm
constexpr double F_mass  = 18.998403163;    // fluorine-19       http://www.ciaaw.org/atomic-masses.htm
constexpr double Cl_mass =      35.4515;    // isotope averaged, midpoint of the [35.446, 35.457] range  http://www.ciaaw.org/atomic-weights.htm
constexpr double Br_mass =       79.904;    // isotope averaged, midpoint of the [79.901, 79.907] range  http://www.ciaaw.org/atomic-weights.htm
constexpr double I_mass  =    126.90447;    // iodine-127        http://www.ciaaw.org/atomic-masses.htm

double getAtomMass(const atom& atm){
	if(std::string(atm.symbol) == "H")  return H_mass;
	if(std::string(atm.symbol) == "C")  return C_mass;
	if(std::string(atm.symbol) == "N")  return N_mass;
	if(std::string(atm.symbol) == "O")  return O_mass;
	if(std::string(atm.symbol) == "F")  return F_mass;
	if(std::string(atm.symbol) == "Cl") return Cl_mass;
	if(std::string(atm.symbol) == "Br") return Br_mass;
	if(std::string(atm.symbol) == "I")  return I_mass;
	std::cout<<"ERROR: atom symbol "<<std::string(atm.symbol)<<" has no mass assigned!"<<std::endl;
	abort();
}

double getFragmentMass(const std::vector<atom>& atomvec){
	double mass = 0.0;
	for(const auto& a : atomvec) mass += getAtomMass(a);
	return mass;
}

template <typename T>
double getFragmentMass(const T& geom){
	return getFragmentMass(geom.atomvec);
}
