#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <cassert>
#include <quadmath.h>

#include "../../common_headers/misc_util.hpp"
#include "../../common_headers/storageClasses.hpp"
#include "../../common_headers/XYZ_IO.hpp"
#include "../../common_headers/atomic_masses.hpp"
#include "../../common_headers/system_types.hpp"

struct VECT{
	__float128 x,y,z;
};

struct QUAT{
	__float128 r, i, j, k;
};

struct ROTMAT{
	__float128 mat[9];
};

QUAT CalcRotQuat(const __float128& X0, const __float128& X1, const __float128& X2){
	///based on Ken Shoemake's "UNIFORM RANDOM ROTATIONS" book chapter in "Graphics Gems III"
	///DOI: 10.1016/B978-0-08-050755-2.50036-1
	const __float128 r1 = sqrtq(1.0q-X0);
	const __float128 r2 = sqrtq(X0);
	const __float128 X1_interm = 2.0q*M_PIq*X1;
	const __float128 X2_interm = 2.0q*M_PIq*X2;
	const __float128 s1 = sinq(X1_interm);
	const __float128 c1 = cosq(X1_interm);
	const __float128 s2 = sinq(X2_interm);
	const __float128 c2 = cosq(X2_interm);
	QUAT rotQ;
	rotQ.r = s1*r1;
	rotQ.i = c1*r1;
	rotQ.j = s2*r2;
	rotQ.k = c2*r2;
	return rotQ;
}

ROTMAT Quat2Rotmat(const QUAT& Q){
	ROTMAT R;
	R.mat[0] = 1.0q - 2.0q*Q.j*Q.j - 2.0q*Q.k*Q.k;
	R.mat[1] = 2.0q*Q.i*Q.j - 2.0q*Q.r*Q.k;
	R.mat[2] = 2.0q*Q.i*Q.k + 2.0q*Q.r*Q.j;
	R.mat[3] = 2.0q*Q.i*Q.j + 2.0q*Q.r*Q.k;
	R.mat[4] = 1.0q - 2.0q*Q.i*Q.i - 2.0q*Q.k*Q.k;
	R.mat[5] = 2.0q*Q.j*Q.k - 2.0q*Q.r*Q.i;
	R.mat[6] = 2.0q*Q.i*Q.k - 2.0q*Q.r*Q.j;
	R.mat[7] = 2.0q*Q.j*Q.k + 2.0q*Q.r*Q.i;
	R.mat[8] = 1.0q - 2.0q*Q.i*Q.i - 2.0q*Q.j*Q.j;
	return R;
}

ROTMAT GenRandRotmat(auto& PRNG, std::uniform_real_distribution<double>& URD){
	const __float128 X0 = URD(PRNG);
	const __float128 X1 = URD(PRNG);
	const __float128 X2 = URD(PRNG);
	const QUAT rotQ = CalcRotQuat(X0, X1, X2); //calculate unit quaternion of rotation (rotQ)
	return Quat2Rotmat(rotQ); //convert the unit quaternion to a rotation matrix
}

void RotateAtom(atom& input, const ROTMAT& R){
	const __float128 x = input.xcoord;
	const __float128 y = input.ycoord;
	const __float128 z = input.zcoord;
	const __float128 x_out = R.mat[0]*x + R.mat[1]*y + R.mat[2]*z;
	const __float128 y_out = R.mat[3]*x + R.mat[4]*y + R.mat[5]*z;
	const __float128 z_out = R.mat[6]*x + R.mat[7]*y + R.mat[8]*z;
	input.xcoord = x_out;
	input.ycoord = y_out;
	input.zcoord = z_out;
}

template<char fragSel='A', auto frag2size=0>
__float128 CalcFragMass(const xyzfile& geom, const std::array<uint32_t, frag2size>* const frag2_idx=nullptr){
	if constexpr(fragSel != 'A') assert(frag2_idx != nullptr);
	static_assert( (fragSel=='A') || (fragSel=='F') || (fragSel=='S') );
	
	__float128 fragMass = 0.0q;
	for(uint32_t j=0; j<geom.natoms; j++){
		if constexpr(fragSel == 'A'){
			fragMass += getAtomMass(geom.atomvec[j]);
		}else{
			if constexpr(fragSel == 'F'){
				if ( std::none_of(frag2_idx->cbegin(), frag2_idx->cend(), [&](uint32_t idx){ return idx == j; }) ){
					fragMass += getAtomMass(geom.atomvec[j]);
				}
			}else{
				if ( std::any_of(frag2_idx->cbegin(), frag2_idx->cend(), [&](uint32_t idx){ return idx == j; }) ){
					fragMass += getAtomMass(geom.atomvec[j]);
				}
			}
		}
	}
	return fragMass;
}

template<char fragSel='A', auto frag2size=0>
VECT CalcCOM(const xyzfile& geom, const std::array<uint32_t, frag2size>* const frag2_idx=nullptr){
	if constexpr(fragSel != 'A') assert(frag2_idx != nullptr);
	if (frag2_idx != nullptr) assert(fragSel != 'A');
	static_assert( (fragSel=='A') || (fragSel=='F') || (fragSel=='S') );
	
	const __float128 fragMass = CalcFragMass<fragSel>(geom, frag2_idx);
	//std::cout<<"mass "<< static_cast<double>(fragMass)<< std::endl;
	VECT COM;
	COM.x = 0.0q;
	COM.y = 0.0q;
	COM.z = 0.0q;
	for(uint32_t j=0; j<geom.natoms; j++){
		if constexpr(fragSel == 'A'){
			const __float128 mass = getAtomMass(geom.atomvec[j]);
			COM.x += mass*geom.atomvec[j].xcoord;
			COM.y += mass*geom.atomvec[j].ycoord;
			COM.z += mass*geom.atomvec[j].zcoord;
		}else{
			if constexpr(fragSel == 'F'){
				if ( std::none_of(frag2_idx->cbegin(), frag2_idx->cend(), [&](uint32_t idx){ return idx == j; }) ){
					const __float128 mass = getAtomMass(geom.atomvec[j]);
					COM.x += mass*geom.atomvec[j].xcoord;
					COM.y += mass*geom.atomvec[j].ycoord;
					COM.z += mass*geom.atomvec[j].zcoord;
				}
			}else{
				if ( std::any_of(frag2_idx->cbegin(), frag2_idx->cend(), [&](uint32_t idx){ return idx == j; }) ){
					const __float128 mass = getAtomMass(geom.atomvec[j]);
					COM.x += mass*geom.atomvec[j].xcoord;
					COM.y += mass*geom.atomvec[j].ycoord;
					COM.z += mass*geom.atomvec[j].zcoord;
				}
			}
		}
	}
	COM.x /= fragMass;
	COM.y /= fragMass;
	COM.z /= fragMass;
	//std::cout<<"COM "<< static_cast<double>(COM.x)<<"  "<< static_cast<double>(COM.y)<<"  "<< static_cast<double>(COM.z)<<"  "<< std::endl;
	return COM;
}

void TranslateOriginToCoords(xyzfile& geom, const VECT& coords){
	///translates the origin of the Cartesian coordinate system to a given point
	for(uint32_t i=0; i<geom.natoms; i++){
		geom.atomvec[i].xcoord -= coords.x;
		geom.atomvec[i].ycoord -= coords.y;
		geom.atomvec[i].zcoord -= coords.z;
	}
}

VECT Scale(const VECT& a, const __float128 factor){
	VECT c;
	c.x = factor * a.x;
	c.y = factor * a.y;
	c.z = factor * a.z;
	return c;
}

__float128 CalcNorm(const VECT& a){
	return sqrtq(a.x*a.x + a.y*a.y + a.z*a.z);
}

VECT Normalize(const VECT& a){
	const __float128 norm = CalcNorm(a);
	return Scale(a, 1.0q/norm);
}

template<char fragSel='A', auto frag2size=0>
void TranslateFrag(xyzfile& geom, const VECT& transVec, const std::array<uint32_t, frag2size>* const frag2_idx=nullptr){
	if constexpr(fragSel != 'A') assert(frag2_idx != nullptr);
	if (frag2_idx != nullptr) assert(fragSel != 'A');
	static_assert( (fragSel=='A') || (fragSel=='F') || (fragSel=='S') );
	
	for(uint32_t j=0; j<geom.natoms; j++){
		if constexpr(fragSel == 'A'){
			geom.atomvec[j].xcoord += transVec.x;
			geom.atomvec[j].ycoord += transVec.y;
			geom.atomvec[j].zcoord += transVec.z;
		}else{
			if constexpr(fragSel == 'F'){
				if ( std::none_of(frag2_idx->cbegin(), frag2_idx->cend(), [&](uint32_t idx){ return idx == j; }) ){
					geom.atomvec[j].xcoord += transVec.x;
					geom.atomvec[j].ycoord += transVec.y;
					geom.atomvec[j].zcoord += transVec.z;
				}
			}else{
				if ( std::any_of(frag2_idx->cbegin(), frag2_idx->cend(), [&](uint32_t idx){ return idx == j; }) ){
					geom.atomvec[j].xcoord += transVec.x;
					geom.atomvec[j].ycoord += transVec.y;
					geom.atomvec[j].zcoord += transVec.z;
				}
			}
		}
	}
}

template<typename T>
double AtomDist(const T& input, const uint32_t geomidx, const uint32_t atom1, const uint32_t atom2){
	/// WIP
	///
	//calculate the distance between two atoms
	const uint32_t NATOMS = GetNatoms();
    return std::sqrt(std::pow(input.xcoordVec[geomidx*NATOMS+atom2]-input.xcoordVec[geomidx*NATOMS+atom1],2) + std::pow(input.ycoordVec[geomidx*NATOMS+atom2]-input.ycoordVec[geomidx*NATOMS+atom1],2) + std::pow(input.zcoordVec[geomidx*NATOMS+atom2]-input.zcoordVec[geomidx*NATOMS+atom1],2));
}

template<typename T>
void CalcDistmat(double* const dm, const T& input, const uint32_t geomidx){
/// WIP
///
	const uint32_t NATOMS = GetNatoms();
	uint32_t ctr = 0;
    for (uint32_t i=1; i<NATOMS; i++){
        for (uint32_t j=0; j<i; j++){
            dm[ctr] = std::exp(AtomDist(input, geomidx, i, j)/-6.0);
            ctr++;
        }
    }
}

template<uint32_t distsPerGeom>
double CompareDistmats(const double*__restrict__ const distmat1, const double*__restrict__ const distmat2){
	///WIP
	//
	//compare two modified distance matrices and calculate their RMS difference
    //this variant avoids the nested loop by operating on a flattened matrix
    //should be faster due to less loop overhead, better vectorization and improved CPU cache use due to sequential access
    double RMSD = 0.0;
    for (uint32_t i=0; i<distsPerGeom; i++){
			RMSD += std::pow(distmat1[i]-distmat2[i],2);
    }
	return std::sqrt(RMSD / distsPerGeom);
}

template <typename T>
void PruneGems(T& input, const double RMSD_TRESH){
	const uint32_t natoms = GetNatoms();
	const uint32_t distsPerGeom = (natoms*(natoms-1))/2;
	assert(natoms == input.natoms);
	
	T prunedInput;
	prunedInput.reserve(input.size(), natoms);
	prunedInput.natoms = input.natoms;
	double tmpDistmat[distsPerGeom];
	std::vector<double> prunedDistmatVec;
	prunedDistmatVec.reserve(static_cast<uint64_t>(input.size()) * distsPerGeom); //to prevent int overflow when this array needs to be very big (>32 GiB)
	double RMSD;
	
	const uint32_t inputsize = input.size();
	const uint32_t progressReportInterval = std::max<uint32_t>(inputsize/10, 1);
	for(uint32_t i = 0; i<inputsize; i++){
		bool dissimilar=true;
		CalcDistmat(tmpDistmat, input, i);
		const uint32_t currPrunedSize = prunedInput.energyVec.size();
		for(uint32_t j = 0; j<currPrunedSize; j++){
			RMSD = CompareDistmats<distsPerGeom>(tmpDistmat, &prunedDistmatVec[distsPerGeom*j]);
			if (RMSD < RMSD_TRESH){
				dissimilar=false;
				break;
			}
		}
		if (dissimilar){
			prunedInput.energyVec.push_back(input.energyVec[i]);
			prunedInput.extralineVec.emplace_back(input.extralineVec[i]); //this should be just a copy constructor call
			prunedInput.symbolVec.insert(prunedInput.symbolVec.end(), input.symbolVec.cbegin()+(i*natoms*3), input.symbolVec.cbegin()+((i+1)*natoms*3));
			prunedInput.xcoordVec.insert(prunedInput.xcoordVec.end(), input.xcoordVec.cbegin()+(i*natoms), input.xcoordVec.cbegin()+((i+1)*natoms));
			prunedInput.ycoordVec.insert(prunedInput.ycoordVec.end(), input.ycoordVec.cbegin()+(i*natoms), input.ycoordVec.cbegin()+((i+1)*natoms));
			prunedInput.zcoordVec.insert(prunedInput.zcoordVec.end(), input.zcoordVec.cbegin()+(i*natoms), input.zcoordVec.cbegin()+((i+1)*natoms));
			if constexpr ((std::is_same<T,multiGEM>::value) || (std::is_same<T,multiDOS>::value)){
				prunedInput.trajNumVec.push_back(input.trajNumVec[i]);
				prunedInput.trajStepVec.push_back(input.trajStepVec[i]);
			}
			if constexpr (std::is_same<T,multiDOS>::value){
				prunedInput.xvelVec.insert(prunedInput.xvelVec.end(), input.xvelVec.cbegin()+(i*natoms), input.xvelVec.cbegin()+((i+1)*natoms));
				prunedInput.yvelVec.insert(prunedInput.yvelVec.end(), input.yvelVec.cbegin()+(i*natoms), input.yvelVec.cbegin()+((i+1)*natoms));
				prunedInput.zvelVec.insert(prunedInput.zvelVec.end(), input.zvelVec.cbegin()+(i*natoms), input.zvelVec.cbegin()+((i+1)*natoms));
			}
			prunedDistmatVec.insert(prunedDistmatVec.end(), tmpDistmat, tmpDistmat+distsPerGeom);
		}
		if ((i%progressReportInterval) == 0) std::cout<<i/progressReportInterval<<"0% of input processed"<<std::endl;
	}
	input = prunedInput;
	input.shrink_to_fit();
}

template <uint32_t NUM_ATOMS>
void randPertGeom(xyzfile& geom_pert, std::mt19937_64& PRNG, std::uniform_real_distribution<double>& RN){
	///Given an initalized pseudorandom number generator engine, and a desired uniform distribution, randomly perturb the coordinates
	///of the given geometry.
	for(uint32_t i=0; i<NUM_ATOMS; i++){
		geom_pert.atomvec[i].xcoord += RN(PRNG);
		geom_pert.atomvec[i].ycoord += RN(PRNG);
		geom_pert.atomvec[i].zcoord += RN(PRNG);
	}
}


int main(){
	std::array<uint32_t,1> frag2_idx = {7}; //atoms of the smaller reactant (count starts at 0)
	constexpr uint32_t NATOMS = 8; //total number of atoms
	constexpr __float128 APPROACH_STEP = 0.02; //translation distance per approach step
	constexpr uint32_t ROTATIONS_PER_APPROACH_STEP = 10000; //number of random rotations per approach step
	constexpr __float128 APPROACH_MIN_COMDIST = 1.0; //stop generating approach steps if the two centers of mass are closer than this
	constexpr double PRUNING_THRESH = 0.0175; //EW-RMSD pruning threshold for pruning too similar geometries
	const std::string REF_GEOM_PATH = "ref_geom_menh2_clrad_rmp2_avdz_opt_qsd_manyhf_v13.xyz"; //reference geometry with optimized free reactants
	
	xyzfile ref_geom = LoadSingleXYZ<xyzfile>(REF_GEOM_PATH, 123456.7, -123456.7, 1.0E+6);
	assert(ref_geom.natoms == NATOMS);
	{
		std::vector<xyzfile> rotatedGeomVec;
		auto PRNG = initPRNG();
		std::uniform_real_distribution<double> URD(0.0, 1.0);
		std::uniform_real_distribution<double> pertDist(-0.1, 0.1);
		while(true){
			for(uint32_t i=0; i<ROTATIONS_PER_APPROACH_STEP; i++){
				xyzfile rotated_geom = ref_geom;
				randPertGeom<NATOMS>(rotated_geom, PRNG, pertDist);
				{ //rotate the smaller reactant
					const ROTMAT R = GenRandRotmat(PRNG, URD);
					TranslateOriginToCoords(rotated_geom, CalcCOM<'S'>(rotated_geom, &frag2_idx));
					for(uint32_t j=0; j<ref_geom.natoms; j++){
						if (std::any_of(frag2_idx.cbegin(), frag2_idx.cend(), [&](uint32_t idx){ return idx == j;})){
							RotateAtom(rotated_geom.atomvec[j], R);
						}
					}
				}
				{ //rotate the larger reactant
					const ROTMAT R = GenRandRotmat(PRNG, URD);
					TranslateOriginToCoords(rotated_geom, CalcCOM<'F'>(rotated_geom, &frag2_idx));
					for(uint32_t j=0; j<ref_geom.natoms; j++){
						if (std::none_of(frag2_idx.cbegin(), frag2_idx.cend(), [&](uint32_t idx){ return idx == j;})){
							RotateAtom(rotated_geom.atomvec[j], R);
						}
					}
				}
				TranslateOriginToCoords(rotated_geom, CalcCOM(rotated_geom));
				rotatedGeomVec.push_back(rotated_geom);
			}
			
			TranslateOriginToCoords(ref_geom, CalcCOM<'F'>(ref_geom, &frag2_idx));
			const VECT frag2vec = CalcCOM<'S'>(ref_geom, &frag2_idx);
			if (const __float128 norm = CalcNorm(frag2vec); norm < APPROACH_MIN_COMDIST){
				break;
			}else{
				std::cout<<static_cast<double>(norm)<<std::endl;
				TranslateFrag<'S'>(ref_geom, Scale(frag2vec, -APPROACH_STEP/norm), &frag2_idx);
			}
		};
		WriteMultiXYZ("approach.xyz",rotatedGeomVec);
	}
	
	std::vector<xyzfile> pruned_geoms;
	{ //prune the generated geometries
		multiXYZ gen_geoms;
		assert(0==LoadMultiXYZ(gen_geoms, "approach.xyz", 123456.7, -123456.7, 1.0E+6));
		const auto PruningStart = std::chrono::steady_clock::now();
		PruneGems(gen_geoms, PRUNING_THRESH);
		const auto PruningFinish = std::chrono::steady_clock::now();
		std::cout << "Done in "<< std::chrono::duration <double> (PruningFinish - PruningStart).count() << " seconds" << std::endl;
		assert(CheckMultiXYZ_InternalConsistence(gen_geoms) == 0);
		gen_geoms.isValidated = true;
		std::cout<<gen_geoms.size()<<" geometries remain after pruning"<<std::endl<<std::endl;
		
		pruned_geoms.reserve(gen_geoms.size());
		for(uint32_t i=0; i<gen_geoms.size(); i++) pruned_geoms.push_back(gen_geoms.getXYZ(i));
	}
	WriteMultiXYZ("approach_pruned.xyz",pruned_geoms);
	
	return 0;
}
