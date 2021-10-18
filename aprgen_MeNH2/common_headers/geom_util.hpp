#include "cblas.h"
#include "lapacke.h"
#include "physical_constants.hpp"
#include "atomic_masses.hpp"
#include "pes_interface_cppside.hpp"

template <auto NUM_COORDS>
void randPertGeom(std::array<double,NUM_COORDS>& geom_pert, std::mt19937_64& PRNG, std::uniform_real_distribution<double>& RN){
	///Given an initalized pseudorandom number generator engine, and a desired uniform distribution, randomly perturb the coordinates
	///of the given geometry.
	for(auto& coord : geom_pert) coord += RN(PRNG);
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

template <bool geomInBohrs, auto NUM_COORDS>
double calcDisplacedE(std::array<double,NUM_COORDS> geom, const uint32_t idx, const double disp){
	///Calculate the energy at a point with one coordinate displaced.
	geom[idx] += disp;
	return EvaluateEnergy<geomInBohrs>(geom);
}

template <bool geomInBohrs, auto NUM_COORDS>
double calcDisplacedE(std::array<double,NUM_COORDS> geom, const uint32_t idx_1, const uint32_t idx_2, const double disp_1, const double disp_2){
	///Calculate the energy at a point with two coordinates displaced.
	geom[idx_1] += disp_1;
	geom[idx_2] += disp_2;
	return EvaluateEnergy<geomInBohrs>(geom);
}

template <bool extraAcc, bool geomInBohrs, auto NUM_COORDS>
void calcGrad(const std::array<double,NUM_COORDS>& geom, std::array<double,NUM_COORDS>& grad, const double disp){
	///Calculate gradient vector via two-sided numerical differentiation.
	///If requested, a more accurate six-point numerical differentiation formula will be used.
	[[maybe_unused]] const double inv_2h = static_cast<double>(1)/(2*disp);
	[[maybe_unused]] const double inv_12h = static_cast<double>(1)/(12*disp);
	[[maybe_unused]] const double inv_60h = static_cast<double>(1)/(60*disp);
	for (uint32_t i=0; i<grad.size(); i++){
		if constexpr (!extraAcc){
			grad[i] = inv_2h*( calcDisplacedE<geomInBohrs>(geom, i, disp) - calcDisplacedE<geomInBohrs>(geom, i, -disp) );
		}else{
			//grad[i] = inv_12h*( calcDisplacedE<geomInBohrs>(geom, i, -2*disp) - 8*calcDisplacedE<geomInBohrs>(geom, i, -disp) + 8*calcDisplacedE<geomInBohrs>(geom, i, disp) - calcDisplacedE<geomInBohrs>(geom, i, 2*disp));
			grad[i] = inv_60h*(
				- calcDisplacedE<geomInBohrs>(geom, i, -3*disp) + 9*calcDisplacedE<geomInBohrs>(geom, i, -2*disp)
				- 45*calcDisplacedE<geomInBohrs>(geom, i, -disp) + 45*calcDisplacedE<geomInBohrs>(geom, i, disp)
				- 9*calcDisplacedE<geomInBohrs>(geom, i, 2*disp) + calcDisplacedE<geomInBohrs>(geom, i, 3*disp)
			);
		}
	}
}
template <bool extraAcc, bool geomInBohrs, auto NUM_COORDS>
std::array<double,NUM_COORDS> calcGrad(const std::array<double,NUM_COORDS>& geom, const double disp){
	///Calculate gradient vector via two-sided numerical differentiation.
	///If requested, a more accurate six-point numerical differentiation formula will be used.
	std::array<double,NUM_COORDS> grad;
	calcGrad<extraAcc, geomInBohrs>(geom, grad, disp);
	return grad;
}

template <bool extraAcc, bool geomInBohrs, auto NUM_COORDS>
void calcHess(const std::array<double,NUM_COORDS>& geom, std::array<double,NUM_COORDS*NUM_COORDS>& hess, const double disp){
	///Calculate Hessian matrix via two-sided numerical differentiation.
	///If requested, more accurate (fourth order) formulas are used. (5 point for diagonal, 8 point for off-diagonal elements)
	///While the Hessian is exactly symmetric in infinite precision, due to finite precision and especially due to numerical differentiation
	///this would not be the case by default. Regardless, the Hessian will be populated to be exactly symmetric, by using (i,j) for (i,j) and (j,i).
	[[maybe_unused]] const double twoEcenter = 2*EvaluateEnergy_copy<geomInBohrs>(geom);
	[[maybe_unused]] const double thirtyEcenter = 30*EvaluateEnergy_copy<geomInBohrs>(geom);
	[[maybe_unused]] const double inv_4sqrh = static_cast<double>(1)/(4*disp*disp);
	[[maybe_unused]] const double inv_12sqrh = static_cast<double>(1)/(12*disp*disp);
	[[maybe_unused]] const double inv_48sqrh = static_cast<double>(1)/(48*disp*disp);
	[[maybe_unused]] const double inv_sqrh = static_cast<double>(1)/(disp*disp);
	//does not touch the diagonal
	for (uint32_t i=1; i<geom.size(); i++){
		for (uint32_t j=0; j<i; j++){
			if constexpr (!extraAcc){
				hess[i*geom.size() + j] = inv_4sqrh*(
					calcDisplacedE<geomInBohrs>(geom, i, j, disp, disp)
					- calcDisplacedE<geomInBohrs>(geom, i, j, disp, -disp)
					- calcDisplacedE<geomInBohrs>(geom, i, j, -disp, disp)
					+ calcDisplacedE<geomInBohrs>(geom, i, j, -disp, -disp) );
			}else{
				hess[i*geom.size() + j] = inv_48sqrh*(
					  - calcDisplacedE<geomInBohrs>(geom, i, j, +2*disp, +2*disp) + 16.0*calcDisplacedE<geomInBohrs>(geom, i, j, +disp, +disp)
					  + calcDisplacedE<geomInBohrs>(geom, i, j, +2*disp, -2*disp) - 16.0*calcDisplacedE<geomInBohrs>(geom, i, j, +disp, -disp)
					  + calcDisplacedE<geomInBohrs>(geom, i, j, -2*disp, +2*disp) - 16.0*calcDisplacedE<geomInBohrs>(geom, i, j, -disp, +disp)
					  - calcDisplacedE<geomInBohrs>(geom, i, j, -2*disp, -2*disp) + 16.0*calcDisplacedE<geomInBohrs>(geom, i, j, -disp, -disp)
				);
			}
			hess[j*geom.size() + i] = hess[i*geom.size() + j];
		}
	}
	//diagonal elements
	for (uint32_t i=0; i<geom.size(); i++){
		//std::cout<<"diag "<<i<<" value "<<hess[i*geom.size() + i]<<std::endl;
		if constexpr (!extraAcc){
			hess[i*geom.size() + i] = inv_sqrh*( calcDisplacedE<geomInBohrs>(geom, i, disp) - twoEcenter + calcDisplacedE<geomInBohrs>(geom, i, -disp) );
		}else{
			hess[i*geom.size() + i] = inv_12sqrh*(
				- calcDisplacedE<geomInBohrs>(geom, i, 2*disp)
				+ 16.0*calcDisplacedE<geomInBohrs>(geom, i, disp)
				- thirtyEcenter
				+ 16.0*calcDisplacedE<geomInBohrs>(geom, i, -disp)
				- calcDisplacedE<geomInBohrs>(geom, i, -2*disp)
			);
		}
	}
}

template <auto NUM_COORDS>
double calcRMS(const std::array<double,NUM_COORDS>& arr){
	///Calculates the Root-Mean-Square of the elements in an array.
	double RMS = 0.0;
	for(uint32_t i=0; i<arr.size(); i++) RMS += arr[i] * arr[i];
	return std::sqrt(RMS/NUM_COORDS);
}

template <auto NUM_COORDS>
double calcMaxAbs(const std::array<double,NUM_COORDS>& arr){
	///Calculates the maximum absolute value of the elements in an array.
	double max = 0.0;
	for(uint32_t i=0; i<arr.size(); i++) if (std::abs(arr[i]) > max) max = std::abs(arr[i]);
	return max;
}

template <auto NUM_COORDS>
std::array<double,NUM_COORDS> arrayAdd(const std::array<double, NUM_COORDS>& input1, const std::array<double, NUM_COORDS>& input2){
	std::array<double, NUM_COORDS> tmp;
	for (uint32_t j=0; j<NUM_COORDS; j++) tmp[j] = input1[j]+input2[j];
	return tmp;
}

template <auto NUM_COORDS>
std::array<double,NUM_COORDS> arrayScale(const std::array<double, NUM_COORDS>& input, const double scalar){
	std::array<double, NUM_COORDS> tmp;
	for (uint32_t j=0; j<NUM_COORDS; j++) tmp[j] = input[j]*scalar;
	return tmp;
}

template <auto NUM_COORDS>
[[nodiscard]] bool diagMat(std::array<double,NUM_COORDS*NUM_COORDS>& mat, std::array<double,NUM_COORDS>& eigenvals){
	///Diagonalizes square and symmetric matrix using the DSYEV LAPACK subroutine, yielding eigenvalues and eigenvectors.
	///The input matrix is overwritten with the eigenvectors, if you need them make sure to mind the Fortran-style column-major indexing.
	///The caller must check the return value of this fn, as the diagonalizer can fail, yielding garbage results.
	const lapack_int N = eigenvals.size();
	eigenvals.fill(0.0);
	if (0 != LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U', N, mat.data(), N, eigenvals.data())){
		std::cout<<"Could not diagonalize matrix, DSYEV failed!"<<std::endl;
		return false;
	}
	return true; //eigenvals should now be filled with eigenvalues in ascending order
}

template <bool silent, auto NUM_COORDS>
[[nodiscard]] bool calcNewtonStep(std::array<double,NUM_COORDS*NUM_COORDS> hess, const std::array<double,NUM_COORDS>& grad_in, std::array<double,NUM_COORDS>& step, const double RCOND){
	///This function calculates the undamped Newton optimization step (i.e. displacement vector) from the gradient and Hessian via SVD.
	///Takes the Hessian by value as it will be overwritten with the SVD details.
	///Takes the gradient by const ref. and make a copy, the copy will be negated and then used in the SVD.
	///If the SVD succeeds, copy the resulting vector into step, otherwise step is not modified.
	///The caller must check the return value, to make sure the SVD did not fail.
	auto grad = grad_in;
	std::array<double,grad.size()> S;
	const lapack_int N = grad.size();
	lapack_int rank = 0;
	S.fill(0.0);
	for(auto& der : grad) der = -der; //turn grad into -grad
	//find a minimum norm solution to the system of linear eqs. H d = -g, where d is the optimization step
	//this uses SVD, as the Hessian should in theory be non-invertible (6 zero eigenvalues, corresponding to translations and rotations)
	//since the potential energy is invariant to those, the second derivative with respect to them is zero
	//in practice tiny eigenvalues creep in, so the matrix may have numerically fully rank
	//regardless, the SVD should be able to handle all cases
	//numerical experiments suggest RCOND should be set to around 1.0E-5 to 1.0E-3
	if (0 != LAPACKE_dgelss(LAPACK_COL_MAJOR, N, N, 1, hess.data(), N, grad.data(), N, S.data(), RCOND, &rank) ){
		std::cout<<"DGELSS failed!"<<std::endl;
		return false;
	}
	if constexpr (!silent) std::cout<<"Hessian rank: "<<rank<<std::endl;
	step = grad;
	return true;
}

template <auto NUM_ATOMS>
void massWeightHess(const std::array<double, NUM_ATOMS>& masses, std::array<double,(NUM_ATOMS*3)*(NUM_ATOMS*3)>& hess){
	const auto NUM_COORDS = NUM_ATOMS*3;
	for (uint32_t i=0; i<NUM_COORDS; i++){
		for (uint32_t j=0; j<NUM_COORDS; j++){
			hess.at(i*NUM_COORDS + j) /= std::sqrt( emass_per_amu * masses.at(i/3) * emass_per_amu * masses.at(j/3) );
		}
	}
}

template <auto NUM_COORDS>
void calcFreq_core(const std::array<double, NUM_COORDS/3>& masses, std::array<double,NUM_COORDS*NUM_COORDS>& hessian, std::array<double, NUM_COORDS>& freqs){
	std::array<double,NUM_COORDS> eigenvals;
	massWeightHess(masses, hessian);
	assert(diagMat(hessian, eigenvals));
	for (uint32_t i=0; i<NUM_COORDS; i++){
		freqs.at(i) = std::sqrt(std::abs(eigenvals.at(i)))*invcm_per_hartree;
		if (eigenvals.at(i) < 0.0) freqs.at(i) *= -1;
	}
}

template <auto NUM_COORDS>
std::array<double,NUM_COORDS*NUM_COORDS> calcFreq(const xyzfile& geom_in, std::array<double, NUM_COORDS>& freqs, const double numderDisp){
	/// Eigenvectors are returned in Fortran-style column-major indexing!
	assert(geom_in.natoms == (NUM_COORDS/3));
	std::array<double,NUM_COORDS/3> masses;
	std::array<double,NUM_COORDS> geom;
	std::array<double,NUM_COORDS*NUM_COORDS> hessian;
	
	ExtractLinearizedCoords(geom_in, geom);
	// Convert atom coordinates to bohrs get the masses of the atoms.
	for (uint32_t j=0; j<NUM_COORDS; j++) geom.at(j) /= angs_per_bohr;
	for (uint32_t j=0; j<geom_in.natoms; j++) masses.at(j) = getAtomMass(geom_in.atomvec.at(j));
	// Calculate the Hessian, use the more accurate derivative formulas.
	calcHess<true, true>(geom, hessian, numderDisp/angs_per_bohr);
	// Calculate the frequencies and the eigenvectors of the mass-weighted Hessian.
	calcFreq_core(masses, hessian, freqs);
	return hessian;
}

template <auto NUM_COORDS>
void matmul(std::array<double,NUM_COORDS*NUM_COORDS> A, std::array<double,NUM_COORDS*NUM_COORDS> B,
  std::array<double,NUM_COORDS*NUM_COORDS>& C, const bool transposeA, const bool transposeB){
	enum CBLAS_TRANSPOSE TransA;
	enum CBLAS_TRANSPOSE TransB;
	if(transposeA){
		TransA = CblasTrans;
	}else{
		TransA = CblasNoTrans;
	}
	if(transposeB){
		TransB = CblasTrans;
	}else{
		TransB = CblasNoTrans;
	}
	const blasint N = NUM_COORDS;
	
	cblas_dgemm(CblasColMajor, TransA, TransB, N, N, N, 1.0, A.data(), N, B.data(), N, 0.0, C.data(), N);
}

template <auto NUM_COORDS>
void matvecmul(std::array<double,NUM_COORDS*NUM_COORDS> A, std::array<double,NUM_COORDS> x,
  std::array<double,NUM_COORDS>& y){
	const blasint N = NUM_COORDS;
	const blasint inc = 1;
	cblas_dgemv(CblasColMajor, CblasNoTrans, N, N, 1.0, A.data(), N, x.data(), inc, 0.0, y.data(), inc);
}

template <auto NUM_ATOMS>
std::array<double,(NUM_ATOMS*3)*(NUM_ATOMS*3)> invertHessian(std::array<double,(NUM_ATOMS*3)*(NUM_ATOMS*3)> hess_in,
  const std::array<double,NUM_ATOMS>& masses, const uint32_t numImagKeep){
	/// Calculate a modified inverse of the Hessian via eigendecomposition.
	/// H^-1 = Q * L^-1 * Q^-1 , where Q is the matrix where the columns are eigenvectors and L is a digonal mtx,
	/// where the diagonal elements are the eigenvalues.
	/// Due to the near-singular nature of the Hessian in Cartesian coordinates, six of the eigenvalues should be
	/// close to zero. These are explicitly set to zero
	const auto NUM_COORDS = NUM_ATOMS*3;
	const double smallFreqThresh = 15.0; //in cm-1
	std::array<double,NUM_COORDS> eigenvals, freqs;
	std::array<double,NUM_COORDS*NUM_COORDS> diagInverse;
	// Calculate harmonic vibrational frequencies to be able to decide how many eigenvalues to null out when
	// constructing the Hessian inverse. We need a copy of the Hessian to keep the original intact.
	auto hess_mw = arrayScale(hess_in, (angs_per_bohr)*(angs_per_bohr));
	calcFreq_core(masses, hess_mw, freqs);
	uint32_t numSmallFreq = 0;
	for (uint32_t j=0; j<NUM_COORDS; j++){
		if(freqs.at(j) < smallFreqThresh) numSmallFreq++;
	}
	std::cout<<"Number of small eigenvalues: "<<numSmallFreq<<std::endl;
	for (uint32_t j=0; j<numSmallFreq; j++) std::cout<<" "<<to_string_prec(freqs.at(j), 1);
	std::cout<<std::endl;
	// Diagonalize the original Hessian.
	assert(diagMat(hess_in, eigenvals));
	diagInverse.fill(0.0);
	// Construct L^-1 from the eigenvalues. Ordinarily the inverse of a diagonal mtx would be just the the reciprocal
	// of its diagonal elements. Here, with the exception of the the lowest numImagKeep eigenvaues, all small
	// eigenvalues are "replaced by infinity", to yield zero when inverted. This corresponds to a physical picture of
	// having infinitely strong springs restraining any pure translation and rotation, or other near-degenerate modes
	// to zero displacement in Newton steps.
	for (uint32_t j=0; j<numImagKeep; j++)	diagInverse.at(j*NUM_COORDS + j) = 1.0/eigenvals.at(j);
	for (uint32_t j=numImagKeep; j<NUM_COORDS; j++){
		if(std::abs(freqs.at(j)) > smallFreqThresh) diagInverse.at(j*NUM_COORDS + j) = 1.0/eigenvals.at(j);
	}
	std::cout<<"Inverted diagonal:";
	for (uint32_t j=0; j<NUM_COORDS; j++) std::cout<<" "<<diagInverse.at(j*NUM_COORDS + j);
	std::cout<<std::endl;
	// Now that we have L^-1, do the matrix multiplications.
	matmul<NUM_COORDS>(diagInverse, hess_in, diagInverse, false, true); //L^-1 * Q^-1 = L^-1 * Q^T
	matmul<NUM_COORDS>(hess_in, diagInverse, diagInverse, false, false); //Q * (L^-1 * Q^-1) = H^-1
	return diagInverse;
}

template <auto NUM_COORDS>
void calcNewtonStep_inverse(const std::array<double,NUM_COORDS*NUM_COORDS>& hess, const std::array<double,NUM_COORDS>& grad,
  std::array<double,NUM_COORDS>& step, const std::array<double,NUM_COORDS/3>& masses, const uint32_t numImagKeep){
	// the pure undamped Newton step would be this vector: d = -(H^(-1) g)
	// (-1 times the inverse of the Hessian matrix times the gradient vector)
	auto invHess = invertHessian(hess, masses, numImagKeep);
	matvecmul(invHess, grad, step);
	for(auto& s : step) s = -s;
}

template <bool seekMinimum, auto NUM_COORDS>
double optimizeStepLen(const std::array<double,NUM_COORDS>& geom, const std::array<double,NUM_COORDS>& step, const double initStepFactor,
  const double rmsDisp_limit, const double maxDisp_limit, const double numderDisp){
	/// Optimize the length of an optimization step. The maximum step length is always enforced, and if a minimum is
	/// sought, the step length is crudely optimized to approximately minimize the energy after the step is taken.
	auto currStepFactor = initStepFactor;
	// Reduce step size until both RMS and max displacements are within limits.
	while ((calcRMS(arrayScale(step, currStepFactor)) > rmsDisp_limit) ||
	  (calcMaxAbs(arrayScale(step, currStepFactor)) > maxDisp_limit) ){
		currStepFactor *= 0.9;
	}
	if constexpr (seekMinimum){
		// If a smaller step yields a lower energy, keep reducing step length until it stops getting better.
		// If the step size becomes miniscule, return immediately.
		while (EvaluateEnergy_copy<false>(arrayAdd(geom, arrayScale(step, currStepFactor))) >
		EvaluateEnergy_copy<false>(arrayAdd(geom, arrayScale(step, 0.9*currStepFactor))) ){
			if (calcRMS(arrayScale(step, currStepFactor)) < 1E-12) return currStepFactor;
			currStepFactor *= 0.9;
		}
		// If a larger step yields a lower energy, keep increasing step length until it stops getting better, or
		// one of the step size limits would be exceeded by a larger step.
		while (EvaluateEnergy_copy<false>(arrayAdd(geom, arrayScale(step, currStepFactor))) >
		EvaluateEnergy_copy<false>(arrayAdd(geom, arrayScale(step, 1.1*currStepFactor))) ){
			if ( (calcRMS(arrayScale(step, 1.1*currStepFactor)) > rmsDisp_limit) ||
			(calcMaxAbs(arrayScale(step, 1.1*currStepFactor)) > maxDisp_limit) ){
				break;
			}else{
				currStepFactor *= 1.1;
			}
		}
	}else{
		// Trying to optimize Newton step length for minimum RMS grad when looking for a TS ruins convergence for 
		// some strange reason. Therefore it is commented out here.
		// while (calcRMS(calcGrad<true, false>((arrayAdd(geom, arrayScale(step, currStepFactor))),numderDisp)) > calcRMS(calcGrad<true, false>((arrayAdd(geom, arrayScale(step, 0.9*currStepFactor))),numderDisp))){
			// currStepFactor *= 0.9;
			// //std::cout<<"grad down!\n";
		// }
		// while (calcRMS(calcGrad<true, false>((arrayAdd(geom, arrayScale(step, currStepFactor))),numderDisp)) > calcRMS(calcGrad<true, false>((arrayAdd(geom, arrayScale(step, 1.1*currStepFactor))),numderDisp))){
			// if ( (calcRMS(arrayScale(step, 1.1*currStepFactor)) > rmsDisp_limit) || (calcMaxAbs(arrayScale(step, 1.1*currStepFactor)) > maxDisp_limit) ) break;
			// currStepFactor *= 1.1;
			// //std::cout<<"grad up!\n";
		// }
	}
	return currStepFactor;
}
