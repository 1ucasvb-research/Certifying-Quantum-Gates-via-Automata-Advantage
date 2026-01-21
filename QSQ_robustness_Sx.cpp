// Code used in "Certifying Quantum Gates via Automata Advantage"
// by Anna Schroeder, Lucas B. Vieira, Jan Nöller, Nikolai Miklin, and Mariami Gachechiladze
// https://arxiv.org/abs/2510.09575
//
// Written by Lucas Vieira
//
// Generates a density plot for the QSQ protocol on the S = sqrt(Z) gate by sampling over arbitrary qubit CPTP maps.
// The code produces a simple .bin file consisting of a raw array of floats of size (DENSITY_PLOT_W, DENSITY_PLOT_H),
// where each entry corresponds to a (pfail, inFid) coordinate with the accumulated samples.
// For simplicity, we rotate to the Sx = sqrt(X) basis.
// Only requires Armadillo and tinyformat (included). Can be compiled directly with g++ using -std=c++17 -larmadillo flags.
//
// License: LGPL v2.1
#include <iostream>
#include <armadillo>
#include <cmath>
#include <complex>
#include <string>
#include <fstream>
#include <vector>
#include <filesystem>
#include "tinyformat.h"

typedef arma::mat    RMatrix;
typedef arma::cx_mat CMatrix;
typedef arma::vec    RVector;
typedef arma::cx_vec CVector;
typedef arma::uvec   UVector;

// tau = 2*pi
constexpr double TAU = 6.283185307179586;

// Complex constants for simplicity
constexpr std::complex C0  = std::complex(0.0, 0.0);
constexpr std::complex C1  = std::complex(1.0, 0.0);
constexpr std::complex CI  = std::complex(0.0, 1.0);
constexpr std::complex C11p = std::complex(1.0, 1.0);
constexpr std::complex C11m = std::complex(1.0, -1.0);

// Standard matrices and vectors
// Identity on a qubit
CMatrix I = arma::eye<CMatrix>( 2, 2 );
// Maximally entangled state for Choi matrix
CMatrix PhiPlus = arma::vectorise(arma::eye<CMatrix>( 2, 2 ));

// Complex magnitude bias in the generation of Kraus operators
// See random_complex_matrix() below for details
constexpr double KRAUS_BIAS = 2.0; // 2.0 seems good, around 4.0 has artifacts

// Resolution of the final density plot
constexpr int DENSITY_PLOT_W = 4000;
constexpr int DENSITY_PLOT_H = 3000;

// Range of pfail and inFid to consider within the plot region
// Essentially, the rectangle width and height, focused at the origin
// Note that we cannot, a priori, sample channels within this rectangle, so smaller rectangles
// lead to many rejected samples. These were the values used in the paper.
constexpr float MAX_PFAIL = 0.35f;
constexpr float MAX_INFID = 0.15f;

// Number of trials for infidelity optimization
constexpr int INFIDELITY_OPT_ATTEMPTS = 5;
// Number of iterations in the infidelity optimization
constexpr int INFIDELITY_OPT_ITERS = 400;
// Power k for the decay curve (1 - t)^k, where t is the iteration fraction from [0,1]
constexpr double INFIDELITY_OPT_POWER = 4.0;

// How many samples to print/save current progress
constexpr int PROGRESS_PRINT = 1'000;
constexpr int PROGRESS_SAVE = 10'000;

// Generates a Haar-random unitary matrix from U(n)
// Reference:
// "How to generate random matrices from the classical compact groups"
// Francesco Mezzadri
// https://arxiv.org/abs/math-ph/0609050
inline CMatrix random_unitary(int n) {
	CMatrix Z = arma::randn<CMatrix>(n, n) / sqrt(2.0);
	CMatrix Q, R;
	arma::qr_econ(Q, R, Z);
	CVector r_diag = R.diag();
	CVector phases = arma::conv_to<CVector>::from(r_diag / arma::abs(r_diag));
	Q.each_row() %= phases.st();
	return Q;
}

// Generates a k-biased random complex matrix of size n x m
// For each entry, generate a uniformly random magnitude r = [0,1) and phase theta = [0,2pi)
// The complex entry is z = r^k exp(i theta)
// We opted for this method as the usual Ginibre ensemble sampling does not sample well the extremal
// behaviors in the QSQ protocol we considered
inline CMatrix random_complex_matrix(int n, int m, double k) {
	return arma::pow(arma::randu<RMatrix>(n, m), k) % arma::exp(CI * TAU * arma::randu<RMatrix>(n, m));
}

// Quick and dirty partial trace of for the action of the Choi state
// Not general, but faster than a loop
// Note: We consider the Hilbert space convention of H_Choi = H_in (x) H_out
inline CMatrix partial_trace_Choi_input(CMatrix M) {
	CMatrix res(2, 2);
	res(0, 0) = M(0,0) + M(2,2);
	res(0, 1) = M(0,1) + M(2,3);
	res(1, 0) = M(1,0) + M(3,2);
	res(1, 1) = M(1,1) + M(3,3);
	return res;
}

// Same for link product
// Note: In principle, we could speed up the code further by using the superoperator representation,
// avoiding the link product altogether.
CMatrix link_product(CMatrix M1, CMatrix M2) {
	CMatrix res(4, 4);
	res(0,0) = M1(0,0) * M2(0,0) + M1(0,1) * M2(0,2) + M1(1,0) * M2(2,0) + M1(1,1) * M2(2,2);
	res(0,1) = M1(0,0) * M2(0,1) + M1(0,1) * M2(0,3) + M1(1,0) * M2(2,1) + M1(1,1) * M2(2,3);
	res(0,2) = M1(0,2) * M2(0,0) + M1(0,3) * M2(0,2) + M1(1,2) * M2(2,0) + M1(1,3) * M2(2,2);
	res(0,3) = M1(0,2) * M2(0,1) + M1(0,3) * M2(0,3) + M1(1,2) * M2(2,1) + M1(1,3) * M2(2,3);
	res(1,0) = M1(0,0) * M2(1,0) + M1(0,1) * M2(1,2) + M1(1,0) * M2(3,0) + M1(1,1) * M2(3,2);
	res(1,1) = M1(0,0) * M2(1,1) + M1(0,1) * M2(1,3) + M1(1,0) * M2(3,1) + M1(1,1) * M2(3,3);
	res(1,2) = M1(0,2) * M2(1,0) + M1(0,3) * M2(1,2) + M1(1,2) * M2(3,0) + M1(1,3) * M2(3,2);
	res(1,3) = M1(0,2) * M2(1,1) + M1(0,3) * M2(1,3) + M1(1,2) * M2(3,1) + M1(1,3) * M2(3,3);
	res(2,0) = M1(2,0) * M2(0,0) + M1(2,1) * M2(0,2) + M1(3,0) * M2(2,0) + M1(3,1) * M2(2,2);
	res(2,1) = M1(2,0) * M2(0,1) + M1(2,1) * M2(0,3) + M1(3,0) * M2(2,1) + M1(3,1) * M2(2,3);
	res(2,2) = M1(2,2) * M2(0,0) + M1(2,3) * M2(0,2) + M1(3,2) * M2(2,0) + M1(3,3) * M2(2,2);
	res(2,3) = M1(2,2) * M2(0,1) + M1(2,3) * M2(0,3) + M1(3,2) * M2(2,1) + M1(3,3) * M2(2,3);
	res(3,0) = M1(2,0) * M2(1,0) + M1(2,1) * M2(1,2) + M1(3,0) * M2(3,0) + M1(3,1) * M2(3,2);
	res(3,1) = M1(2,0) * M2(1,1) + M1(2,1) * M2(1,3) + M1(3,0) * M2(3,1) + M1(3,1) * M2(3,3);
	res(3,2) = M1(2,2) * M2(1,0) + M1(2,3) * M2(1,2) + M1(3,2) * M2(3,0) + M1(3,3) * M2(3,2);
	res(3,3) = M1(2,2) * M2(1,1) + M1(2,3) * M2(1,3) + M1(3,2) * M2(3,1) + M1(3,3) * M2(3,3);
	return res;
}

// Computes the probabilities for Lambda^2k, k = 0,...,L in terms of the Choi of Lambda
// rho0 is the initial state, Ms is the POVM
// Note that we assume the correct output alternates between M0 and M1, as expected for Sx
inline void get_probabilities(int L, double* probs, CMatrix rho0, CMatrix* Ms, CMatrix Choi) {
	CMatrix rho = rho0;
	CMatrix Choi2 = link_product(Choi, Choi);
	for (int i = 0; i <= L; ++i) {
		probs[i] = std::real(arma::trace(rho * Ms[i % 2]));
		rho = partial_trace_Choi_input(Choi2 * arma::kron(rho.st(), I));
	}
}

// Generates a random unitary and return its Choi state
inline CMatrix get_general_unitary_channel() {
	CMatrix U = random_unitary(2);
	CMatrix v = arma::kron(I, U) * PhiPlus;
	return v * v.t();
}

// Generates a random quantum channel with Kraus rank between min_rank and max_rank (inclusive)
// Reference:
// "Generating random quantum channels"
// Ryszard Kukulski, Ion Nechita, Łukasz Pawela, Zbigniew Puchała, Karol Życzkowski
// https://arxiv.org/abs/2011.02994
//
// Note, however, that we replace the original Ginibre ensemble with our biased sampler
// See random_complex_matrix() for more details
inline CMatrix get_general_noisy_channel(int min_rank, int max_rank) {
	int nk = min_rank + int(std::floor(arma::randu() * (double)(max_rank - min_rank + 1)));
	std::vector<CMatrix> Ks;
	
	// For nk = 1, we just have a unitary so generate that directly instead
	if (nk == 1) return get_general_unitary_channel();
	
	// If nk >= 1, we have CPTP maps, so we need to put in some effort
	for (int i = 0; i < nk; ++i) {
		Ks.push_back(random_complex_matrix(2, 2, KRAUS_BIAS));
	}
	CMatrix E = CMatrix(2, 2, arma::fill::zeros);
	for (int i = 0; i < nk; ++i) {
		E += Ks[i].t() * Ks[i];
	}
	E = arma::inv(arma::sqrtmat(E));
	CMatrix Choi = CMatrix(4, 4, arma::fill::zeros);
	for (int i = 0; i < nk; ++i) {
		Ks[i] = Ks[i] * E;
		CMatrix v = arma::kron(I, Ks[i]) * PhiPlus;
		Choi += v * v.t();
	}
	return Choi;
}

// Computes the average channel infidelity between target gate G with unitary gauge U and noisy channel
//   inFid = (d/d+1)(1 - Tr[Choi[U.G.U'].Choi[Lambda]]/d^2)
// Reference:
// "Theory of Quantum System Certification"
// Martin Kliesch, and Ingo Roth
// https://journals.aps.org/prxquantum/pdf/10.1103/PRXQuantum.2.010201
float infidelity_with_gauge_unitary(CMatrix &G, CMatrix &ChoiNoisy, CMatrix U) {
	CMatrix ChoiIdeal_gauged = arma::kron(I, U * G * U.t()) * PhiPlus;
	ChoiIdeal_gauged = ChoiIdeal_gauged * ChoiIdeal_gauged.t();
	return (2.0f/3.0f)*(1.0f - std::real(arma::trace(ChoiNoisy * ChoiIdeal_gauged)) / 4.0f);
};

// Minimizes worst-case model infidelity based on states, measurement and gate under unitary gauge.
// Optimization over unitaries is done via a modified version of the Luus-Jaakola method,
// where we use a polynomial decay (as opposed to exponential) for better control.
//
// Reference:
// "Optimization by direct search and systematic reduction of the size of search region"
// Rein Luus and T. H. I. Jaakola
// https://aiche.onlinelibrary.wiley.com/doi/10.1002/aic.690190413
//
// Note: since here we consider pure states and measurements diagonal in the computational basis
// to analyze the Sx = sqrt(X) gate, it is easy to show that the state and measurement contribution
// to the model infidelity reduce to (1 - |W_00|^2) and sqrt(1 - |W_00|^2), for W the gauge unitary
float get_infidelity_full(CMatrix &G, CMatrix &ChoiNoisy, int iterations) {
	using namespace std;
	float best = 1.0f;
	float v;
	CMatrix W; // the gauge unitary
	CMatrix bestW = random_unitary(2);
	for (int i = 0; i < iterations; ++i) {
		// Polynomial decay (1 - t)^4 leads to a more controlled convergence rate
		W = bestW * arma::powmat(random_unitary(2), pow(1.0 - (double)i/(double)iterations, INFIDELITY_OPT_POWER));
		v = max(
			infidelity_with_gauge_unitary(G, ChoiNoisy, W),
			(float)sqrt(1.0 - pow(abs(W(0,0)), 2.0))
		);
		if (v < best) {
			best = v;
			bestW = W;
		}
	}
	return best;
}

// Saves the accumulated plot data into the file path given by fn
// The output is raw binary file of W*H floats
void save_plot(const char* fn, float* plot) {
	std::ofstream file(fn, std::ios::binary);
	std::cout << tfm::format("\nSaving to '%s'... ", fn);
	file.write(reinterpret_cast<const char*>(plot), DENSITY_PLOT_W*DENSITY_PLOT_H*sizeof(float));
	file.close();
	std::cout << "Done\n";
}

// Loads the plot data from path given by fn, allowing us to resume computation
bool load_plot(const char* fn, float* plot) {
	std::ifstream file(fn, std::ios::binary);
	if (!file) {
		return false;
	}
	std::cout << tfm::format("\nLoading from '%s'... ", fn);
	file.read(reinterpret_cast<char*>(plot), DENSITY_PLOT_W*DENSITY_PLOT_H*sizeof(float));
	file.close();
	std::cout << "Done\n";
	return true;
}

// Computes the total value within the plot, essentially the number of samples collected so far.
// Since we perform sub-pixel interpolation, this sum is not exact, hence the rounding.
uint64_t sum_plot(float* plot) {
	double sum = 0.0;
	for (int i = 0; i < DENSITY_PLOT_W*DENSITY_PLOT_H; ++i) {
		sum += plot[i];
	}
	return (uint64_t)std::round(sum);
}

int main(int argc, char* argv[]) {
	if (argc <= 3) {
		std::string prog_name = std::filesystem::path(argv[0]).filename().string();
		std::cout << "Usage: " << prog_name << " <L> <num_samples> <filename>" << std::endl;
		std::cout << std::flush;
		exit(0);
	}
	
	int L = (int)atoi(argv[1]);
	uint64_t NUM_SAMPLES = (uint64_t)atoll(argv[2]);
	char* name = argv[3];
	
	// Initialize Sx matrix as our target gate G
	CMatrix G = {{C11p, C11m}, {C11m, C11p}};
	G /= 2.0;
	
	// Initialize POVM {|0><0|, |1><1|} and initial state |0><0|
	CMatrix Ms[2] = {
		{{C1,C0},{C0,C0}},
		{{C0,C0},{C0,C1}}
	};
	CMatrix rho = {{C1,C0},{C0,C0}};
	
	// Vector of probabilities for Lambda^2k, k=0,...,L
	double* probs = (double*)calloc(L+1, sizeof(double));
	
	// Initialize plot data array
	std::string plotfn = tfm::format("data_%dx%d_L%d_%s.bin", DENSITY_PLOT_W, DENSITY_PLOT_H, L, name);
	float* plot = (float*)calloc(DENSITY_PLOT_W*DENSITY_PLOT_H, sizeof(float));
	// Load existing data if filename already exists, and resume
	if (load_plot(plotfn.c_str(), plot)) {
		std::cout << tfm::format("Successfully loaded '%s'", plotfn) << std::endl;
		std::cout << std::flush;
	}
	// Current number of samples
	uint64_t s = sum_plot(plot);
	
	// Some utility counters
	uint64_t k = 0;
	uint64_t ksave = 0;
	
	float inFid, pfail;
	std::cout << "Started..." << std::endl;
	while (s < NUM_SAMPLES) {
		// Sample a random noisy channel
		// Kraus rank from 1 (unitary) to 4 cover all channels
		CMatrix ChoiNoisy = get_general_noisy_channel(1, 4);
		
		// Compute probabilities and final pfail
		get_probabilities(L, probs, rho, Ms, ChoiNoisy);
		pfail = 0.0f;
		for (int i = 0; i <= L; ++i) {
			pfail += probs[i];
		}
		pfail /= (float)(L+1);
		pfail = 1.0f - pfail;
		
		// If outside the target rectangle, reject and sample again
		if (pfail > MAX_PFAIL) continue;
		
		{
			// x coordinate within plot
			int x = (int)std::floor((float)(DENSITY_PLOT_W-1) * (pfail / MAX_PFAIL));
			if (x < 0 || x >= DENSITY_PLOT_W) continue; // safety check
			
			// Optimize infidelity multiple times, keep track of best value
			float inFid = 1.0;
			for (int i = 0; i < INFIDELITY_OPT_ATTEMPTS; ++i) {
				inFid = std::min(inFid, get_infidelity_full(
					G, ChoiNoisy,
					INFIDELITY_OPT_ITERS
				));
			}
			
			// Skip this sample if outside target rectangle
			if (inFid > MAX_INFID) continue;
			
			// y coordinate within the plot
			int y = (int)std::floor((float)(DENSITY_PLOT_H-1) * (1.0 - inFid / MAX_INFID));
			if (y < 0 || y >= DENSITY_PLOT_H) continue; // safety check
			
			// To improve resolution, we apply a simple sub-pixel interpolation of the data point
			float fx = (float)(DENSITY_PLOT_W-1) * (pfail / MAX_PFAIL) - (float)x;
			float fy = (float)(DENSITY_PLOT_H-1) * (1.0 - inFid / MAX_INFID) - (float)y;
			// We accumulate the value accordingly
			plot[y*DENSITY_PLOT_W + x] += (1.0f - fx)*(1.0f - fy);
			plot[y*DENSITY_PLOT_W + (x+1)] += fx*(1.0f - fy);
			plot[(y+1)*DENSITY_PLOT_W + x] += (1.0f - fx)*fy;
			plot[(y+1)*DENSITY_PLOT_W + (x+1)] += fx*fy;
		}
		
		// Print progress every PROGRESS_PRINT samples and save every PROGRESS_SAVE samples
		++k;
		if (k == PROGRESS_PRINT) {
			std::cout << tfm::format("\r%0.06f", (double)(s+1) / (double)NUM_SAMPLES);
			std::cout << std::flush;
			k = 0;
		}
		++ksave;
		if (ksave == PROGRESS_SAVE) {
			save_plot(plotfn.c_str(), plot);
			ksave = 0;
		}
		
		s++;
	}
	
	// Final save and complete
	save_plot(plotfn.c_str(), plot);
	std::cout << "Completed " << tfm::format("'%s'", plotfn) << std::endl;
	return 0;
}