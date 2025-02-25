#include "emswrap_nf2ff.h"
#include <nf2ff.h>

#include <emswrap_compat.h>
#include <emswrap_csx.h>
#include <emswrap_csproperties.h>

#include <hdf5.h>
#include <highfive/highfive.hpp>

#include <stdexcept>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>

// """
// Create an nf2ff recording box. The nf2ff can either record in time-domain
// or frequency-domain. Further more certain directions and boundary condition
// mirroring can be enabled/disabled.

// :param name: str -- Name for this recording box.
// :param start/stop: (3,) array -- Box start/stop coordinates.
// :param directions: (6,) bool array -- Enable/Disables directions.
// :param mirror: (6,) int array -- 0 (Off), 1 (PEC) or 2 (PMC) boundary mirroring
// :param frequency: array like -- List of frequencies (FD-domain recording)
// """
_nf2ff::_nf2ff(ContinuousStructure* CSX, const std::string& name, const std::array<double, 3>& start, const std::array<double, 3>& stop,
               const std::vector<bool>& directions, const std::vector<int>& mirror, const std::vector<double>& frequency)
    : csx(CSX), name(name), start(start), stop(stop), directions(directions), mirror(mirror), freq(frequency)
{
    {
        int dump_type = 0;
        if (!frequency.empty()) {
            dump_type = 10;
        }
        int dump_mode = 1;

        std::string e_file = name + "_E";
        std::string h_file = name + "_H";

        // Create dumps (assuming CSX->AddDump returns a pointer to a dump object)
        e_dump = _CSX(csx).AddDump(e_file, dump_type, dump_mode, 1, frequency);
        h_dump = _CSX(csx).AddDump(h_file, dump_type + 1, dump_mode, 1, frequency);

        // Add recording boxes for each face (three dimensions)
        for (int ny = 0; ny < 3; ++ny) {
            int pos = 2 * ny;
            if (directions[pos]) {
                std::array<double, 3> l_start = start;
                std::array<double, 3> l_stop  = stop;
                l_stop[ny] = start[ny];
                _CSProperties(e_dump).AddBox(l_start, l_stop);
                _CSProperties(h_dump).AddBox(l_start, l_stop);
            }
            if (directions[pos + 1]) {
                std::array<double, 3> l_start = start;
                std::array<double, 3> l_stop  = stop;
                l_start[ny] = stop[ny];
                _CSProperties(e_dump).AddBox(l_start, l_stop);
                _CSProperties(h_dump).AddBox(l_start, l_stop);
            }
        }
    }
}

_nf2ff_results _nf2ff::CalcNF2FF(const std::string &sim_path,
                                 const std::vector<double> &freq,
                                 const std::vector<double> &theta,
                                 const std::vector<double> &phi,
                                 const std::string &outfile,
                                 bool read_cached,
                                 int verbose)
{
    return CalcNF2FF(sim_path, freq, theta, phi, 1.0, {0, 0, 0}, outfile, read_cached, verbose);
}

// """ CalcNF2FF(sim_path, freq, theta, phi, center=[0,0,0], outfile=None, read_cached=True, verbose=0):

// Calculate the far-field after the simulation is done.

// :param sim_path: str -- Simulation path
// :param freq: array like -- list of frequency for transformation
// :param theta/phi: array like -- Theta/Phi angles to calculate the far-field
// :param radius: float -- Radius to calculate the far-field (default is 1m)
// :param center: (3,) array -- phase center, must be inside the recording box
// :param outfile: str -- File to save results in. (defaults to recording name)
// :param read_cached: bool -- enable/disable read already existing results (default off)
// :param verbose: int -- set verbose level (default 0)

// :returns: nf2ff_results class instance
// """
#pragma optimize("",off)
_nf2ff_results _nf2ff::CalcNF2FF(const std::string &sim_path, const std::vector<double> &freq, const std::vector<double> &theta, const std::vector<double> &phi,
                                 double radius, const std::array<double, 3> &center, const std::string &outfile, bool read_cached, int verbose)
{
    // Determine the output file name
    std::string fn;
    if (outfile.empty())
        fn = sim_path + "/" + name + ".h5";
    else
        fn = sim_path + "/" + outfile;

    // Decide whether to perform the calculation or rely on a cached file.
    bool perform_calc = !read_cached;
    {
        std::ifstream in(fn.c_str());
        if (read_cached && !in.good())
            perform_calc = true;
    }

    // Convert theta & phi from degrees to radians
    std::vector<float> theta_rad, phi_rad;
    theta_rad.reserve(theta.size());
    phi_rad.reserve(phi.size());
    const double deg2rad = PI / 180.0;
    for (double t : theta)
        theta_rad.push_back(t * deg2rad);
    for (double p : phi)
        phi_rad.push_back(p * deg2rad);

    // Create a far-field calculator instance.
    std::vector<float> ffreq;
    for (double f : freq)
        ffreq.push_back(f);

    std::vector<float> fcenter = {static_cast<float>(center[0]), static_cast<float>(center[1]), static_cast<float>(center[2])};

    // nf2ff(vector<float> freq, vector<float> theta, vector<float> phi, vector<float> center, unsigned int numThreads=0);
    nf2ff ffcalc(ffreq, theta_rad, phi_rad, fcenter );
    ffcalc.SetVerboseLevel(verbose);

    // Set mirror conditions for each of the three dimensions.
    for (int ny = 0; ny < 3; ++ny) {
        int mirrorStart = mirror[2 * ny];
        int mirrorEnd = mirror[2 * ny + 1];

        if (mirrorStart != 0)
            ffcalc.SetMirror(mirrorStart, ny, start[ny]);

        if (mirrorEnd != 0)
            ffcalc.SetMirror(mirrorEnd, ny, stop[ny]);
    }

    // Set the radius for the far-field calculation.
    ffcalc.SetRadius(radius);

    // Process the six dump files (one for each face) if they are available.
    for (int n = 0; n < 6; ++n) {
        std::string fn_e = sim_path + "/" + name + "_E_" + std::to_string(n) + ".h5";
        std::string fn_h = sim_path + "/" + name + "_H_" + std::to_string(n) + ".h5";
        {
            std::ifstream in_e(fn_e.c_str()), in_h(fn_h.c_str());
            if (in_e.good() && in_h.good())
            {
                if (!ffcalc.AnalyseFile(fn_e, fn_h))
                {
                    throw std::runtime_error("CalcNF2FF: Unable to analyse files!");
                }
            }
        }
    }

    // Write the far-field results to an HDF5 file.
    ffcalc.Write2HDF5(fn);

    // Create the result object from the output file.
    _nf2ff_results result(fn);

    // Validate the results if data is available.
    if (!result.phi.empty()) {
        if (std::fabs((result.r - radius) / radius) >= 1e-6) {
            throw std::runtime_error("Radius does not match. Did you read an invalid cached result? Try read_cached=false");
        }
        // Assuming utilities::Check_Array_Equal takes vectors of doubles.
        // Convert result.theta and result.phi from radians to degrees for comparison.
        std::vector<double> res_theta_deg, res_phi_deg;
        res_theta_deg.reserve(result.theta.size());
        res_phi_deg.reserve(result.phi.size());
        for (double rt : result.theta)
            res_theta_deg.push_back(rt * 180.0 / PI);
        for (double rp : result.phi)
            res_phi_deg.push_back(rp * 180.0 / PI);

        if (!Check_Array_Equal(res_theta_deg, theta, 1e-4)) {
            throw std::runtime_error("Theta array does not match. Did you read an invalid cached result? Try read_cached=false");
        }
        if (!Check_Array_Equal(res_phi_deg, phi, 1e-4)) {
            throw std::runtime_error("Phi array does not match. Did you read an invalid cached result? Try read_cached=false");
        }
        if (!Check_Array_Equal(result.freq, freq, 1e-6, true)) {
            throw std::runtime_error("Frequency array does not match. Did you read an invalid cached result? Try read_cached=false");
        }
    }

    return result;

}


//////////////////////////////////////////////////////////////////////
// _nf2ff_results implementation
//////////////////////////////////////////////////////////////////////


_nf2ff_results::_nf2ff_results(const std::string& fn_) : fn(fn_) {
    load_hdf5_data();
}

void _nf2ff_results::load_hdf5_data() {
    using namespace HighFive;
    File file(fn, File::ReadOnly);

    // Load mesh data (phi, theta, r)
    {
        Group mesh = file.getGroup("Mesh");
        mesh.getDataSet("phi").read(phi);
        mesh.getDataSet("theta").read(theta);
        std::vector<double> rvec;
        mesh.getDataSet("r").read(rvec);
        if (!rvec.empty()) {
            r = rvec[0];
        }
    }

    // Load nf2ff attributes: Frequency, Dmax, Prad
    {
        Group nf2ff_group = file.getGroup("nf2ff");
        nf2ff_group.getAttribute("Frequency").read(freq);
        nf2ff_group.getAttribute("Dmax").read(Dmax);
        nf2ff_group.getAttribute("Prad").read(Prad);
    }

    // Determine dimensions from mesh data
    const size_t nTheta = theta.size();
    const size_t nPhi = phi.size();

    // Prepare storage for far-field data for each frequency point
    size_t numFreq = freq.size();
    E_theta.clear();
    E_phi.clear();
    E_norm.clear();
    E_cprh.clear();
    E_cplh.clear();
    P_rad.clear();

    // Loop through each frequency index and load datasets
    for (size_t n = 0; n < numFreq; n++) {
        // Construct dataset names (real and imaginary parts)
        std::string baseEtheta = "/nf2ff/E_theta/FD/f" + std::to_string(n);
        std::string baseEphi   = "/nf2ff/E_phi/FD/f"   + std::to_string(n);
        std::string dsEthetaRealName = baseEtheta + "_real";
        std::string dsEthetaImagName = baseEtheta + "_imag";
        std::string dsEphiRealName   = baseEphi   + "_real";
        std::string dsEphiImagName   = baseEphi   + "_imag";
        std::string dsPradName       = "/nf2ff/P_rad/FD/f" + std::to_string(n);

        // Read E_theta and E_phi real/imag parts as flat vectors.
        DataSet dsetEthetaReal = file.getDataSet(dsEthetaRealName);
        DataSet dsetEthetaImag = file.getDataSet(dsEthetaImagName);
        DataSet dsetEphiReal   = file.getDataSet(dsEphiRealName);
        DataSet dsetEphiImag   = file.getDataSet(dsEphiImagName);
        DataSet dsetPrad       = file.getDataSet(dsPradName);

        std::vector<double> EthetaRealRaw, EthetaImagRaw;
        std::vector<double> EphiRealRaw,   EphiImagRaw;
        std::vector<double> PradRaw;

        dsetEthetaReal.read(EthetaRealRaw);
        dsetEthetaImag.read(EthetaImagRaw);
        dsetEphiReal.read(EphiRealRaw);
        dsetEphiImag.read(EphiImagRaw);
        dsetPrad.read(PradRaw);

        // Expecting datasets to have shape (nPhi, nTheta).
        // We'll transpose them into (nTheta, nPhi).
        size_t totalSize = nTheta * nPhi;
        if (EthetaRealRaw.size() != totalSize || EthetaImagRaw.size() != totalSize ||
            EphiRealRaw.size()   != totalSize || EphiImagRaw.size()   != totalSize ||
            PradRaw.size()       != totalSize) {
            throw std::runtime_error("Dataset dimensions do not match expected sizes.");
        }

        // Temporary storage for transposed data
        std::vector<std::complex<double>> EthetaTrans(totalSize);
        std::vector<std::complex<double>> EphiTrans(totalSize);
        std::vector<double> PradTrans(totalSize);

        for (size_t i = 0; i < nPhi; i++) {
            for (size_t j = 0; j < nTheta; j++) {
                // Original index: (i, j) -> i * nTheta + j
                // Transposed index: (j, i) -> j * nPhi + i
                size_t origIdx = i * nTheta + j;
                size_t transIdx = j * nPhi + i;
                EthetaTrans[transIdx] = std::complex<double>(EthetaRealRaw[origIdx], EthetaImagRaw[origIdx]);
                EphiTrans[transIdx]   = std::complex<double>(EphiRealRaw[origIdx],   EphiImagRaw[origIdx]);
                PradTrans[transIdx]   = PradRaw[origIdx];
            }
        }

        // Reshape flat vectors into 2D (vector of rows, each row has nPhi elements)
        std::vector<std::complex<double>> EthetaRow;
        std::vector<std::complex<double>> EphiRow;
        std::vector<double> PradRow;
        std::vector<std::vector<std::complex<double>>> EthetaMatrix;
        std::vector<std::vector<std::complex<double>>> EphiMatrix;
        std::vector<std::vector<double>> PradMatrix;
        for (size_t j = 0; j < nTheta; j++) {
            size_t offset = j * nPhi;
            EthetaRow.assign(EthetaTrans.begin() + offset, EthetaTrans.begin() + offset + nPhi);
            EphiRow.assign(EphiTrans.begin() + offset, EphiTrans.begin() + offset + nPhi);
            PradRow.assign(PradTrans.begin() + offset, PradTrans.begin() + offset + nPhi);
            EthetaMatrix.push_back(EthetaRow);
            EphiMatrix.push_back(EphiRow);
            PradMatrix.push_back(PradRow);
        }

        // Compute E_norm, E_cprh, E_cplh for each cell.
        std::vector<std::vector<double>> E_normMatrix;
        std::vector<std::vector<std::complex<double>>> E_cprhMatrix;
        std::vector<std::vector<std::complex<double>>> E_cplhMatrix;
        const double normFactor = 1.0 / std::sqrt(2.0);
        for (size_t j = 0; j < nTheta; j++) {
            std::vector<double> normRow;
            std::vector<std::complex<double>> cprhRow, cplhRow;
            for (size_t i = 0; i < nPhi; i++) {
                // Calculate norm = sqrt(|E_theta|^2 + |E_phi|^2)
                double norm_val = std::sqrt(
                    std::norm(EthetaMatrix[j][i]) + std::norm(EphiMatrix[j][i])
                );
                normRow.push_back(norm_val);

                // Use phi value corresponding to current column
                double phi_val = phi[i];
                std::complex<double> phaseFactor_pos = std::complex<double>(std::cos(phi_val), std::sin(phi_val));
                std::complex<double> phaseFactor_neg = std::complex<double>(std::cos(phi_val), -std::sin(phi_val));
                std::complex<double> sumE = EthetaMatrix[j][i] + std::complex<double>(0,1)*EphiMatrix[j][i];
                std::complex<double> diffE = EthetaMatrix[j][i] - std::complex<double>(0,1)*EphiMatrix[j][i];
                cprhRow.push_back(phaseFactor_pos * sumE * normFactor);
                cplhRow.push_back(phaseFactor_neg * diffE * normFactor);
            }
            E_normMatrix.push_back(normRow);
            E_cprhMatrix.push_back(cprhRow);
            E_cplhMatrix.push_back(cplhRow);
        }

        // Append data for this frequency index.
        E_theta.push_back(EthetaMatrix);
        E_phi.push_back(EphiMatrix);
        E_norm.push_back(E_normMatrix);
        E_cprh.push_back(E_cprhMatrix);
        E_cplh.push_back(E_cplhMatrix);
        P_rad.push_back(PradMatrix);
    }
}
#pragma optimize("",on)