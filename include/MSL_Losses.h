#pragma once

#include <emswrap.h>
#include <thread>

// Ported from MSL_Losses.m
void MSL_Losses() {
    std::string Sim_Path = "MSL_Losses";

    // physical constants and simulation parameters
    const double unit = 1e-6;          // all lengths in [um]
    const double MSL_length = 10000;
    const double MSL_port_dist = 5000;
    const double MSL_width = 225;
    const double MSL_conductivity = 41e6;
    const double MSL_thickness = 35e-6;
    const double substrate_thickness = 250;
    const double substrate_epr = 9.8;

    const double f_start = 0e9;
    const double f_stop  = 25e9;
    const double f0 = 0.5*(f_start+f_stop);
    const double fc = 0.5*(f_stop-f_start);
    const double lambda = C0/f_stop;

    // initialize FDTD simulation
    _openEMS FDTD;
    FDTD.SetEndCriteria(1e-4);
    FDTD.SetGaussExcite(f0, fc);
    std::vector<std::string> BC = { "PML_8", "PML_8", "PML_8", "PML_8", "PEC", "PML_8" };
    FDTD.SetBoundaryCond(BC);

    // setup CSX structure and mesh grid
    _CSX CSX(new ContinuousStructure());
    FDTD.SetCSX(CSX);
    _CSRectGrid mesh(CSX.GetGrid());
    mesh.SetDeltaUnit(unit);
    double resolution = C0 / (f_stop*sqrt(substrate_epr)) / unit / 20;

    // mesh in x-direction
    mesh.SmoothMeshLines('x', { -MSL_length * 0.5 - MSL_port_dist, 0, MSL_length * 0.5 + MSL_port_dist }, resolution, 1.3);

    // mesh in y-direction
    num_vector half_y_lines = SmoothMesh::SmoothMeshLines2({0, MSL_width / 2}, resolution / 6, 1.3);
    num_vector y_lines(-0.5 * lambda / unit, -half_y_lines, half_y_lines, 0.5 * lambda / unit);
    mesh.SmoothMeshLines('y', y_lines, resolution, 1.4);

    // mesh in z-direction: substrate divided into 10 parts, plus an extra point at 0.5*lambda/unit
    num_vector z_lines(linspace(0.0, substrate_thickness, 10), 0.5 * lambda / unit);
    mesh.SmoothMeshLines('z', z_lines, resolution);


    // add substrate material and geometry
    auto substrate = CSX.AddMaterial("RO4350B");
    substrate.SetMaterialProperty("epsilon", substrate_epr);
    std::array<double, 3> sub_start = { mesh.GetLines('x').front(), mesh.GetLines('y').front(), 0 };
    std::array<double, 3> sub_stop  = { mesh.GetLines('x').back(),  mesh.GetLines('y').back(),  substrate_thickness };
    substrate.AddBox(0, sub_start, sub_stop);

    // add conducting sheet for microstrip line
    auto goldSheet = CSX.AddConductingSheet("gold", MSL_conductivity, MSL_thickness);

    // define MSL port 1 (excited port)
    std::array<double, 3> port1_start = { mesh.GetLines('x').front(), -MSL_width/2, substrate_thickness };
    std::array<double, 3> port1_stop  = { mesh.GetLines('x').front() + MSL_port_dist, MSL_width/2, 0 };

    auto port1 = CSX.AddMSLPort(999, 1, goldSheet, port1_start, port1_stop, 0,
                                {0, 0, -1},
                                true, 10 * resolution, MSL_port_dist);

    // define MSL port 2 (measurement port)
    std::array<double, 3> port2_start = { mesh.GetLines('x').back(), -MSL_width/2, substrate_thickness };
    std::array<double, 3> port2_stop  = { mesh.GetLines('x').back() - MSL_port_dist, MSL_width/2, 0 };
    auto port2 = CSX.AddMSLPort(999, 2, goldSheet, port2_start, port2_stop, 0,
                                {0, 0, -1},
                                false, 0.0, MSL_port_dist);

    // add the microstrip line (between the two ports)
    std::array<double, 3> line_start = { mesh.GetLines('x').front() + MSL_port_dist, -MSL_width/2, substrate_thickness };
    std::array<double, 3> line_stop  = { mesh.GetLines('x').back() - MSL_port_dist, MSL_width/2, substrate_thickness };
    _CSProperties(goldSheet).AddBox(500, line_start, line_stop);

    auto myd1 = CSX.AddDump("EDump", 0, 0, 0);
    _CSProperties(myd1).AddBox(mesh.GetSimArea()[0], mesh.GetSimArea()[1]);

    auto myd2 = CSX.AddDump("HDump", 1, 0, 0);
    _CSProperties(myd2).AddBox(mesh.GetSimArea()[0], mesh.GetSimArea()[1]);

    auto myd3 = CSX.AddDump("CurrentDump", 2, 0, 0);
    _CSProperties(myd3).AddBox(mesh.GetSimArea()[0], mesh.GetSimArea()[1]);

    // write simulation files and run openEMS
    FDTD.DebugCSX();
    FDTD.DebugBox();
    FDTD.DebugMaterial();
    FDTD.DebugOperator();
    FDTD.DebugPEC();
    FDTD.SetVerboseLevel(3);

    auto openemsArgs = {std::string("numThreads=") + std::to_string(std::thread::hardware_concurrency())};
    FDTD.Run(Sim_Path, false, false, openemsArgs);

    // post-processing: calculate S-parameters
    std::vector<double> f;
    int numF = 1601;
    for (int i = 0; i < numF; i++)
        f.push_back(f_start + (f_stop - f_start) * i / double(numF-1));

    std::vector<Port*> ports = { &port1, &port2 };
    calcPort(ports, Sim_Path, {f}, 50);

    // print S21 attenuation in dB versus frequency
    for (size_t i = 0; i < f.size(); i++) {
        // compute S21 = port2 uf.ref / port1 uf.inc (assuming uf holds complex vector data)
        std::complex<double> s21 = ports[1]->uf_ref[i] / ports[0]->uf_inc[i];
        double s21_dB = -20 * log10(std::abs(s21));
        std::cout << f[i]/1e9 << " GHz: " << s21_dB << " dB" << std::endl;
    }

    // write simulation data to CSV file for plotting
    {
        // CSV file for FDTD simulated S21 attenuation data
        std::string csvFilename = std::format("{0}/msl_data.csv", Sim_Path);
        std::ofstream csvFile(csvFilename);
        if (!csvFile)
            throw std::runtime_error("Failed to open file: " + csvFilename);
        csvFile << "frequency_GHz,s21_dB\n";
        for (size_t i = 0; i < f.size(); i++) {
            std::complex<double> s21 = ports[1]->uf_ref[i] / ports[0]->uf_inc[i];
            double s21_dB = -20 * log10(std::abs(s21));
            csvFile << (f[i] / 1e9) << "," << s21_dB << "\n";
        }
        csvFile.close();

        // CSV file for loss model data
        std::string csvLossFilename = std::format("{0}/loss_model.csv", Sim_Path);
        std::ofstream csvLoss(csvLossFilename);
        if (!csvLoss)
            throw std::runtime_error("Failed to open file: " + csvLossFilename);
        csvLoss << "frequency_GHz,loss_dB\n";
        std::vector<double> model_f = { 1, 2, 2.5, 3, 4, 5, 7.5, 10, 12.5, 15, 17.5, 20, 25 };
        std::vector<double> model_loss = { 3.0, 4.2, 4.7, 5.2, 5.9, 6.6, 8.1, 9.38, 10.5, 11.5, 12.4, 13.2, 14.65 };
        for (size_t i = 0; i < model_f.size(); i++) {
            double lossVal = model_loss[i] * MSL_length * unit;
            csvLoss << model_f[i] << "," << lossVal << "\n";
        }
        csvLoss.close();

        PyPlot("msl_plot", std::format(R"(
# load simulation and loss model data
data = pd.read_csv('{0}/msl_data.csv')
loss = pd.read_csv('{0}/loss_model.csv')

plt.figure()
plt.plot(data['frequency_GHz'], data['s21_dB'], 'r--', linewidth=2, label='FDTD simulated attenuation')
plt.plot(loss['frequency_GHz'], loss['loss_dB'], 'k-', linewidth=1, label='t=35um, loss model by E. Hammerstad & F. Bekkadal')
plt.grid(True)
plt.xlabel('frequency (GHz)')
plt.ylabel('-|S_{{21}}| (dB)')
plt.legend(loc='upper left')
plt.savefig('{0}/plot.png')
    )", Sim_Path));

    }
}