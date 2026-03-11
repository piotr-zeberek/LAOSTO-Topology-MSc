#include <iostream>
#include <iomanip>
#include <fstream>

#include "ToyModel.h"
#include "LAOSTO.h"

#include "System2DCalculations.h"
#include "System2DCalculationsPrinter.h"

int main(int argc, char *argv[])
{
    // ToyModel sys;
    LAOSTO sys;

    System2DCalculations calc(sys);
    System2DCalculationsPrinter printer(calc);

    double LAOSTO_low_band = -47.5031092363;    // xy, Bzc = 0.24 T
    double LAOSTO_bottom_band = -2.83022409703; // yz, Bzc = 0.73 T
    double LAOSTO_top_band = 3.33333333333;     // xz, Bzc = 0.14 T

    std::vector<double> LAOSTO_mus = {
        LAOSTO_low_band,
        LAOSTO_bottom_band,
        LAOSTO_top_band};

    std::vector<std::string> LAOSTO_mus_labels = {
        "xy",
        "yz",
        "xz"};

    sys.mu = meV2au(0.0);
    sys.Bx = T2au(0.0);
    sys.By = T2au(0.0);
    sys.Bz = T2au(0.0);

    // example
    printer.printBandStructureSlice_normal_orbital_type("data/dispersion_normal.dat",
                                                        Eigen::VectorXd::LinSpaced(2001, -0.6, 0.6), 0);
    return 0;
}