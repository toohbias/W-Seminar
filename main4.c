#include <math.h>
#include "plot.h"

#define MEGA 1E6
#define KILO 1E3
#define MILLI 1E-3
#define MIKRO 1E-6
#define NANO 1E-9
#define PIKO 1E-12

#define DT (1 * NANO)
// #define TF (1 * MIKRO)
#define TF (10 * MILLI)

// Randwerte
const double R1 = 5 * KILO;
const double R2 = 100 * KILO;
const double R3 = 1 * MEGA;

const double C1 = 10 * NANO;
const double C2 = 1 * NANO;
const double C3 = 1 * NANO;

const double VD = 0.7;

const double Vnf_0 = 3.3;
const double Fnf = 1 * KILO;
double Vnf(double t) { return sin(2 * M_PI * Fnf * t); }
const double Vhf_0 = 5.0;
const double Fhf = 1.2 * MEGA;
double Vhf(double t) { return sin(2 * M_PI * Fhf * t); }
double Vam(double t) { return (Vhf_0 + Vnf_0 * Vnf(t)) * Vhf(t); }

double V1_n, V1_n1 = 0;
double V2_n, V2_n1 = 0;
double V3_n, V3_n1 = 0;

double Ir1, Ic1 = 0;

int main(void) {

    char x_desc[] = "Zeit t (ms)";

    Plot plotV0 = initPlot(0.0, TF, -10, 10, x_desc, NULL, NULL, "v0.pbm");
    Plot plotV1 = initPlot(0.0, TF, -10, 10, x_desc, NULL, NULL, "v1.pbm");
    Plot plotV2 = initPlot(0.0, TF, -10, 10, x_desc, NULL, NULL, "v2.pbm");
    Plot plotV3 = initPlot(0.0, TF, -10, 10, x_desc, NULL, NULL, "v3.pbm");

    for(double t = 0; t < TF; t += DT) {
        setPoint(t, Vam(t), &plotV0);

        if(Vam(t) - V1_n > VD) {
            
            // ON MODE

            V1_n1 = Vam(t) - VD;
            V2_n1 = V2_n + (((V1_n1 - V2_n) / R2) - (V3_n / R3)) / C2 * DT;
            V3_n1 = V3_n + (((C2 + C3) * (-V3_n / R3) + C3 * ((V1_n1 - V2_n1) / R2)) / (C2 * C3)) * DT;

        } else {

            // OFF MODE

            V1_n1 = V1_n + ((-V1_n / R1 - (V1_n - V2_n) / R2) / C1) * DT;
            V2_n1 = V2_n + ((((V1_n1 - V2_n) / R2) - (V3_n / R3)) / C2) * DT;
            V3_n1 = V3_n + (((C2 + C3) * (-V3_n / R3) + C3 * ((V1_n1 - V2_n1) / R2)) / (C2 * C3)) * DT;

        }

        V1_n = V1_n1;
        V2_n = V2_n1;
        V3_n = V3_n1;

        setPoint(t, V1_n, &plotV1);
        setPoint(t, V2_n, &plotV2);
        setPoint(t, V3_n, &plotV3);

    }

    savePlot(&plotV0);
    savePlot(&plotV1);
    savePlot(&plotV2);
    savePlot(&plotV3);
}
