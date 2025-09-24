#include <math.h>
#include "plot.h"

#define MEGA 1E6
#define KILO 1E3
#define MILLI 1E-3
#define MIKRO 1E-6
#define NANO 1E-9
#define PIKO 1E-12

#define DT (1 * NANO)
// #define TF (100 * MIKRO)
#define TF (2 * MILLI)

// Randwerte
const double R1 = 10 * KILO;
const double R2 = 1 * KILO;
const double R3 = 10 * KILO;
const double R4 = 1 * KILO;
const double R5 = 10 * KILO;

const double Rth = 1 / (1/R1 + 1/R2 + 1/R3);

const double C1 = 100 * NANO;
const double C2 = 1.8 * NANO;

const double L1 = 10 * MIKRO;

const double beta = 200;
const double VBE = 0.7;
const double VCE = 0.2;
const double VCB = -0.5;

const double Vcc = 5.0;
const double Vnf_0 = 3.0;
const double Fnf = 1 * KILO;
double Vnf(double t) { return Vnf_0 * sin(2.0 * M_PI * Fnf * t); }
const double Vhf_0 = 1.5;
const double Fhf = 1.2 * MEGA;
double Vhf(double t) { return Vhf_0 * sin(2.0 * M_PI * Fhf * t); }

double Vth(double t) { return (Vnf(t)/R1 + Vhf(t)/R2) * Rth; }

double V1_n, V1_n1 = 0;
double V2_n, V2_n1 = 0;
double V3_n, V3_n1 = Vcc;
double V4_n, V4_n1 = 0;
double dV4 = 0;

double I1, I2, I3u, I4l = 0;

int main(void) {

    char x_desc[] = "Zeit t (ms)";

    Plot plotV0 = initPlot(0.0, TF, -3, 3, x_desc, "Spannung am Punkt 1 (V)", NULL, "v0.pbm");
    Plot plotV1 = initPlot(0.0, TF, -3, 3, x_desc, NULL, NULL, "v1.pbm");
    Plot plotV2 = initPlot(0.0, TF, -2.0, 2.0, x_desc, NULL, NULL, "v2.pbm");
    Plot plotV3 = initPlot(0.0, TF, -1.0, 5.0, x_desc, NULL, NULL, "v3.pbm");
    Plot plotV3b = initPlot(0.0, TF, -5, 5, x_desc, NULL, NULL, "v3b.pbm");
    Plot plotV4 = initPlot(0.0, TF, -2, 2, x_desc, "Spannung an V4 (V)", "AM-Moduliertes Signal (Spannung an C2)", "v4.pbm");
    Plot plotV5 = initPlot(0.0, TF, -5.0, 5.0, x_desc, NULL, NULL, "v5.pbm");
    Plot plotV6 = initPlot(0.0, TF, -5.0, 5.0, x_desc, NULL, NULL, "v6.pbm");

    Plot plotV11 = initPlot(0.0, TF, -5, 5, x_desc, NULL, NULL, "v11.pbm");
    Plot plotV12 = initPlot(0.0, TF, -5, 5, x_desc, NULL, NULL, "v12.pbm");
    Plot plotV13 = initPlot(0.0, TF, -5, 5, x_desc, NULL, NULL, "v13.pbm");

    Plot plotI1 = initPlot(0.0, TF, -5 * MIKRO, 5 * MIKRO, x_desc, NULL, NULL, "i1.pbm");
    Plot plotI2 = initPlot(0.0, TF, -1 * MILLI, 1 * MILLI, x_desc, NULL, NULL, "i2.pbm");
    Plot plotI3u = initPlot(0.0, TF, -1 * MILLI, 1 * MILLI, x_desc, NULL, NULL, "i3u.pbm");
    Plot plotI4l = initPlot(0.0, TF, -20 * MILLI, 20 * MILLI, x_desc, NULL, NULL, "i4l.pbm");

    for(double t = 0; t < TF; t += DT) {
        setPoint(t, Vth(t), &plotV0);

        // ENFORCE ACTIVE / CUTOFF MODE

        I1 = (Vth(t) - VBE) / (Rth + (beta + 1) * R4);
        I2 = (beta + 1) * I1;
        I3u = beta * I1;

        if(I1 < 0) {
            I1 = I2 = I3u = 0;
        }

        setPoint(t, Vcc - (R5 * I3u), &plotV3b);

        V1_n1 = Vth(t) - Rth * I1;
        V2_n1 = 0 + R4 * I2;

        I4l += (V4_n / L1) * DT;
        dV4 = (1 / C2) * (-I3u + (Vcc - V3_n) / R5 - I4l) * DT;
        V4_n1 = V4_n + dV4;

        V3_n1 = V3_n + ((1 / C1) * (-I3u - (V3_n - Vcc) / R5) * DT + dV4);

        setPoint(t, V3_n1 - V2_n1, &plotV5);

        // CHECK MODE / ENFORCE SATURATION MODE

        if(V3_n1 - V2_n1 <= VBE && I1 > 0) {
            I3u = (V3_n - VCE) / R4 - (Vth(t) - (V3_n - VCB)) / Rth;

            I4l = I4l;
            dV4 = (1 / C2) * (-I3u + (V3_n - Vcc) / R5 - I4l) * DT;
            V4_n1 = V4_n + dV4;

            V3_n1 = V3_n + ((1 / C1) * (-I3u - (V3_n - Vcc) / R5) * DT + dV4);
            
            V1_n1 = V3_n1 - VCB;
            V2_n1 = V3_n1 - VCE;

            I1 = (Vth(t) - V1_n1) / Rth;
            I2 = (V2_n1 - 0) / R4;

        }

        V1_n = V1_n1;
        V2_n = V2_n1;
        V3_n = V3_n1;
        V4_n = V4_n1;

        setPoint(t, V1_n, &plotV1);
        setPoint(t, V2_n, &plotV2);
        setPoint(t, V3_n, &plotV3);
        setPoint(t, V4_n, &plotV4);

        setPoint(t, V3_n - V4_n, &plotV6);

        setPoint(t, I1, &plotI1);
        setPoint(t, I2, &plotI2);
        setPoint(t, I3u, &plotI3u);
        setPoint(t, I4l, &plotI4l);
    }

    savePlot(&plotV0);
    savePlot(&plotV1);
    savePlot(&plotV2);
    savePlot(&plotV3);
    savePlot(&plotV3b);
    savePlot(&plotV4);
    savePlot(&plotV5);
    savePlot(&plotV6);

    savePlot(&plotI1);
    savePlot(&plotI2);
    savePlot(&plotI3u);
    savePlot(&plotI4l);
}
