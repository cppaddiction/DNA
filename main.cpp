#include <iostream>
#include <fstream>
#include <vector>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef vector<complex<double>> state_type;

const complex<double> I(0.0, 1.0);
const complex<double> realI(1.0, 0.0);
const complex<double> zero(0.0, 0.0);

complex<double> x0 = complex<double>(4.0, 0.0);
const complex<double> omega_0_squared = complex<double>(0.0001, 0.0);
const complex<double> omega_0_shtrix = complex<double>(0.006, 0.0);
const complex<double> nu0 = complex<double>(124.0 / 6.582119514, 0.0);
const complex<double> nu01 = complex<double>(1.92, 0.0);

complex<double> x1 = complex<double>(4.0, 0.0);
const complex<double> omega_1_squared = complex<double>(0.0001, 0.0);
const complex<double> omega_1_shtrix = complex<double>(0.006, 0.0);
const complex<double> nu10 = complex<double>(1.92, 0.0);
const complex<double> nu1 = complex<double>(169.0 / 6.582119514, 0.0);
const complex<double> nu12 = complex<double>(1.92, 0.0);

complex<double> x2 = complex<double>(4.0, 0.0);
const complex<double> omega_2_squared = complex<double>(0.0001, 0.0);
const complex<double> omega_2_shtrix = complex<double>(0.006, 0.0);
const complex<double> nu21 = complex<double>(1.92, 0.0);
const complex<double> nu2 = complex<double>(190.0 / 6.582119514, 0.0);
const complex<double> nu23 = complex<double>(1.92, 0.0);

complex<double> x3 = complex<double>(4.0, 0.0);
const complex<double> omega_3_squared = complex<double>(0.0001, 0.0);
const complex<double> omega_3_shtrix = complex<double>(0.006, 0.0);
const complex<double> nu32 = complex<double>(1.92, 0.0);
const complex<double> nu3 = complex<double>(124.0 / 6.582119514, 0.0);
const complex<double> nu34 = complex<double>(1.92, 0.0);

complex<double> x4 = complex<double>(4.0, 0.0);
const complex<double> omega_4_squared = complex<double>(0.0001, 0.0);
const complex<double> omega_4_shtrix = complex<double>(0.006, 0.0);
const complex<double> nu43 = complex<double>(1.92, 0.0);
const complex<double> nu4 = complex<double>(190.0 / 6.582119514, 0.0);
const complex<double> nu45 = complex<double>(1.92, 0.0);

complex<double> x5 = complex<double>(4.0, 0.0);
const complex<double> omega_5_squared = complex<double>(0.0001, 0.0);
const complex<double> omega_5_shtrix = complex<double>(0.006, 0.0);
const complex<double> nu54 = complex<double>(1.92, 0.0);
const complex<double> nu5 = complex<double>(124.0 / 6.582119514, 0.0);
const complex<double> nu56 = complex<double>(1.92, 0.0);

complex<double> x6 = complex<double>(4.0, 0.0);
const complex<double> omega_6_squared = complex<double>(0.0001, 0.0);
const complex<double> omega_6_shtrix = complex<double>(0.006, 0.0);
const complex<double> nu65 = complex<double>(1.92, 0.0);
const complex<double> nu6 = complex<double>(124.0 / 6.582119514, 0.0);
const complex<double> nu67 = complex<double>(1.92, 0.0);

complex<double> x7 = complex<double>(4.0, 0.0);
const complex<double> omega_7_squared = complex<double>(0.0001, 0.0);
const complex<double> omega_7_shtrix = complex<double>(0.006, 0.0);
const complex<double> nu76 = complex<double>(1.92, 0.0);
const complex<double> nu7 = complex<double>(124.0 / 6.582119514, 0.0);

struct diff_system
{
    template<class State>
    void operator()(const State& x, State& dxdt, double t)
    {
        dxdt[0] = (-I) * ((nu0 + x0 * omega_0_squared * x[8]) * x[0] + nu01 * x[1]);
        dxdt[1] = (-I) * ((nu1 + x1 * omega_1_squared * x[10]) * x[1] + nu10 * x[0] + nu12 * x[2]);
        dxdt[2] = (-I) * ((nu2 + x2 * omega_2_squared * x[12]) * x[2] + nu21 * x[1] + nu23 * x[3]);
        dxdt[3] = (-I) * ((nu3 + x3 * omega_3_squared * x[14]) * x[3] + nu32 * x[2] + nu34 * x[4]);
        dxdt[4] = (-I) * ((nu4 + x4 * omega_4_squared * x[16]) * x[4] + nu43 * x[3] + nu45 * x[5]);
        dxdt[5] = (-I) * ((nu5 + x5 * omega_5_squared * x[18]) * x[5] + nu54 * x[4] + nu56 * x[6]);
        dxdt[6] = (-I) * ((nu6 + x6 * omega_6_squared * x[20]) * x[6] + nu65 * x[5] + nu67 * x[7]);
        dxdt[7] = (-I) * ((nu7 + x7 * omega_7_squared * x[22]) * x[7] + nu76 * x[6]);
        dxdt[8] = x[9];
        dxdt[9] = -omega_0_shtrix * x[9] - omega_0_squared * x[8] - norm(x[0]);
        dxdt[10] = x[11];
        dxdt[11] = -omega_1_shtrix * x[11] - omega_1_squared * x[10] - norm(x[1]);
        dxdt[12] = x[13];
        dxdt[13] = -omega_2_shtrix * x[13] - omega_2_squared * x[12] - norm(x[2]);
        dxdt[14] = x[15];
        dxdt[15] = -omega_3_shtrix * x[15] - omega_3_squared * x[14] - norm(x[3]);
        dxdt[16] = x[17];
        dxdt[17] = -omega_4_shtrix * x[17] - omega_4_squared * x[16] - norm(x[4]);
        dxdt[18] = x[19];
        dxdt[19] = -omega_5_shtrix * x[19] - omega_5_squared * x[18] - norm(x[5]);
        dxdt[20] = x[21];
        dxdt[21] = -omega_6_shtrix * x[21] - omega_6_squared * x[20] - norm(x[6]);
        dxdt[22] = x[23];
        dxdt[23] = -omega_7_shtrix * x[23] - omega_7_squared * x[22] - norm(x[7]);
    }
};

struct streaming_observer
{
    ostream& m_out;

    streaming_observer(ostream& out) : m_out(out) { }

    template<class State>
    void operator()(const State& x, double t) const
    {
        m_out << t * pow(10, -14) << " ";
        for (int i = 0; i < 8; i++)
        {
            m_out << norm(x[i]) << " ";
        }
        m_out << "\n";
    }
};

void ChooseKappa()
{
    double kappa_span = 26.0;  // choosing kappa from 4.00 to 30.00

    for (int i = 0; i <= 2600; i++)
    {
        double kappa = kappa_span * i / 2600 + 4;

        state_type x(24, zero);
        x[0] = realI;

        x0 = complex<double>(kappa, 0.0);
        x1 = complex<double>(kappa, 0.0);
        x2 = complex<double>(kappa, 0.0);
        x3 = complex<double>(kappa, 0.0);
        x4 = complex<double>(kappa, 0.0);
        x5 = complex<double>(kappa, 0.0);
        x6 = complex<double>(kappa, 0.0);
        x7 = complex<double>(kappa, 0.0);

        auto stepper = make_dense_output<runge_kutta_dopri5<state_type>>(1.0e-8, 1.0e-8);
        stepper.initialize(x, 0.0, 0.01);

        vector<double> values;

        while (true)
        {
            auto time_interval = stepper.do_step(diff_system());
            stepper.calc_state((time_interval.first + time_interval.second) / 2.0, x);

            if (time_interval.second >= 2000 && time_interval.second <= 2100)
            {
                values.push_back(norm(x[0]));
            }
            else if (time_interval.second > 2100)
            {
                double avg = accumulate(values.begin(), values.end(), 0.0) / values.size();
                cout << kappa << " " << avg << "\n";

                if (avg < 0.03)
                {
                    return;
                }
                else
                {
                    break;
                }
            }
        }
    }
}

void PerformCalculations(ostream& out)
{
    state_type x(24, zero);
    x[0] = realI;

    size_t num_of_steps = integrate_const(make_dense_output<runge_kutta_dopri5<state_type>>(1.0e-8, 1.0e-8),
        diff_system(), x, 0.0, 5000.0, 0.01,
        streaming_observer(out));

    cout << num_of_steps << "\n";
}

int main()
{
    ofstream out("output.txt");

    ChooseKappa();

    PerformCalculations(out);

    out.close();

    return 0;
}