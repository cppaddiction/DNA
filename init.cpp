#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

int main()
{
    ofstream out("output.txt");
    ostringstream equations;

    string chain; cin >> chain;

    for (int i = 0; i < chain.size(); i++)
    {
        out << "complex<double> x" << i << " = complex<double>(4.0, 0.0);" << "\n";
        out << "const complex<double> omega_" << i << "_squared" << " = complex<double>(0.0001, 0.0);" << "\n";
        out << "const complex<double> omega_" << i << "_shtrix" << " = complex<double>(0.006, 0.0);" << "\n";
        if (i == 0)
        {
            out << "const complex<double> nu" << i << " = complex<double>(" << (chain[i] == 'G' ? "124.0/6.582119514, 0.0);" : (chain[i] == 'A' ? "169.0/6.582119514, 0.0);" : "190.0/6.582119514, 0.0);")) << "\n";
            out << "const complex<double> nu" << i << i + 1 << " = complex<double>(1.92, 0.0);" << "\n\n";
            equations << "dxdt[" << i << "] = (-I) * ((nu" << i << " + x" << i << " * omega_" << i << "_squared" << " * x[" << chain.size() + 2 * i << "]) * x[" << i << "] + nu" << i << i + 1 << " * x[" << i + 1 << "]);" << "\n";
        }
        else if (i == chain.size() - 1)
        {
            out << "const complex<double> nu" << i << i - 1 << " = complex<double>(1.92, 0.0);" << "\n";
            out << "const complex<double> nu" << i << " = complex<double>(" << (chain[i] == 'G' ? "124.0/6.582119514, 0.0);" : (chain[i] == 'A' ? "169.0/6.582119514, 0.0);" : "190.0/6.582119514, 0.0);")) << "\n\n";
            equations << "dxdt[" << i << "] = (-I) * ((nu" << i << " + x" << i << " * omega_" << i << "_squared" << " * x[" << chain.size() + 2 * i << "]) * x[" << i << "] + nu" << i << i - 1 << " * x[" << i - 1 << "]);" << "\n";
        }
        else
        {
            out << "const complex<double> nu" << i << i - 1 << " = complex<double>(1.92, 0.0);" << "\n";
            out << "const complex<double> nu" << i << " = complex<double>(" << (chain[i] == 'G' ? "124.0/6.582119514, 0.0);" : (chain[i] == 'A' ? "169.0/6.582119514, 0.0);" : "190.0/6.582119514, 0.0);")) << "\n";
            out << "const complex<double> nu" << i << i + 1 << " = complex<double>(1.92, 0.0);" << "\n\n";
            equations << "dxdt[" << i << "] = (-I) * ((nu" << i << " + x" << i << " * omega_" << i << "_squared" << " * x[" << chain.size() + 2 * i << "]) * x[" << i << "] + nu" << i << i - 1 << " * x[" << i - 1 << "]" << " + nu" << i << i + 1 << " * x[" << i + 1 << "]);" << "\n";
        }
    }

    for (int i = chain.size(), k = 0; i < 3 * chain.size(); i += 2, k++)
    {
        equations << "dxdt[" << i << "] = x[" << i + 1 << "];" << "\n";
        equations << "dxdt[" << i + 1 << "] = -omega_" << k << "_shtrix * x[" << i + 1 << "] - omega_" << k << "_squared" << " * x[" << i << "] - norm(x[" << k << "]);" << "\n";
    }

    out << equations.str();
    out.close();

    return 0;
}