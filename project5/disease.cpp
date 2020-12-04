#include "disease.h"

void disease::initialize(double t_total, int S_init, int I_init, int R_init, double a, double b, double c, double d, double d_I, double e, double A, double w, double f_0, double fT)
{
  SIR(0) = S_init;
  SIR(1) = I_init;
  SIR(2) = R_init;
  SIR(3) = 0;
  SIR(4) = 0;
  SIR(5) = 400;
  time = 0;

  total = SIR(0) + SIR(1) + SIR(2);

  a_0 = a;
  this->a = a;
  this->b = b;
  this->c = c;
  this->d = d;
  this->d_I = d_I;
  this->e = e;
  this->A = A;
  this->w = w;
  this->f = 0;
  this->f_0 = f_0;
  this->vaccineTime = fT;


  dt = std::min({ 4 / (a * total), 1 / (b * total), 1 / (c * total), 1 / (d * total), 1 / (d_I * total), 1 / (f * total) });
  std::cout << "Time step: " << dt << std::endl;
  this->mcc = (int) (t_total / dt);

  // Initialize transition probabilities
  transitionProb();
}

void disease::transitionProb()
{

  total = SIR(0) + SIR(1) + SIR(2);

  // Probabiltiies of advancing in SIRS
  tProb(0, 1) = (a * SIR(0) * SIR(1) / total) * dt;
  tProb(1, 2) = (b * SIR(1)) * dt;
  tProb(2, 0) = (c * SIR(2)) * dt;

  // Probabilities of death by unrelated causes
  tProb(0, 3) = (d * SIR(0)) * dt;
  tProb(1, 3) = (d * SIR(1)) * dt;
  tProb(2, 3) = (d * SIR(2)) * dt;

  // Probability of death by infection
  tProb(1, 4) = (d_I * SIR(1)) * dt;

  // Probability of birth
  tProb(5, 0) = (e * total) * dt;

  // Probability of vaccination
  tProb(0, 2) = (f * SIR(0)) * dt;

  //tProb.print();
}

void disease::transmissionRate(double time)
{
  a = A * std::cos(w * time) + a_0;
}

void disease::vaccinationRate(double time)
{
  if (time > vaccineTime)
  {
    f = std::min({ f_0, f_0 * 0.1 * (time - vaccineTime)});
  }
}

void disease::monteCarlo()
{
  // Initialize rng
  std::random_device rd;
  std::mt19937_64 gen(rd());

  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

  for (int cycles = 1; cycles <= mcc; cycles++)
  {
    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        if (RandomNumberGenerator(gen) < tProb(i, j))
        {
          SIR(i)--;
          SIR(j)++;
          transitionProb();
        }
      }
    }
    output(SIR, out_val);
    EXP += SIR;
    output(EXP/cycles, out_exp);

    time += dt;
    transmissionRate(time);
    vaccinationRate(time);
    //std::cout << "Cycle #" << cycles << " completed" << std::endl;
  }
}

void disease::output(arma::vec sir, std::ofstream& out)
{
  for (int n = 0; n <  sir.n_elem; n++)
  {
    out << std::setw(25) << std::setprecision(8) << sir(n);
  }
  out << std::endl;
}

void disease::solver(std::string filename)
{
  // Prepare files
  std::string filename_val = filename + "_val.txt";
  std::string filename_exp = filename + "_exp.txt";

  out_val.open(filename_val);
  out_exp.open(filename_exp);

  monteCarlo();

  out_val.close();
  out_exp.close();

  std::cout << mcc << " Cycles Total." << std::endl;

}
