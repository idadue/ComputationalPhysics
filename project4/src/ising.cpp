#include "ising.h"

//initialize energy, magnetization
void ising::initialize(int n_spins, arma::mat & spin_mat, double& E, double& M, int random_start)
{
  // Initialize rng
  std::random_device rd;
	std::mt19937_64 gen(rd());

	std::uniform_int_distribution<int> RandomInt(0,1);

  //setup spin matrix
  for (arma::sword i = 0; i < n_spins; i++){
    for (arma::sword j = 0; j < n_spins; j++){
      spin_mat(i,j) = 1 - 2 * random_start * RandomInt(gen);
      M += (float) spin_mat(i,j);
    }
  }
  E = compute_energy(spin_mat);
}

// Calculates the energy of a lattice
double ising::compute_energy(arma::mat spin_mat)
{
  double energy = 0;
  arma::sword n_spins = spin_mat.n_rows;
  for (arma::sword i = 0; i < n_spins; i++){
      for (arma::sword j = 0; j < n_spins; j++){
        energy -= spin_mat(i,j) * spin_mat(pbc(i + 1, n_spins), j) + spin_mat(i, j) * spin_mat(i, pbc(j + 1, n_spins));
      }
    }
  return energy;
}

// implement periodic boundary conditions
arma::sword ising::pbc(arma::sword coordinate, arma::sword max_coordinate)
{
  if (coordinate >= max_coordinate)
  {
    coordinate = coordinate - max_coordinate;
  }
  else if (coordinate < 0)
  {
    coordinate = coordinate + max_coordinate;
  }
  return coordinate;
}

// Returns the transition probability for positive energy differences
arma::vec ising::transitionProb(double temperature)
{
  int n_ediff = 5;
  arma::vec probability = arma::zeros<arma::vec>(n_ediff);
  for (arma::sword i = 0; i < 5; i++)
  {
    //std::cout << -4*(i-2.0) << std::endl;
    probability(i) = std::exp(-4*(i-2.0)/temperature);
  }
  return probability;
}

void ising::mc_temp(int n_spins, int mcc, double temperature_start, double temperature_end, int n_temperature, std::string filename, int random_start=0)
{
  // Prepare files
  std::string filename_energy = filename + "_energy.txt";
  std::string filename_temp = filename + "_temp.txt";
  out_energy.open(filename_energy);
  out_temp.open(filename_temp);

  // Coung configurations
  arma::vec config_acceptance = arma::zeros<arma::vec>(n_temperature);

  double dT = (temperature_end - temperature_start) + 2;

  if (n_temperature > 1)
  {
    dT = (temperature_end - temperature_start)/(n_temperature - 1);
  }

  int i = 0;
  arma::vec thermo_temp = arma::zeros<arma::vec>(4);

  for (double temp = temperature_start; temp <= temperature_end; temp += dT)
  {
    std::cout << "T = " << temp << std::endl;
    std::string filename_thermo = filename + "_thermo_" + std::to_string(i) + ".txt";
    out_thermo.open(filename_thermo);
    // Metropolis algorithm to find thermo values
    thermo_temp = metropolis(n_spins, mcc, temp, filename, random_start);
    // Save final values for temperature plot
    output(thermo_temp, out_temp);

    config_acceptance(i) = flip_count;

    out_thermo.close();
    i++;
  }

  output(config_acceptance, out_mcc);

  out_energy.close();
  out_temp.close();

}

void ising::mc_mcc(int n_spins, int mcc, double temperature_start, double temperature_end, int n_temperature, std::string filename, int random_start=0)
{
  // Prepare files
  std::string filename_mcc = filename + "_mcc.txt";
  out_mcc.open(filename_mcc);

  int i = 0;

  for (int mcc_loop = 10; mcc_loop <= mcc; mcc_loop *= 10)
  {
    // mc_temp to find accepted configs as function of temperature
    std::cout << "MCC = " << mcc_loop << std::endl;
    mc_temp(n_spins, mcc_loop, temperature_start, temperature_end, n_temperature, filename, random_start);
  }

  out_mcc.close();

}


arma::vec ising::metropolis(int n_spins, int mcc, double temperature, std::string filename, int random_start=0)
{
  //initializing values
  double energy = 0.;
  double dEnergy;
  double mag_mom = 0.;
  arma::vec exp = arma::zeros<arma::vec>(5);
  arma::mat spin_mat = arma::zeros<arma::mat>(n_spins, n_spins);
  arma::vec energy_level_count = arma::zeros<arma::vec>(n_spins*n_spins + 1);
  bool equilibrium = true;
  int equilibrium_cycles = 0;
  int equilibrium_time = 0;
  flip_count = 0;
  arma::vec thermo_values = arma::zeros<arma::vec>(4);

  // Initialize rng
  std::random_device rd;
  std::mt19937_64 gen(rd());

  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

  // Initialize spin matrix
  initialize(n_spins, spin_mat, energy, mag_mom, random_start);

  // Printing initial state pre-mcc
  arma::vec init = arma::zeros<arma::vec>(5);
  init(0) = energy;  init(1) = energy*energy;  init(2) = mag_mom;  init(3) = mag_mom*mag_mom;  init(4) = std::abs(mag_mom);
  thermo_values = data_conversion(n_spins, 1, init, temperature);
  output(thermo_values, out_thermo);

  // Set up transition probability vector
  arma::vec probVec = transitionProb(temperature);

  for (int cycles = 1; cycles <= mcc; cycles++)
  {
    //std::cout << "Cycle #" << cycles << std::endl;
    for (arma::uword i = 0; i < n_spins; i++)
    {
      for (arma::uword j = 0; j < n_spins; j++)
      {
        //std::cout << "Coordinates " << i << ", " << j << std::endl;
        arma::sword i_flip = (arma::sword) (RandomNumberGenerator(gen)*(double)n_spins);
        arma::sword j_flip = (arma::sword) (RandomNumberGenerator(gen)*(double)n_spins);

        //std::cout << i_flip << ", " << j_flip << std::endl;

        dEnergy = 2 * spin_mat(i_flip, j_flip) * (spin_mat(pbc(i_flip + 1, n_spins), j_flip) + spin_mat(pbc(i_flip - 1, n_spins), j_flip) + spin_mat(i_flip, pbc(j_flip + 1, n_spins)) + spin_mat(i_flip, pbc(j_flip - 1, n_spins)));

        arma::sword energyCoordinate = (arma::sword) dEnergy/4 + 2;
        //std::cout << "denergy = " << dEnergy << std::endl;
        //std::cout << "energy coordinate = " << energyCoordinate << std::endl;
        //std::cout << "transitionProb = " << probVec(energyCoordinate) << std::endl;

        if (RandomNumberGenerator(gen) < probVec(energyCoordinate))
        {
          flip_count += 1;
          spin_mat(i_flip, j_flip) *= -1;
          energy += dEnergy;
          mag_mom += 2*spin_mat(i_flip, j_flip);
          //std::cout << "Flipped " << i_flip << ", " << j_flip << std::endl;
        }
      }
    }

    if (cycles > mcc*0.0)
    {
      exp(0) += energy;
      exp(1) += energy*energy;
      exp(2) += mag_mom;
      exp(3) += mag_mom*mag_mom;
      exp(4) += std::abs(mag_mom);

      thermo_values = data_conversion(n_spins, cycles, exp, temperature);
      output(thermo_values, out_thermo);
    }

    equilibrium_time += 1;

    //equilibrium = is_equilibrium();

    if (equilibrium)
    {
      arma::uword energy_coordinate = (arma::uword) ((compute_energy(spin_mat) + 2 * n_spins * n_spins)/4);
      energy_level_count(energy_coordinate) += 1.0;
      equilibrium_cycles += 1;
    }
  }
  //exp /= (mcc*1.0);

  if (equilibrium == true)
  {
    output(energy_level_count, out_energy);
    //std::cout << "Energy Distribution: " << std::endl;
    energy_level_count /= equilibrium_cycles;
    /*for (int i = 0; i < n_spins*n_spins + 1; i++)
    {
      std::cout << std::setw(5) << i * 4 - n_spins * n_spins * 2 << " J: " << energy_level_count(i) << std::endl;
    }*/
  }
  std::cout << (float) (100.0*flip_count/(mcc*n_spins*n_spins)) << "% accepted configurations." << std::endl;
  std::cout << "E, |M|, C_V, X" << std::endl;
  thermo_values.print();
  std::cout << "Fission Mailed Successfully" << std::endl;
  return thermo_values;
}

bool is_equilibrium()
{
  bool equilibrium = true;

  return equilibrium;
}

arma::vec ising::data_conversion(int n_spins, int cycles, arma::vec exp, double temperature)
{
  exp /= cycles;
  double exp_energy = exp(0);
  double var_energy = (exp(1) - exp(0)*exp(0));
  double exp_magmom = exp(2);
  double var_magmom = (exp(3) - exp(2)*exp(2));
  double exp_magmom_abs = exp(4);

  double spec_heat = var_energy / temperature / temperature;
  double mag_suscept = var_magmom / temperature;

  double L2 = n_spins * n_spins;

  arma::vec thermo_values = arma::zeros<arma::vec>(4);
  thermo_values(0) = exp_energy/L2;
  thermo_values(1) = exp_magmom_abs/L2;
  thermo_values(2) = spec_heat/L2;
  thermo_values(3) = mag_suscept/L2;

  return thermo_values;
}

void ising::output(arma::vec thermo_values, std::ofstream& out)
{
  for (int n = 0; n <  thermo_values.n_elem; n++)
  {
    out << std::setw(15) << std::setprecision(8) << thermo_values(n);
  }
  out << std::endl;

  /*
  out << std::fixed << std::setprecision(15) << exp_energy/L2;
  out << std::fixed << std::setprecision(15) << exp_magmom_abs/L2;
  out << std::fixed << std::setprecision(15) << spec_heat/L2;
  out << std::fixed << std::setprecision(15) << mag_suscept/L2 << std::endl;*/
  /*
  std::cout << "E(E)/L² = " << exp_energy/L2 << ", Var(E)/L² = " << var_energy/L2 << std::endl;
  std::cout << "E(M)/L² = " << exp_magmom/L2 << ", Var(M)/L² = " << var_magmom/L2 << std::endl;
  std::cout << "C_V/L² = " << spec_heat/L2 << ", X/L² = " << mag_suscept/L2 << std::endl;
  std::cout << "E(|M|)/L² = " << exp_magmom_abs/L2 << std::endl;
  */
  /*
  std::cout << "E(E) = " << exp_energy << ", Var(E) = " << var_energy << std::endl;
  std::cout << "E(M) = " << exp_magmom << ", Var(M) = " << var_magmom << std::endl;
  std::cout << "C_V = " << spec_heat << ", X = " << mag_suscept << std::endl;
  */
}
