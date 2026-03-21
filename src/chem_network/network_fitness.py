import numpy as np
from scipy.stats import truncnorm
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import pytest

#
#  define network fitness class
#

class network_fitness():
    def __init__(self, mol_set, mol_mass, distr_p, max_fitness):
        self.nmol = len(mol_set)
        self.molecule_indices = np.array(mol_set)
        self.molecule_mass = np.array(mol_mass)
        self.molec_fitness_distr = None
        self.molecules_fitness = None
        self.tot_mass = 0.
        self.set_molecules_fitness(distr_p, max_fitness)
    def set_molecules_fitness(self, distr_p, max_fitness=1.0):
        mean = distr_p[0]
        std_dev = distr_p[1]
        # Calculate the truncation points for the standard normal distribution
        # (a and b are in terms of standard deviations from the mean)
        min_val = 0
        max_val = self.nmol
        a = (min_val - mean) / std_dev
        b = (max_val - mean) / std_dev
        # Create a truncated normal distribution object
        truncated_dist = truncnorm(a=a, b=b, loc=mean, scale=std_dev)
        # Evaluate PDF at each molecule index (0 to nmol - 1)
        raw_fitness = truncated_dist.pdf(self.molecule_indices)
        # Check if any values > 1, and normalize only if needed
        if np.any(raw_fitness > max_fitness):
            normalized_fitness = raw_fitness / np.max(raw_fitness) * max_fitness
        else:
            normalized_fitness = raw_fitness  # already safe
        # Store the fitness values
        self.molec_fitness_distr = truncated_dist
        self.molecules_fitness = normalized_fitness
        print(self.molecules_fitness)
    def show_fitness_distr(self):
        # 2. Define the range of x-values for plotting
        # It's good practice to create a range that covers most of the distribution's mass.
        # For a normal distribution, the mean ± 4 standard deviations is a common choice.
        x_values = np.linspace(self.molec_fitness_distr.ppf(0.001), self.molec_fitness_distr.ppf(0.999), 100)
        # 3. Calculate the PDF for each x-value using the frozen object
        pdf_values = self.molec_fitness_distr.pdf(x_values)
        # 4. Plot the results
        plt.figure(figsize=(10, 6))
        plt.plot(x_values, pdf_values, label='Frozen Normal PDF', color='blue')
        plt.title('PDF of Molecular Fitness Distribution')
        plt.xlabel('X-values')
        plt.ylabel('Fitness')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.show()
        plt.close()
    def check_total_mass_conserved(self, states_t):
        tot_mass0 = 0.
        for mol_id in self.molecule_indices:
            mass = self.molecule_mass[mol_id-1]
            tot_mass0 += mass * states_t[0, mol_id-1]
        tot_mass = 0.
        for mol_id in self.molecule_indices:
            mass = self.molecule_mass[mol_id-1]
            tot_mass += mass * states_t[-1, mol_id-1]
        assert tot_mass == pytest.approx(tot_mass0, 1.e-4)
        self.tot_mass = tot_mass
    def compute_mass_score(self, states_t):
        """
        Compute fitness based on mass for target molecules.
        """
        total_score = 0.
        for mol_id in self.molecule_indices:
            nmol_profile = states_t[:, mol_id-1]
            final_nmol = nmol_profile[-1]
            total_score += final_nmol * self.molecule_mass[mol_id-1] * self.molecules_fitness[mol_id-1] / self.tot_mass
        self.mass_score = total_score
    def compute_avg_react_times(self, t_grid, states_t, plot=True):
        """
        Compute fitness based on speed for target molecules.
        use time to treshold
        s_i = max(d n_mol(i) / dt)
        if t_gr[-1] < eps:
            nothing happens -> s_i = 0
        else:
            s_i = max(d n_mol(i) / dt)
        s_i^n = s_i / max(s_i)
        """
        self.avg_times = {}
        print('len: ', len(t_grid))
        if len(t_grid) > 1:
            # compute max slope each molecule
            for mol_id in self.molecule_indices:
                nmol_profile = states_t[:, mol_id-1]
                # inst. derivative
                dnmol_dt = np.gradient(savgol_filter(nmol_profile, 5, 3), t_grid)
                pos_slope = np.maximum(dnmol_dt, 0)
                # weighted time : sum(t * dN/dt) / sum(dN/dt)
                total_production = np.sum(pos_slope)
                if total_production > 0.:
                    t_avg = np.sum(pos_slope * t_grid) / total_production
                else:
                    t_avg = np.inf
                self.avg_times[mol_id] = t_avg
                if plot:
                    plt.figure(figsize=(8,4))
                    plt.plot(t_grid, nmol_profile, label='Profile')
                    plt.plot(t_grid, pos_slope, label='Positive derivative')
                    plt.axvline(t_avg, color='r', linestyle='--', label=f'Weighted avg t = {t_avg:.2f}')
                    plt.xlabel('Time')
                    plt.ylabel('Molecule count / Production rate')
                    plt.title(f'Molecule {mol_id} weighted avg reaction time')
                    plt.legend()
                    plt.grid(True, alpha=0.5, linestyle='--')
                    plt.show()
        else:
            for mol_id in self.molecule_indices:
                self.avg_times[mol_id] = np.inf
    def weight_react_times(self):
        """
        Compute time-based weights for molecules based on average reaction times.

        Parameters
        ----------
        mode : str
            'inverse'  -> weight = 1 / t_avg (fast = high weight)
            'linear'   -> weight = 1 - normalized(t_avg)
        normalize : bool
        Whether to normalize weights to [0, 1] or sum to 1.
        """
        mol_ids = list(self.avg_times.keys())
        t_values = np.array(list(self.avg_times.values()), dtype=float)
        print(mol_ids, t_values)
        # Replace NaNs with infs (in case something weird slipped in)
        t_values[np.isnan(t_values)] = np.inf
        # Detect finite entries
        finite_mask = np.isfinite(t_values)
        # Case 1: all molecules have infinite times → no reactions occurred
        if not np.any(finite_mask):
            self.reaction_time_weights = {mid: 0.0 for mid in mol_ids}
        else:
            finite_t = t_values[finite_mask]
            inv_t = 1.0 / finite_t
            norm_inv_t = inv_t / np.sum(inv_t)
            # store weights
            weights = np.zeros_like(t_values)
            weights[finite_mask] = norm_inv_t
            self.reaction_time_weights = {mid: w for mid, w in zip(mol_ids, weights)}
        print(self.reaction_time_weights)
    def compute_react_t_score(self):
        # Compute weighted time score
        self.weighted_time_score = sum(
            self.molecules_fitness[mid - 1] * w for mid, w in self.reaction_time_weights.items()
        )
        print(self.weighted_time_score)
    def set_fitness(self, t_grid, states_t, fitness_p):
        self.check_total_mass_conserved(states_t)
        self.compute_mass_score(states_t)
        # reactions times
        self.compute_avg_react_times(t_grid, states_t)
        self.weight_react_times()
        self.compute_react_t_score()
        fitness = fitness_p[0] * self.mass_score + fitness_p[1] * self.weighted_time_score
        return fitness

'''

def compute_fitness(self, kinetic_solver):
    """
    Compute fitness based on speed and amount of production for target molecules.
    """
    target_ids = p.target_molecules  # indices of target molecules
    t_grid = kinetic_solver.t_grid
    avg_states = kinetic_solver.avg_states_t  # shape: [time, species]

    total_fitness = 0.0

    for mol_id in target_ids:
        conc_profile = avg_states[:, mol_id]  # concentration over time
        final_conc = conc_profile[-1]         # amount at end
        max_conc = max(conc_profile)

        # --- Speed: time to reach 50% of max (or final) ---
        half_max = 0.5 * max_conc
        time_to_half = next(
            (t for t, c in zip(t_grid, conc_profile) if c >= half_max),
            t_grid[-1]  # fallback if never reaches
        )

        # Normalize time-to-half (lower is better)
        speed_score = 1.0 / (time_to_half + 1e-6)  # avoid div by 0

        # Normalize final amount (assuming max is 1.0 for scaling)
        amount_score = final_conc

        # Combine (weights can be tuned)
        fitness = 0.7 * amount_score + 0.3 * speed_score
        total_fitness += fitness

    # Store the fitness value
    self.fitness = total_fitness


'''
