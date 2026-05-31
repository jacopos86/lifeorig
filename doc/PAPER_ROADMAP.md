# Paper Roadmap

## Paper 1: Coupled Atmosphere-Ocean Chemistry Solver

Build and validate the planetary environment engine.

Core contribution:

- hydrostatic layered atmosphere;
- wavelength-dependent radiative-convective temperature solver;
- equilibrium and disequilibrium atmospheric chemistry;
- coupling to ocean chemistry through gas dissolution, surface fluxes, pH, and feedstock delivery.

Main validation:

- hydrostatic and radiative-convective benchmarks;
- EasyChem equilibrium benchmarks;
- simple atmospheric escape / photochemistry checks;
- ocean chemistry benchmarks such as Henry-law dissolution and carbonate/pH behavior.

Scientific angle:

- a predictive boundary-condition model for prebiotic environments from planet/star inputs.

## Paper 2: Earth-Like Volcanic Rock Environments

Use the coupled atmosphere-ocean/rock solver to study early Earth-like volcanic settings.

Core question:

- which atmospheric states and stellar UV conditions generate useful surface chemistry in volcanic rock pore networks?

Key outputs:

- surface temperature and pressure;
- dissolved atmospheric species;
- pH and redox conditions;
- wet-dry or liquid-level forcing;
- prebiotic feedstock availability near mineral surfaces.

## Paper 3: Ocean-Like Hydrothermal Vent Worlds

Extend the environment model to hydrothermal vent settings.

Core question:

- how do atmospheric composition, ocean chemistry, and vent gradients jointly constrain origin-of-life chemistry?

Key outputs:

- dissolved gas inventory;
- pH, redox, and temperature gradients;
- mineral-interface chemistry;
- availability of carbon, nitrogen, sulfur, and hydrogen feedstocks.

## Paper 4: Titan-Like Environments

Develop the Titan branch as a distinct low-temperature organic chemistry environment.

Core question:

- how do known Titan-like atmospheric and surface compositions constrain prebiotic chemistry in hydrocarbon-rich environments?

Key outputs:

- prescribed or observationally anchored molecular abundances;
- surface liquid composition;
- organic feedstock fluxes;
- comparison with Earth-like aqueous pathways.

## Paper 5: Probabilistic Origin-of-Life Study Across Planet Classes

Use the validated environment modules to run large parameter studies.

Planet classes:

- Earth-like terrestrial planets;
- Super-Earths;
- mini-Neptunes / atmosphere-loss transition objects;
- Titan-like worlds;
- hydrothermal ocean worlds.

Core question:

- which planet/star/environment combinations maximize plausible origin-of-life chemical pathways?

Method:

- sample planetary parameters, stellar spectra, volatile inventories, atmosphere loss, and surface environments;
- propagate uncertainties through atmosphere, ocean, and chemical network modules;
- produce probability maps for prebiotic feedstock availability and reaction-network viability.

## Paper 6: Full Biology / Bioinformatics-Inspired Origin Model

Add a biological emergence layer on top of the planetary chemistry model.

Core contribution:

- connect prebiotic chemical networks to protocell-like dynamics;
- use bioinformatics-inspired methods to compare emergent chemical networks with motifs from extant metabolism;
- quantify transitions from geochemistry to evolvable biochemical organization.

Core question:

- which planetary environments produce chemical networks closest to plausible ancestral biochemical structure?

Long-term synthesis:

- planet/star physics -> atmosphere/ocean chemistry -> prebiotic networks -> protocell dynamics -> biological emergence metrics.

## Additional Follow-Up Directions

### Atmospheric Escape and Secondary Atmosphere Formation

Study how H/He-rich primordial atmospheres evolve into secondary terrestrial atmospheres.

Core question:

- how do mini-Neptune or sub-Neptune atmospheres transition toward terrestrial CO2/N2/H2O atmospheres after envelope loss?

Key outputs:

- H/He loss timescales;
- volatile retention;
- surface pressure and temperature after escape;
- prebiotic surface conditions after secondary atmosphere formation.

### Stellar Type and Prebiotic Chemistry

Compare M dwarfs, K dwarfs, G stars, and active young stars.

Core question:

- which stellar spectra maximize useful photochemistry without sterilizing surface environments?

Key outputs:

- UV-driven HCN production;
- CH4/NH3 photolysis;
- atmospheric shielding;
- surface radiation dose.

### Volatile Inventory Phase Diagram

Map atmospheric and prebiotic outcomes across volatile parameter space.

Possible axes:

- C/O ratio;
- N/O ratio;
- H inventory;
- planet mass;
- stellar flux;
- surface pressure.

Key outputs:

- regions producing CO2/N2, H2/CH4, HCN-rich, ocean-compatible, or chemically poor atmospheres.

### Ocean-Atmosphere Redox Evolution

Track coupled atmospheric and ocean redox evolution.

Core question:

- when do atmospheres remain reducing enough for prebiotic synthesis?

Key processes:

- H2 escape;
- Fe2+/Fe3+ ocean buffering;
- volcanic gas redox state;
- CO/CO2/CH4 balance.

### Prebiotic Feedstock Delivery Fluxes

Focus on deposition and delivery rates rather than only concentrations.

Target species:

- HCN;
- formaldehyde;
- CO;
- CH4;
- NH3;
- nitrates and nitrites;
- sulfur species.

Key outputs:

- surface/ocean deposition rates under different atmospheres and stellar spectra.

### Open-Source Methods and Benchmark Paper

Publish the solver as a validated software framework.

Requirements:

- documentation;
- tests;
- benchmark cases;
- example notebooks;
- comparison with known atmosphere, chemistry, and ocean models.

Purpose:

- make later science papers easier to reproduce and cite.

### Impact-Driven Atmosphere Chemistry

Model transient post-impact atmospheric states.

Core question:

- how do impacts create reducing atmospheres, high temperatures, vapor plumes, and short-lived prebiotic feedstock pulses?

Key outputs:

- post-impact atmospheric composition;
- HCN/feedstock production;
- cooling and chemical relaxation timescales.

### Planetary Classifier for Origin-of-Life Potential

Use large simulation ensembles to build a statistical or machine-learning classifier.

Inputs:

- planet/star parameters;
- volatile inventory;
- atmosphere/ocean state;
- surface environment type.

Outputs:

- origin-of-life favorability score;
- dominant chemical pathway class;
- uncertainty maps.

### Comparative Solvent Worlds

Extend the environment framework beyond water.

Candidate solvents:

- methane/ethane Titan-like liquids;
- ammonia-water mixtures;
- supercritical CO2;
- high-salinity brines.

Core question:

- how does solvent identity reshape prebiotic chemistry and molecular stability?

### Abiotic Biosignature False Positives

Use disequilibrium chemistry to study biosignature ambiguity.

Core question:

- which abiotic pathways mimic biosignatures under prebiotic or weakly inhabited conditions?

Examples:

- CH4 + CO2 coexistence;
- O2/O3 buildup;
- N2O production;
- organic haze formation.

Connection:

- this bridges origin-of-life modeling with exoplanet biosignature interpretation.
