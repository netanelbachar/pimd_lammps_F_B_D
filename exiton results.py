
##### Boson3 ##########
print ("This are the results for the bosonic energies of Moire One Well in meV, Density 3(so divide by 3), at different temperatures")
Benergy10 = [-273.80588814799756, -256.7142983400086, -251.46671274677112, -249.45408859182115, -247.10178069683235, -246.54960287398632]
Benergy6 = [-260.07426405614274, -250.55983845109313, -248.12931074310595, -247.31526776332356, -246.58926793833965, -246.0627986151288]
Benergy4 =[-252.21363518100978, -247.0369132205172, -245.94483284465275, -245.36965644986068, -245.1490359658642, -245.07692567130405]



##### Boson10 ##########





############ temperatures for excitons ########### below
kb = 3.16683 * 10 ** (-6)    # Boltzman in Hartree/K
_kb = 1/kb
temperatures_in_lammps = [310.186, 103.396, 77.547, 51.698, 31.019, 10.340, 5.170]
beta_for_graphs = [_kb/temperatures_in_lammps[0], _kb/temperatures_in_lammps[1], _kb/temperatures_in_lammps[2],
                   _kb/temperatures_in_lammps[3], _kb/temperatures_in_lammps[4], _kb/temperatures_in_lammps[5],
                   _kb/temperatures_in_lammps[6]]
beta = [1018.012374625699, 3054.017432353737, 4072.0232431383165, 6108.034864707474, 10179.992470281022, 30538.992885459094, 61077.98577091819]
beta_normalized = [1, 3, 4, 6, 10, 30, 60]
print (beta_for_graphs)
############ temperatures for excitons ########### above