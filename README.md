# rate22_cse_code

The [UDfA Rate22](http://umistdatabase.net/) circumstellar chemical kinetics model describes the gas-phase chemistry in an AGB outflow with a constant mass-loss rate and outflow velocity.
This release includes the effects of
- A clumpy outflow using the porosity formalism [(Van de Sande et al. 2018)](https://ui.adsabs.harvard.edu/abs/2018A&A...616A.106V/abstract)
- Stellar UV photons [(Van de Sande & Millar 2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...873...36V/abstract)
- Close-by stellar companion UV photons [(Van de Sande & Millar 2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.510.1204V/abstract)



## Contents
- The FORTRAN77 code, located in the folder `code/`.
- A makefile, `my_makefile`.
- A `.rates` file, containing the chemical reaction network: `rate22_revised.rates`.
- A `.specs` file, containing all species and parent abundances: `rate22_revised.specs`. 
- A perl script to compile new ODEs: `rate12cse.pl`.
- A test input file, `test_input.txt`.

A `.rates` file excluding the reactions identified by [Tinnaci et al. (2021)](https://ui.adsabs.harvard.edu/abs/2023ApJS..266...38T/abstract) as endothermic (done by putting their rates to zero) can be found on the [website](http://umistdatabase.net/).


### Compiling the model
Running `./my_makefile` compiles the code to the executable `csmodel`. Note that a fortran compiler (e.g., gfortran) is necessary to do so.


## Running a model

The model takes input from an input file.
The command `./csmodel (inputfile)' calculates your desired model.



### Input parameters 

The input parameters and their units are listed in `test_input.txt`. It's essential that you keep the basic format of the input file.

Special care needs to be taken with
1. The clumping parameters
  - `CLUMPMODE` should be either `SMOOTH` for a smooth outflow (classical CSE model) or `POROSITY` for a clumpy (porous) model.
  - `FVOL` is the clump volume filling factor, setting raction of the total volume occupied by the clumps. 
  Therefore, 0 < `FVOL` < 1. `FVOL = 0` will result in an error, `FVOL = 1` is equivalent to a smooth outflow.
  - `FIC` is the interclump density contrast between the interclump component and the mean density.
  Again, 0 < `FIC` < 1. `FIC = 0` will result in an error, `FIC = 1` is equivalent to a smooth outflow.
  - `L` is the size of the clumps at the stellar surface. Therefore, `L` should be smaller than the stellar radius. 
  
The porosity formalism implemented in the model assumes a constant `FVOL`, which results in uniformly expanding clumps. 
More information can be found in [Van de Sande et al. 2018](https://ui.adsabs.harvard.edu/abs/2018A&A...616A.106V/abstract)

2. The inner photon parameters
  - `ISTELLAR` and `IBIN` turn stellar and companion UV photons on and off.
  - When turning on `IBIN`, the type of stellar companion needs to be specified. Make sure to match the companion's radius to its temperature:
    - Red dwarf companion: `RBIN` = 1.53e10 cm, `TBIN` = 4000 K
    - Solar-like companion: `RBIN` = 8.14e10 cm, `TBIN` = 6000 K
    - White dwarf companion: `RBIN` = 6.96e8 cm, `TBIN` = 10000 K
  - `RDUST` is the dust condensation radius. Dust formation isn't included in the model, this is the radius where the dust is assumed to have fully formed. 
  The starting radius `R_INNER_CHEM` cannot be smaller or equal to `RDUST`. For best results, choose the initial radius as close to `RDUST` as possible.

### Parent species
The parent species are listed at the bottom of the `.specs` file.
Their units are fractional abundance relative to H.

## Adding inner photon reactions

The photoreactions from an internal stellar and companion UV source can be found on the [website](http://umistdatabase.net/).
If you want to include these, please follow these steps:
1. Add the desired reactions to a new `.rates` file
2. Write a new ODE file using `./rate12cse.pl (name of your new rates file) -o acodes.f`
3. Move `acodes.f` to the `code/` folder and recompile the model
4. The perl script also writes a new `.specs` file. It works in mysterious ways, more likely than not the order of the species is changed. Make sure to use this `.specs` file!



## Contact

If you have any comments or isues running the model, please contact Marie Van de Sande at "mvdsande at strw.leidenuniv.nl".

## Acknowledgements

The code is free to use. Please cite the [the Rate22 paper](http://umistdatabase.net/)
- [Insert citation later]
When including the effects of a clumpy outflow or internal photons, please also cite the relevant papers. 
- Clumping: [Van de Sande et al., 2018, A&A, 631, A106]((https://ui.adsabs.harvard.edu/abs/2018A&A...616A.106V/abstract))
- Stellar UV photons: [Van de Sande & Millar, 2019, ApJ, 873, 36](https://ui.adsabs.harvard.edu/abs/2019ApJ...873...36V/abstract)
- Companion UV photons [Van de Sande & Millar, 2022, MNRAS, 510, 1204](https://ui.adsabs.harvard.edu/abs/2022MNRAS.510.1204V/abstract)

