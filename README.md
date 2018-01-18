# QE-GIPAW
This is the official Git repository of the GIPAW code for Quantum-Espresso.

## Features
* NMR shielding tensors, EFG tensors
* EPR g-tensor, hyperfine couplings
* Mössbauer
* Norm-conserving (g-tensor only), USPP and PAW
* LDA and GGA functionals
* isolated and periodic systems
* integration with Quantum-Environment (solvent effects)


## Authors and contributors
D. Ceresoli, A. P. Seitsonen, U. Gerstmann, E. Kucukbenli, S. de Gironcoli, P. Giannozzi, N. Varini


## Build instructions:
### From Quantum-Espresso distribution
1. Configure and compile QE as usual, then:
2. ```make gipaw```


### Stand-alone 
1. Configure and compile QE as usual, then:
2. ```git clone https://github.com/dceresoli/qe-gipaw.git```
3. ```cd qe-gipaw```
4. ```./configure --with-qe-source="QE folder containing make.inc"```
5. ```make```

### Source releases
Official source releases are found [here](https://github.com/dceresoli/qe-gipaw/releases).


### News
* Jan 18, 2018: added ```tools/NMR_tensor-convention.py``` to convert NMR shielding tensors anysotropy in various conventions
* Jan 18, 2018: added ```tools/EPR_diagonalize_g-tensor``` to output the principal components and directions of the g-tensor, accordin to Weyl, Bolton, *Principles of electron paramagnetic resonance*, ch. 4.4
* Jan 18, 2018: fixed bug in symmetrization of g-tensor




