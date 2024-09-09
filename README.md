# fconv, dsx and friends

* Back in my PhD-time at the Philipps Universität Marburg, I made some tool available as open source, namely fconv, dsx and hotspotsX.
* The tools were available under www.agklebe.de, but that's no longer the case (https://agklebe.pharmazie.uni-marburg.de/?id=11 only links publications, but not the tools).
* As I got some requests to make it available again, I decided to put it here on GitHub. 
  * All code is provided as of 2012. It's horrible, unmaintainable spaghetti-code (don't look at it, or at least don't judge me on it), but it works fast and rock-solid.
  * All code is under GPL license.
* Some other more or less related tools from the old time (e.g. DSFP) are still missing, but I might add them here once I recovered the latest sources... 

## Instructions

* Shared libs are located in the folder C_libs
* The folder fconv holds:
  * Separate README, INSTALL and LICENSE files (as of 2012)
  * bin folder with statically-linked 32bit and 64bit executables for linux (precompiled MacOs- and Windows-versions were available on the old www.agklebe.de but I currently had no suitable machine to generate these, so for now you have to rely on building from source, if you need it)
  * src folder with makefile (see INSTALL.txt)
* The folder dsx holds:
  * Separate README, INSTALL and LICENSE files (as of 2012)
  * bin folder with statically-linked 32bit and 64bit executables for linux (precompiled MacOs- and Windows-versions were available on the old www.agklebe.de but I currently had no suitable machine to generate these, so for now you have to rely on building from source, if you need it)
  * src folder with makefile (see INSTALL.txt)
  * pdb_pot_0511 folder with the PDB-derived potentials (as used for the publication)
  * pdb_pot_pharmacophoric_0511 folder with PDB-derived potentials on a atom-type set reduced to pharmacophoric characteristics
  * NOTE: CSD-derived potentials are not publicly available. At least back in 2012, these required a valid CSD license. If you want them, please contact the CCDC (maybe they are OK with making them freely available, cause they're based on such an outdated version of the CSD).
* The folder hsx holds:
  * Separate README, INSTALL and LICENSE files (as of 2012)
  * bin folder with statically-linked 32bit and 64bit executables for linux (precompiled MacOs- and Windows-versions were available on the old www.agklebe.de but I currently had no suitable machine to generate these, so for now you have to rely on building from source, if you need it)
  * src folder with makefile (see INSTALL.txt)
  * An example mapping file ACC_DON_AnD_HYD_ARO_map.def