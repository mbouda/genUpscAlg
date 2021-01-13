# genUpscAlg

The code in this repository (run from the driveGenAlgTest.m script) performs a series of tests comparing upscaled solutions to full-network solutions for several simple, synthetic networks.

One funciton called but not included is popJCMP.m, which is contained in the dropbox folder exactKcomp/martinCode, maintained by Jan Vanderborght at Juelich.

Generality and efficiency issues are noted in the issues section of the repository. It can be shown theoretically that all generality issues can be overcome.

An important analytical (theoretical) result contained in the mbouda/synth repository (in the file "symCombPsiS.m") is a proof that the present upscaling approach can be extended to cases where soil water potential is not uniform within a soil 'layer' (or region, more broadly). In this case, so long as a (weighted\*) average soil water potential is used, the upscaling is possible (numerical code not included here). If performed, the resulting calculations with upscaled parameters will yield correct values for (weighted\*) mean xylem water potentials in a region and the resulting flows.

\*weights are the total radial conductance of each root segment, so Kr\*S, where Kr is the soil-root conductance per unit length and S is the total segment length.
