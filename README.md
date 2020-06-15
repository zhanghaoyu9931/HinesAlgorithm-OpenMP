# HinesAlgorithm-OpenMP

Parallel version of Hines algorithm (a core step in neuron simulation, see: [https://doi.org/10.1016/j.procs.2017.05.145](https://doi.org/10.1016/j.procs.2017.05.145)) written in OpenMP.

Source codes are contained in directory `src/`. `serial.cpp` is the serial program, while `parallel.cpp` the parallel. Small-scale and mediate-scale datasets are provided in `data/`.

- To compile the program: Move to directory `script/` and execute `compile.sh(.bat)` for Linux(Windows) OS.

- To run the serial program: Move to directory `script/` and execute `srun_all.sh(.bat)` for Linux(Windows) OS. Results will be written into `sresult/`. If the executable `src/serial` does not yet exist, this script will compile it for you.

- To run the parallel program: Move to directory `script/` and execute `prun_all.sh(.bat)` for Linux(Windows) OS. Results will be written into `presult/`. If the executable `src/parallel` does not yet exist, this script will compile it for you.
