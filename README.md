# FCCee_beamgas_xSections

Repository for computing FCC-ee beam-gas interaction cross sections using BDSIM.

In the `scripts` directory, scripts for running BDSIM simulations for the gas species of interest for FCC-ee (H, CO, CO2) can be found.
To run the scripts, execute the following commands:

```bash
mkdir output
bdsim --file=scripts/scriptname.gmad --batch --ngenerate={n_particles} --outfile=output/outputname
```

This will generate a root file (BDSIM output) named `outputname.root` in the `output` directory.

With the `output_analysis.py` script in the `toolkit` directory, it is possible to analize the BDSIM output and produce output files in the `toolkit` containing:
  - The interaction processes that the particles have experienced
  - The cross sections for each process

To do so, execute the following command from the `toolkit` directory:
```bash
python3 output_analysis.py --root_file=../output/outputname.root --gas={gas}
```

All possible interaction processes, along with the corresponding ProcessType (PT) and ProcessSubType (PST), can be found in the `g4_processes.csv` file.

The `main` branch is intended to provide an overview of all the possible interaction processes for the case under exam and on the corresponding cross sections. The other branches are dedicated to study more in detail individual interaction processes, using proper physics biasings.

Helpful Links:
  - BDSIM documentation: http://www.pp.rhul.ac.uk/bdsim/manual/
  - pybdsim documentation: http://www.pp.rhul.ac.uk/bdsim/pybdsim/

**NOTE**: Currently, the repository is configured for beam-gas interaction studies for FCC-ee, Z operation mode, B1 (positrons). Future updates will include the study of other operation modes (W, H, ttbar) and B2 (electrons).
