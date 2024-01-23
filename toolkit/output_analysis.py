"""
Produces a .json output containing information on Coulomb scattering cross section and angular kicks
--------------------------------------------------------------------------------

Using the root file (BDSIM output), this script produces an output file containing information:
- Coulomb scattering cross section
- Coulomb scattering angular kicks, delta_px and delta_py [rad]

*--Required--*
- **root_file** *(str)*: Path to the root file (BDSIM output).
- **gas** *(str)*: Name of the gas (H, CO, CO2).

*--Optional--*


"""
import ROOT
import pybdsim
import numpy as np
import scipy.constants
import pandas as pd
import json
from generic_parser import EntryPointParameters, entrypoint

# Constants --------------------------------------------------------------------
P0C = 45.6e9 # [eV], FCC-ee Z operation mode

N_A = scipy.constants.N_A # Avogadro's number

MM_H = 1.00784 # Molar mass H [g/mol]
MM_CO = 28.01 # Molar mass CO [g/mol]
MM_CO2 = 44.01 # Molar mass CO2 [g/mol]

RHO_SOLID_H = RHO_SOLID_CO = RHO_SOLID_CO2 = 8.96 # [g/cm^3]

L_TARGET = 0.1 # [cm]


# Script arguments -------------------------------------------------------------
def get_params():
    params = EntryPointParameters()
    params.add_parameter(
        name="root_file",
        type=str,
        required=True,
        help="Path to the root file.",
    )
    params.add_parameter(
        name="gas",
        type=str,
        required=True,
        help="Name of the gas (H, CO, CO2).",
    )

    return params


# Entrypoint -------------------------------------------------------------------
@entrypoint(get_params(), strict=True)
def main(inp):
    df_g4_processes = pd.read_csv('../g4_processes.csv')
    output_data = pybdsim.Data.Load(inp.root_file)
    samplerData = pybdsim.Data.SamplerData(output_data, 1).data

    n_primaries = get_n_primaries(output_data)
    primaries_dict = {'n_primaries': n_primaries, 'Part_ID': -11, 'p0c': P0C} 

    weights = samplerData['weight']

    trajectories = get_trajectories(output_data, n_primaries)

    mask_CoulombScat = get_CoulombScat_processes(trajectories)

    # For elastic Coulomb scattering the conversion BDSIM --> Xsuite is not necessary
    # as delta=0 for elastic processes and px = xp*(1+delta) => px = xp
    delta_px, delta_py = get_angular_kicks(output_data, mask_CoulombScat, 'xsuite')

    number_surface_density = get_number_surface_density(inp.gas) # [cm^-2]

    n_CoulombScat = sum(weights[mask_CoulombScat])
    cross_section_CoulombScat = n_CoulombScat / (n_primaries * number_surface_density)
    cross_section = {'CoulombScat': cross_section_CoulombScat}

    angular_kicks_dict = {'n_CoulombScat': len(weights[mask_CoulombScat]), 'delta_px': delta_px.tolist(), 'delta_py': delta_py.tolist()}
    output_analysis = {'Primaries': primaries_dict, 'Cross section': cross_section, 'Angular kicks': angular_kicks_dict}

    return save_to_json(inp.gas, output_analysis)


def get_n_primaries(output_data):

    return len(pybdsim.Data.SamplerData(output_data, 0).data['n'])


def get_trajectories(output_data, n_primaries):

    trajectories = []
    keys_of_interest = ['px', 'py', 'postPT', 'postPST']

    for i in range(n_primaries):
        trajectory = pybdsim.Data.TrajectoryData(output_data, i).trajectories[0]
        trajectory_reduced = subset_dict = {key: trajectory[key] for key in keys_of_interest if key in trajectory}
        trajectories.append(trajectory_reduced)

    return np.array(trajectories)


def get_CoulombScat_processes(trajectories):

    mask_CoulombScat = [False] * len(trajectories)

    for index, trajectory in enumerate(trajectories):

        # I am assuming that taking the post process and subprocess type is correct
        # to double check
        processType = trajectory['postPT']
        processSubType = trajectory['postPST']

        if np.any((processType == 2) & (processSubType == 1)): # CoulombScat has PT=2 and PST=1
            mask_CoulombScat[index] = True

    return np.array(mask_CoulombScat)


def get_angular_kicks(output_data, mask_CoulombScat, format):

    samplerData = pybdsim.Data.SamplerData(output_data, 1).data

    delta_xp = samplerData['xp'][mask_CoulombScat] - 0
    delta_yp = samplerData['yp'][mask_CoulombScat] - 0

    if format == 'bdsim':
        return delta_xp, delta_yp

    elif format == 'xsuite':
        delta = (samplerData['p'][mask_CoulombScat]*1e9) / P0C - 1
        delta_px = delta_xp * (1+delta)
        delta_py = delta_yp * (1+delta)

        return delta_px, delta_py
    
    else:
        raise ValueError(f'Format {format} not recognized.')


def get_number_surface_density(gas):

    if gas == 'H':
        return L_TARGET * N_A * RHO_SOLID_H / MM_H # [cm^-2]
    elif gas == 'CO':
        return L_TARGET * N_A * RHO_SOLID_CO / MM_CO # [cm^-2]
    elif gas == 'CO2':
        return L_TARGET * N_A * RHO_SOLID_CO2 / MM_CO2 # [cm^-2]
    else:
        raise ValueError(f'Gas {gas} not recognized.')
    

def save_to_json(gas, output_analysis):

    with open(f'FCCee_z_CoulombScat_{gas}.json', 'w') as jsonfile:
        json.dump(output_analysis, jsonfile, indent=4)

    return
        

# Script Mode ------------------------------------------------------------------
if __name__ == '__main__':
    main()