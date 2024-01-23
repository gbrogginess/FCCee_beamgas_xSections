"""
Produces .json output files containing information on interaction processes and cross sections
--------------------------------------------------------------------------------

Using the root file (BDSIM output), this script produces output files containing information on:
- interaction processes
- cross sections

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

    n_primaries = get_n_primaries(output_data)
    primaries_dict = {'n_primaries': n_primaries, 'Part_ID': -11, 'p0c': P0C} 

    trajectories = get_trajectories(output_data, n_primaries)

    interaction_processes = get_interaction_processes(trajectories, df_g4_processes)
    interaction_processes = {'Primaries': primaries_dict, 'Processes': interaction_processes}

    number_surface_density = get_number_surface_density(inp.gas)

    cross_sections = {key: value / (n_primaries * number_surface_density) for key, value in interaction_processes['Processes'].items()}
    total_cross_section = {'Total': sum(cross_sections.values())}
    cross_sections = {**total_cross_section, **cross_sections}
    cross_sections = {'Primaries': primaries_dict, 'Cross sections': cross_sections}

    return save_to_json(inp.gas, interaction_processes, cross_sections)


def get_n_primaries(output_data):

    return len(pybdsim.Data.SamplerData(output_data, 0).data['n'])


def get_trajectories(output_data, n_primaries):

    trajectories = []

    for i in range(n_primaries):
        trajectory = pybdsim.Data.TrajectoryData(output_data, i).trajectories[0]
        trajectories.append(trajectory)

    return trajectories


def get_interaction_processes(trajectories, df_g4_processes):

    process_names = df_g4_processes['TypeName'].values
    interaction_processes = {process_name: 0 for process_name in process_names}

    for trajectory in trajectories:

        # I am assuming that taking the post process and subprocess type is correct
        # to double check
        processType = trajectory['postPT']
        processSubType = trajectory['postPST']

        mask_transport = np.abs(processType) == 1
        processType = processType[~mask_transport]
        processSubType = processSubType[~mask_transport]

        for i, j in zip(processType, processSubType):

            mask_process = (i == df_g4_processes['ProcessType']) & (j == df_g4_processes['SubType'])
            process_name = df_g4_processes[mask_process]['TypeName'].iloc[0]

            if process_name in interaction_processes:
                interaction_processes[process_name] += 1
            else:
                raise ValueError(f'Process {process_name} not found in list of processes.')
            
    interaction_processes = dict(sorted(interaction_processes.items(), key=lambda x: x[1], reverse=True))

    return interaction_processes


def get_number_surface_density(gas):

    if gas == 'H':
        return L_TARGET * N_A * RHO_SOLID_H / MM_H # [cm^-2]
    elif gas == 'CO':
        return L_TARGET * N_A * RHO_SOLID_CO / MM_CO # [cm^-2]
    elif gas == 'CO2':
        return L_TARGET * N_A * RHO_SOLID_CO2 / MM_CO2 # [cm^-2]
    else:
        raise ValueError(f'Gas {gas} not recognized.')
    

def save_to_json(gas, interaction_processes, cross_sections):

    with open(f'{gas}_interaction_processes.json', 'w') as jsonfile:
        json.dump(interaction_processes, jsonfile, indent=4)

    with open(f'{gas}_cross_sections.json', 'w') as jsonfile:
        json.dump(cross_sections, jsonfile, indent=4)

    return
        

# Script Mode ------------------------------------------------------------------
if __name__ == '__main__':
    main()