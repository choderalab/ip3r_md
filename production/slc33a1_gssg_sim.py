"""
Script name: slc33a1_gssg_sim.py

Intended Use:
This script is intended to run a 500 ns molecular dynamics simulation of SLC33A1 
with or without GSSG bound in a lipid bilayer with waters and ions under CHARMMFF. 
It can also be adapted for different types of simulations by changing the input 
parameters in the YAML configuration.

Requirements:
- OpenMM (Version 7.7.0, maybe be incompatible with earlier or later versions)
- MDTraj for trajectory output (DCD format)
- YAML for configuration parsing
- CHARMM parameter and topology files
- JSON

Ensure the appropriate libraries are installed and that the paths to input files 
are correctly set in `sim_params.yaml` before running the script.

"""


from sys import stdout
import mdtraj as md 
import openmm as mm
import openmm.app as app
from openmm import LangevinMiddleIntegrator
from openmm.app import CharmmPsfFile, PDBFile, PDBxFile, CharmmParameterSet
import simtk.unit as unit 
import os, time, yaml, bz2
import json


### Open YAML file containing simulation parameters

with open('sim_params.yaml') as file: # Load YAML containing location of projects
        
    params = yaml.load(file,Loader=yaml.FullLoader) 

### Load PSF Topology:
    
input_filepath=params['input_paths']['input_filepath']
    
psf = CharmmPsfFile(input_filepath+'step5_input.psf')
pdb = PDBFile(input_filepath+'step5_input.pdb')

state_file = input_filepath+params['input_paths']['state_filepath']

### Load paramter files
param_filepath = params['input_paths']['param_filepath']
param_filenames = params['input_paths']['param_filenames']
sim_params=params['sim_params']
restraint_params=params['restraint_params']
output_filenames=params['output_paths']

param_paths = [os.path.join(param_filepath, file) for file in param_filenames]
params = CharmmParameterSet(*param_paths)

### Set system parameters (load from sysinfo.dat)
sys_info_path = '../input/openmm/sysinfo.dat'
with open(sys_info_path) as sysdata:
	data = json.load(sysdata)
	x,y,z=map(float,data['dimensions'][:3]) * unit.angstroms


psf.setBox(x,y,z) # Matches box size set by CHARMM

nonbonded_method = app.PME
constraints = app.HBonds
hydrogen_mass = sim_params['hydrogen_mass'] * unit.amu
temperature = sim_params['temperature']*unit.kelvin
friction = sim_params['friction']/unit.picosecond
time_step = sim_params['time_step']*unit.picoseconds
pressure = sim_params['pressure']*unit.bar
surface_tension = sim_params['surface_tension'] 

### Set up the system 

system = psf.createSystem(params,
                              nonbondedMethod=nonbonded_method,
                              constraints=constraints,
                              removeCMMotion=False,
                              hydrogenMass=hydrogen_mass,
                              )

integrator = LangevinMiddleIntegrator(temperature,        
                                     friction,
                                     time_step)

barostat = mm.MonteCarloMembraneBarostat(pressure,
                                         surface_tension, 
                                         temperature,
                                         mm.MonteCarloMembraneBarostat.XYIsotropic, 
                                         mm.MonteCarloMembraneBarostat.ZFree
                                        )
barostat.setFrequency(50)    
                           
system.addForce(barostat)

simulation = app.Simulation(psf.topology, system, integrator)
    
with bz2.open(state_file, 'rb') as infile:
        state = mm.XmlSerializer.deserialize(infile.read().decode())

simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
simulation.context.setPositions(state.getPositions())
simulation.context.setVelocities(state.getVelocities())

### Simulation parameters

nsteps= sim_params['nsteps'] 
report_freq= sim_params['report_freq'] 
chk_freq= sim_params['chk_freq'] 
traj_freq= sim_params['traj_freq']  # every 1ns at 4fs / step

current_time = simulation.context.getState().getTime() / unit.nanoseconds
total_simulation_time = nsteps*time_step / unit.nanoseconds
simulation_time = total_simulation_time - current_time
steps_left = round(simulation_time*unit.nanoseconds/time_step)

traj_freq_time = traj_freq * time_step/ unit.nanoseconds
report_freq_time = report_freq*time_step / unit.nanoseconds

chk_freq_time = chk_freq * time_step / unit.nanoseconds

### Set output stuff

# Set file names
integrator_xml_filename = output_filenames['integrator_xml_filename']
state_xml_filename = output_filenames['state_xml_filename']
state_pdb_filename = output_filenames['state_pdb_filename']
state_pdbx_filename = output_filenames['state_pdbx_filename']
system_xml_filename = output_filenames['system_xml_filename']
checkpoint_filename = output_filenames['checkpoint_filename']
traj_output_filename = output_filenames['traj_output_filename']

# write limited state information to standard out:
simulation.reporters.append(
    app.StateDataReporter(
        stdout,
        reportInterval=report_freq,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        speed=True,
        progress=True,
        remainingTime=True,
        totalSteps=steps_left,
        separator="\t",
    )
)

# Write to checkpoint files regularly:

simulation.reporters.append(app.CheckpointReporter(
    file=checkpoint_filename,
    reportInterval=chk_freq
    )
)

# Write out the trajectory

simulation.reporters.append(md.reporters.DCDReporter(
    file=traj_output_filename,
    reportInterval=traj_freq
    )
)

# Run NPT dynamics
print("Running dynamics in the NPT ensemble...")
initial_time = time.time()
simulation.step(steps_left)
elapsed_time = (time.time() - initial_time) * unit.seconds
simulation_time = nsteps * time_step
print('    Simulation took %.3f s for %.3f ns (%8.3f ns/day)' % (elapsed_time / unit.seconds, simulation_time / unit.nanoseconds, simulation_time / elapsed_time * unit.day / unit.nanoseconds))

# Save and serialize the final state
print("Serializing state to %s" % state_xml_filename)
state = simulation.context.getState(
    getPositions=True,
    getVelocities=True,
    getEnergy=True,
    getForces=True
)

state_to_write = simulation.context.getState(
    
    getPositions=True,
    
    getVelocities=True,
    
    getEnergy = True,
    
    getForces = True
    
    )

with bz2.open(state_xml_filename, "wt") as outfile:
    xml = mm.XmlSerializer.serialize(state_to_write)
    outfile.write(xml)

# Save the final state as a PDB
print("Saving final state as %s" % state_pdb_filename)
with open(state_pdb_filename, "wt") as outfile:
    PDBxFile.writeFile(
        simulation.topology,
        simulation.context.getState(
            getPositions=True,
            enforcePeriodicBox=True).getPositions(),
            file=outfile,
            keepIds=True
    )
    
# Save the final state as a PDBx file (.cif)

with open(state_pdbx_filename,'wt') as outfile:
    
    PDBxFile.writeFile(
        simulation.topology,
        
        simulation.context.getState(
            
            getPositions=True,
            
            enforcePeriodicBox=True).getPositions(),
        
        file = outfile,
        keepIds=True

        )

# Save and serialize system
print("Serializing system to %s" % system_xml_filename)
system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
with bz2.open(system_xml_filename, "wt") as outfile:
    xml = mm.XmlSerializer.serialize(simulation.system)
    outfile.write(xml)
    
    
integrator_to_write = simulation.context.getIntegrator()

# Save and serialize integrator
print("Serializing integrator to %s" % integrator_xml_filename)
with bz2.open(integrator_xml_filename, "wt") as outfile:
    xml = mm.XmlSerializer.serialize(integrator_to_write)
    outfile.write(xml)

