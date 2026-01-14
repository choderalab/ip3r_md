#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import stdout
import mdtraj as md 
import openmm as mm
import openmm.app as app
from openmm import LangevinMiddleIntegrator
from openmm.app import CharmmPsfFile, PDBFile, PDBxFile, CharmmParameterSet
import simtk.unit as unit 
import os, time, yaml, bz2, argparse, logging
from openmm import CustomExternalForce # Will be used to create custom restraint force
import json

### Set up logging apparatus

parser = argparse.ArgumentParser()
parser.add_argument('-logFile', dest='logFile', required=True)
args = parser.parse_args()

# Create and configure logger

logging.basicConfig(filename=args.logFile,
                    format='%(asctime)s %(message)s',
                    filemode='w')
logger = logging.getLogger()

### Load PSF Topology:
    
input_filepath='../equil/openmm/'
    
psf = CharmmPsfFile(input_filepath+'step3_input.psf')
pdb = PDBFile(input_filepath+'step3_input.pdb')

state_file = input_filepath+'step4_equilibration.rst.bz2'

### Load paramter files
param_filepath = '../equil/toppar/'
param_filenames = ['par_all36m_prot.prm', 'top_all36_prot.rtf',
                       'par_all36_lipid.prm', 'top_all36_lipid.rtf',
                       'toppar_water_ions.str','par_all36_cgenff.prm','top_all36_cgenff.rtf','par_all36_na.prm','top_all36_na.rtf','toppar_all36_na_nad_ppi.str']

param_paths = [os.path.join(param_filepath, file) for file in param_filenames]
params = CharmmParameterSet(*param_paths)

### Set system parameters

sys_info_path = '../equil/openmm/sysinfo.dat'

with open(sys_info_path) as sysdata:
	data = json.load(sysdata)
	x,y,z=map(float,data['dimensions'][:3])*unit.angstroms

psf.setBox(x,y,z) # Matches box size set by CHARMM

nonbonded_method = app.PME
constraints = app.HBonds
hydrogen_mass = 4.0 * unit.amu
temperature = 303.15*unit.kelvin
friction = 1/unit.picosecond
time_step = 0.004*unit.picoseconds
pressure = 1*unit.bar
surface_tension = 0 # From Alex 

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

logger.info('Getting state from state file')
    
with bz2.open(state_file, 'rb') as infile:
        state = mm.XmlSerializer.deserialize(infile.read().decode())

logger.info('Setting state from state file...')
simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
simulation.context.setPositions(state.getPositions())
simulation.context.setVelocities(state.getVelocities())

# simulation.context.setPositions(pdb.positions)
# simulation.context.setVelocitiesToTemperature()


### Simulation stuff

nsteps= 125000000 # 500 ns, 100*10^6 fs
report_freq= 25000 # every 0.1 ns
chk_freq= 12500000 # every 50 ns
traj_freq= 250000  # every 1ns at 4fs / step

#nsteps= 250000 # 1 ns, 1*10^6 fs
#report_freq= 25000 # every 0.1 ns
#chk_freq= 125000 # every 0.5 ns
#traj_freq= 25000  # every 0.1 ns at 4fs / step

current_time = simulation.context.getState().getTime() / unit.nanoseconds
total_simulation_time = nsteps*time_step / unit.nanoseconds
simulation_time = total_simulation_time - current_time
steps_left = round(simulation_time*unit.nanoseconds/time_step)

traj_freq_time = traj_freq * time_step/ unit.nanoseconds
report_freq_time = report_freq*time_step / unit.nanoseconds

chk_freq_time = chk_freq * time_step / unit.nanoseconds

logger.info(f'\t{steps_left:.0f} steps, {simulation_time:.0f} ns, will be run for a total simulation time of {total_simulation_time:.3f} ns \n'
            f'\tSaving a frame to traj_file every {traj_freq} step(s), or every {traj_freq_time:.3f} ns \n'
            f'\tWriting checkpoint file every {chk_freq} step(s), or every {chk_freq_time:.3f} ns \n'
            f'\tWriting state info to logfile every {report_freq} step(s), or every {report_freq_time:.3f} ns \n'
            )

### Apply custom restraint forces (for equilibration)

restraint_type = 'backbone'
restraint_k=10 # 20 kJ/mol/nm^2

posresPROT = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2;')
posresPROT.addGlobalParameter('k', restraint_k)
posresPROT.addPerParticleParameter('x0')
posresPROT.addPerParticleParameter('y0')
posresPROT.addPerParticleParameter('z0')

crd = simulation.context.getState(getPositions=True).getPositions()
system = simulation.context.getSystem()

atoms = [atom for atom in simulation.topology.atoms()]

for i, atom_crd in enumerate(crd):
    
    if (atoms[i].residue.id == ('1909') or 
        atoms[i].residue.id == ('1910') or 
        atoms[i].residue.id == ('1911') or
        atoms[i].residue.id == ('1912') or
        atoms[i].residue.id == ('1913') or
        atoms[i].residue.id == ('1914') or
        atoms[i].residue.id == ('1915') or 
        atoms[i].residue.id == ('2214') or 
        atoms[i].residue.id == ('2215') or 
        atoms[i].residue.id == ('2216') or 
        atoms[i].residue.id == ('2549') or
        atoms[i].residue.id == ('2550') or
        atoms[i].residue.id == ('2551') or 
        atoms[i].residue.id == ('2634') or
        atoms[i].residue.id == ('2635') or
        atoms[i].residue.id == ('2636')
        ) and (atoms[i].residue.chain == atoms[0].residue.chain or # Checks if the residue is in chain A.
               
               atoms[i].residue.chain == atoms[5000].residue.chain # Checks if the residue is in chain B.
               
               ) and (atoms[i].name in ('CA','C','N')): # Applies the restraints to backbone atoms only.
        
        print(i,crd)
        posresPROT.addParticle(i,atom_crd.value_in_unit(unit.nanometers))
        
force_id = system.addForce(posresPROT)

simulation.context.reinitialize(preserveState=True)

forces = [force for force in simulation.context.getSystem().getForces()]

print(f'Successfully added force: {forces[force_id]}')

### Set output stuff

# Set file names
integrator_xml_filename = "integrator.xml.bz2"
state_xml_filename = "state.xml.bz2"
state_pdb_filename = "output.pdb"
state_pdbx_filename = "output.cif"
system_xml_filename = "system.xml.bz2"
checkpoint_filename = "output.chk"
traj_output_filename = "output.dcd"

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


