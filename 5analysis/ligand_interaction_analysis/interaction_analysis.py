from matplotlib import pyplot as plt
import mdtraj as md 
import numpy as np
from scipy import stats
import seaborn as sb
from collections import Counter
import pandas as pd
#############################
# User‐defined parameters
#############################

# List of ligand names. (Ensure these match the residue names in your topology files.)
ligands = ['ATP', 'ADP', 'AMP', 'CMP','RNUA']

# Distance cutoff (in nm); 0.4 nm is 4 Å.
distance_cutoff = 0.30  

# Interaction frequency cutoff for visualization (in percent)
interaction_cutoff = 40.0  

# Set stride value

stride=1

#############################
# Load trajectories for each ligand (trim first, then stride)
#############################

n_replicates = 10          # number of 500-frame segments
frames_per_rep = 500       # total frames per segment
trim_frames_per_rep = 100  # remove first 100 frames from each segment

ligand_trajectories = {}

for ligand in ligands:
    traj_path = f'/path/to/{ligand}/joined.dcd'
    top_path  = f'/path/to/{ligand}/step3_input.pdb'
    
    try:
        traj = md.load_dcd(traj_path, top=top_path)  # load *without* stride
        print(f"Loaded {traj.n_frames} frames for ligand {ligand}")
    except Exception as e:
        print(f"Error loading {traj_path} for {ligand}: {e}")
        continue

    if traj.n_frames == 0:
        print(f"No frames for {ligand}; skipping.")
        continue

    # ---- Trim first 100 frames of every 500-frame block ----
    n_frames = traj.n_frames
    frames_per_rep = n_frames // n_replicates

    trimmed_segments = []
    for i in range(n_replicates):
        start = i * frames_per_rep + trim_frames_per_rep
        end = (i + 1) * frames_per_rep
        trimmed_segments.append(traj[start:end])

    trimmed_traj = trimmed_segments[0].join(trimmed_segments[1:])
    print(f"After trimming: {trimmed_traj.n_frames} frames retained for {ligand}")

    # ---- Apply stride following trimming ----
    if stride > 1:
        trimmed_traj = trimmed_traj[::stride]
        print(f"After striding ({stride}): {trimmed_traj.n_frames} frames kept")

    ligand_trajectories[ligand] = trimmed_traj

#############################
# Compute interaction frequencies
#############################

# We will store the protein interaction counts (number of frames a residue interacts) per ligand.
# For each ligand we build a Counter keyed by the protein residue index.
interaction_counts = {ligand: Counter() for ligand in ligand_trajectories.keys()}

# Process each ligand’s trajectory
for ligand, traj in ligand_trajectories.items():
    print(f"Processing ligand {ligand}...")
    
    # Select ligand atoms (assuming the ligand’s residue name is exactly ligand)
    ligand_indices = traj.topology.select(f'resname {ligand}')
    if ligand_indices.size == 0:
        print(f"No atoms found for ligand {ligand} using selection 'resname {ligand}'.")
        continue

    # Select all protein atoms (this query remains the same across ligands)
    protein_indices = traj.topology.select('protein')
    
    # Loop over frames in the trajectory and count protein residues within the cutoff of any ligand atom.
    for frame in traj:
        # md.compute_neighbors returns, for a trajectory (or frame), an array of indices from the haystack
        # (protein atoms) that lie within the cutoff distance of any query atom (ligand atoms).
        # Since we are iterating frame‐by‐frame (each frame is a trajectory slice with one frame),
        # we take the first (and only) element of the returned list.
        neighbors = md.compute_neighbors(frame, distance_cutoff, ligand_indices, haystack_indices=protein_indices)[0]
        # Get the unique protein residue indices corresponding to the neighbor atoms.
        interacting_residues = {frame.topology.atom(idx).residue.index for idx in neighbors}
        # Update counts
        interaction_counts[ligand].update(interacting_residues)
        
    print(f"Done processing {ligand}: interactions counted over {traj.n_frames} frames.")

#############################
# Prepare data for heatmap
#############################

# Use one example topology (they should all be identical) to map residue index to a label.
example_topology = list(ligand_trajectories.values())[0].topology
# For each residue, we can create a label such as 'LYS45'
#residue_labels = {res.index: f'{res.name}{res.index}' for res in example_topology.residues}
residue_labels={res.index: str(res) for res in example_topology.residues}
# Build a dictionary where for each ligand we have a mapping of protein residue label -> % interaction frequency.
# (Interaction frequency is calculated as: (frames with interaction / total frames) * 100)
data = {}
for ligand, counter in interaction_counts.items():
    freq_dict = {}
    total_frames = ligand_trajectories[ligand].n_frames
    for res_index, count in counter.items():
        label = residue_labels[res_index]
        freq_dict[label] = (count / total_frames) * 100
    data[ligand] = freq_dict

# Create a DataFrame.
# In this DataFrame rows correspond to protein residues and columns to ligands.
df = pd.DataFrame.from_dict(data, orient='columns').fillna(0)
df.to_pickle('int_data.pkl')
print("Interaction frequency DataFrame (in %):")
print(df)

# Optionally, filter the DataFrame so that only cells with interaction frequency above the cutoff are shown.
filtered_df = df.where(df >= interaction_cutoff)

#############################
# Plot heatmap
#############################

plt.figure(figsize=(10, 20))
# Since rows are protein residues (y-axis) and columns are ligands (x-axis), we can directly use imshow.
# Note: imshow expects a 2D array; filtered_df will have NaNs where the frequency is below cutoff.
im = plt.imshow(filtered_df, cmap='YlOrRd', aspect='equal', interpolation='nearest')
plt.colorbar(im, label='Interaction Frequency (%)')

# Set tick labels
plt.xticks(ticks=np.arange(len(filtered_df.columns)), labels=filtered_df.columns, rotation=90)
plt.yticks(ticks=np.arange(len(filtered_df.index)), labels=filtered_df.index)

plt.title('Protein Residue Interaction Frequency with Different Ligands')
plt.xlabel('Ligands')
plt.ylabel('Protein Residues')
plt.tight_layout()

# Save the figure at high resolution
plt.savefig('ligand_protein_interaction_new.png', dpi=300)
#plt.show()
