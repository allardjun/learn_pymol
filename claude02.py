import pymol
from pymol import cmd

def load_and_color_split(pdb_file):
   # Clear everything
    cmd.delete('all')
    cmd.reinitialize()


    # Load the PDB file and show all particles
    cmd.load(pdb_file, "protein")
    cmd.show("cartoon", "all")
    #cmd.show("spheres", "all")

    # Change the background to white
    cmd.bg_color("white")

    # Get the total number of atoms
    total_atoms = cmd.count_atoms("all")
    
    # Calculate the midpoint
    midpoint = total_atoms // 2

    print(f"Total atoms: {total_atoms}")

    # Color the first half red and the second half blue
    cmd.select("first_half", f"index 1-{midpoint}")
    cmd.select("second_half", f"index {midpoint+1}-{total_atoms}")


    cmd.select("centerline", f"index 1-389")
    cmd.select("tp_strand", f"index 390-778")
    cmd.select("fp_strand", f"index 779-1171")


    cmd.color("red", "centerline")
    cmd.color("blue", "tp_strand")
    cmd.color("green", "fp_strand")


    cmd.color("yellow", "/protein//A/SSN`0/A4")


    # Clear selections to clean up
    cmd.deselect()

    # Center the view
    #cmd.zoom()

    view = (
        -0.199920803,   -0.976015866,    0.086165324,
        -0.679188073,    0.201427713,    0.705781817,
        -0.706212342,    0.082578234,   -0.703167737,
        0.000000000,    0.000000000,  -69.845001221,
        -22.193262100,    9.702558517,   17.049539566,
        55.066329956,   84.623672485,  -20.000000000 
    )
    cmd.set_view(view)

    print(f"Loaded {pdb_file} with {total_atoms} atoms.")
    print(f"First {midpoint} atoms colored red, remaining {total_atoms - midpoint} atoms colored blue.")

# Replace 'your_structure.pdb' with the path to your actual PDB file
load_and_color_split('data/two_nucs.pdb')