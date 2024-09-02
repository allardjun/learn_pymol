import pymol
from pymol import cmd

def load_and_color_split(pdb_file):
   # Clear everything
    cmd.delete('all')
    cmd.reinitialize()


    # Load the PDB file and show all particles
    cmd.load(pdb_file, "chromosomal_dna")
    #cmd.show("spheres", "all")

    # Change the background to white
    cmd.bg_color("white")

    # Get the total number of atoms
    total_atoms = cmd.count_atoms("all")
    
    print(f"Total atoms: {total_atoms}")

    # Centerline
    cmd.select("centerline", f"index 1-389")
    cmd.hide("everything", "centerline")

    # DNA
    cmd.select("tp_strand", f"index 390-778")
    cmd.select("fp_strand", f"index 779-1171")

    cmd.hide("licorice", "tp_strand")
    cmd.hide("licorice", "fp_strand")

    cmd.show("sticks", "tp_strand")
    cmd.show("sticks", "fp_strand")

    DNA_COLOR = "0x21409A"

    cmd.color(DNA_COLOR, "tp_strand")
    cmd.color(DNA_COLOR, "fp_strand")

    cmd.set("stick_radius", 0.15, "all")

    # Nucleosomes
    cmd.select("nucleosomes", "/chromosomal_dna//A/SSN`0/A4")
    cmd.show("sphere", "nucleosomes")
    cmd.color("0xEA0A8C", "nucleosomes")
    cmd.set("sphere_scale", 1.7, "nucleosomes")



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

# Replace 'your_structure.pdb' with the path to your actual PDB file
load_and_color_split('data/two_nucs.pdb')

if True:
    cmd.set("ray_trace_mode", 1)  # Higher quality ray-tracing
    cmd.ray()
    cmd.png("figs/high_quality_image.png", dpi=300, ray=1)
