from phonopy import load

from vibes.brillouin import get_bands_and_labels
from vibes.structure.convert import to_Atoms

phonon = load("phonopy.yaml")

primitive = to_Atoms(phonon.primitive)

# this is the default collection of bandpaths, can be anything else
paths = primitive.cell.get_bravais_lattice().special_path.split(",")
bands, labels = get_bands_and_labels(primitive, paths, latex=False)

# compute the new bandstructure and save to file
phonon.run_band_structure(bands, labels=labels)
phonon.write_yaml_band_structure(filename="custom_band.yaml")

