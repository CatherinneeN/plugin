from chimerax.core.commands import CmdDesc
from chimerax.atomic import AtomsArg

# imports for sequence alignment
from chimerax.alignment_algs import SmithWaterman
from chimerax import sim_matrices

__version__ = "0.1.0"

hello_world_desc = CmdDesc()


class DistanceCalculator:
    def __init__(self, session):
        self.session = session
        self.distances = []

    # Add matrix
    def calc_distances(self, a1, a2):
        matrix = sim_matrices.matrix("BLOSUM-62", self.session.logger)
        protein_seq = "QCXTSI"
        # TSICS
        # code for extracting the proteins and their sequences
        for model in self.session.models:
            model_id = model.id_string

            # # Add this if we want to do seq align on specific model
            # # this will be added to UI
            # if model_id != 1:
            #     continue

            self.session.logger.info(model_id)

            # Sequence alignment

            if hasattr(model, 'chains') and model.chains:
                for chain in model.chains:
                    # chain_id = chain.id_string figure out to print chain_ID as in the /A
                    chain_sequence = chain.characters
                    # self.session.logger.info(chain_id)
                    self.session.logger.info(chain_sequence)
                    # Parameters provided on chimeraX example
                    score = SmithWaterman.align(protein_seq, chain_sequence, matrix, 1, 0.333)
                    self.session.logger.info(str(score))

            else:
                self.session.logger.info("Not found")

        # search with sequence alignment for CCTSIC we prefer local search engines for SA
        # Adding pseudo-bonds

        # Make a new group
        pbg = self.session.pb_manager.get_group("crosslinks")

        if pbg.id is None:  # If the group is new, add it to the session
            self.session.models.add([pbg])

        pbg.display = True

        # Assign atom pairs to each atom
        a1 = AtomsArg.parse(a1, self.session)
        a2 = AtomsArg.parse(a2, self.session)
        pb = pbg.new_pseudobond(a1[0][0], a2[0][0])
        pb.radius = 0.1
        self.session.logger.info(str(pb))   # result: /A ASN 18 CA ↔ GLN 5 CA

        # Measure distance
        distance = self.session.pb_dist_monitor.add_group(pbg)
        # Print distance
        self.session.logger.info(str(distance))  # result: /A ASN 18 CA ↔ GLN 5 CA


def hello_world(session):
    calculator = DistanceCalculator(session)

    # # Define atom pairs (MULTIPLE)
    # atom_pairs = [
    #     ("#1/A:52@ca", "#1/A:44@ca"), ("#1/A:125@ca", "#1/D:138@ca"),
    #     ("#1/G:148@ca", "#1/G:157@ca"), ("#1/G:183@ca", "#1/G:176@ca"),
    #     ("#1/B:58@ca", "#1/E:127@ca"), ("#1/K:133@ca", "#1/H:127@ca"),
    #     ("#1/B:337@ca", "#1/K:334@ca"), ("#1/C:35@ca", "#1/F:62@ca"),
    #     ("#1/I:159@ca", "#1/L:212@ca"), ("#1/C:88@ca", "#1/L:196@ca")
    # ]

    # Define atom pairs
    atom_pairs = [
        ("/A:18@ca", "/A:5@ca")
    ]

    # Calculate distances for each atom pair
    for atom1, atom2 in atom_pairs:
        calculator.calc_distances(atom1, atom2)

