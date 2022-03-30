from pprint import pprint

from utils.utils import *
from Bio.PDB import PDBParser
import os
from collections import defaultdict
import numpy as np

def align_residues(proteins):
    """
    Takes in a list of *protein pairs* [[wt, mut], [wt, mut], ...]
    Saves an aligned fasta in /data

    returns alignments and masks
    ret = {protein: [alignment, mask] ...}
    """
    residue_dict = {'PHE': 'F', 'LEU': 'L', 'SER': 'S', 'TYR': 'Y', 'TER': 'X', 'CYS': 'C', 'TRP': 'W',
                    'PRO': 'P', 'HIS': 'H', 'GLN': 'Q', 'ARG': 'R', 'ILE': 'I', 'MET': 'M', 'THR': 'T', 'ASN': 'N',
                    'LYS': 'K', 'VAL': 'V', 'ALA': 'A', 'ASP': 'D', 'GLU': 'E', 'GLY': 'G'}

    residue_dict = defaultdict(lambda: '', residue_dict)  # Default val gets rid of ligands

    pdbs = rel_to_abs('data/experimental-structures', os.path.dirname(__file__))
    fasta = rel_to_abs('data/proteins.fasta', os.path.dirname(__file__))
    tmp = rel_to_abs('tmp/pair.fasta', os.path.dirname(__file__))
    tmp_align = rel_to_abs('tmp/align.al', os.path.dirname(__file__))
    clustal = rel_to_abs('scripts/clustalo', os.path.dirname(__file__))

    ret = {}

    with open(fasta, 'w') as f:
        for pair in proteins:
            parser = PDBParser(PERMISSIVE=True, QUIET=True)
            protein1, protein2 = pair
            res1 = ''.join(map(lambda k: residue_dict[k.get_resname()], parser.get_structure(protein1, os.path.join(pdbs, '{}.pdb'.format(protein1))).get_residues()))
            res2 = ''.join(map(lambda k: residue_dict[k.get_resname()], parser.get_structure(protein2, os.path.join(pdbs, '{}.pdb'.format(protein2))).get_residues()))

            with open(tmp, 'w') as g:
                # Write to the temp file
                g.writelines(['>{}'.format(protein1), '\n',
                              res1, '\n'])
                g.writelines(['>{}'.format(protein2), '\n',
                              res2, '\n'])

            system_call("{} -i {} -o {} --outfmt clu --force".format(clustal, tmp, tmp_align))
            from Bio import AlignIO
            align = AlignIO.read(tmp_align, "clustal")

            a1, a2 = list(align)

            ret[protein1] = [a1, generate_mask(a1, a2)]
            ret[protein2] = [a2, generate_mask(a1, a2)]

            mask = []
            for i, char in enumerate(a1.seq) if len(res1) < len(res2) else enumerate(a2.seq):
                if char != "-":
                    mask.append(i)

            res1_aligned = "".join(map(a1.seq.__getitem__, mask))
            res2_aligned = "".join(map(a2.seq.__getitem__, mask))

            # Fill missing intermediate residues in experimental data with corresponding residue from other paired seq

            for i, (r1, r2) in enumerate(zip(res1_aligned, res2_aligned)):
                if r1 == '-':
                    assert r2 != '-'
                    res1_aligned = res1_aligned[:i] + r2 + res1_aligned[i+1:]

                if r2 == '-':
                    assert r1 != '-'
                    res2_aligned = res2_aligned[:i] + r1 + res2_aligned[i+1:]

            # Write to the fasta file
            f.writelines(['>{}'.format(protein1), '\n',
                          res1_aligned, '\n'])
            f.writelines(['>{}'.format(protein2), '\n',
                          res2_aligned, '\n'])

            assert len(res1_aligned) == len(res2_aligned)

        os.remove(tmp)
        os.remove(tmp_align)

    return ret


def generate_mask(a1, a2):
    return np.array([0 if r1 == '-' or r2 == '-' else 1 for r1, r2 in zip(a1.seq, a2.seq)])


if __name__ == '__main__':
    pairings = \
        [['3otf', '4hbn'],
         ['1bz4', '1gs9'],
         ['1bz4', '1gs9'],
         ['1or3', '1gs9'],
         ['3s4m', '3s5e'],
         ['2z3y', '5l3d'],
         ['3k8y', '3oiw'],
         ['1p5f', '2rk4'],
         ['4wxq', '6ou0'],
         ['3vf0', '5l0h'],
         ['1n5u', '1hk3'],
         ['4ye6', '4ye9'],
         ['4aoh', '4ahm'],
         ['3sd6', '4gjf'],
         ['1r3q', '1jph'],
         ['6s2m', '5n4p']]

    pprint(align_residues(pairings))