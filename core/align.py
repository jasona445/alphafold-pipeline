from Bio.PDB import PDBParser
from utils.utils import *
from collections import defaultdict
import numpy as np

residue_dict = {'PHE': 'F', 'LEU': 'L', 'SER': 'S', 'TYR': 'Y', 'TER': 'X', 'CYS': 'C', 'TRP': 'W',
                'PRO': 'P', 'HIS': 'H', 'GLN': 'Q', 'ARG': 'R', 'ILE': 'I', 'MET': 'M', 'THR': 'T', 'ASN': 'N',
                'LYS': 'K', 'VAL': 'V', 'ALA': 'A', 'ASP': 'D', 'GLU': 'E', 'GLY': 'G'}


def align(pairings):
    """
    :arg pairings: list of protein pairs (WT, MUT) in PDB format.

    saves fasta of seqs for alphafold, returns masks
    """

    masks = {}
    pdbs = rel_to_abs('data/experimental-structures', os.path.dirname(__file__))
    fasta = rel_to_abs('data/proteins.fasta', os.path.dirname(__file__))

    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    with open(fasta, 'w') as f:
        for p1, p2 in pairings:
            s1, s2 = parser.get_structure(p1, os.path.join(pdbs, '{}.pdb'.format(p1))), \
                     parser.get_structure(p2, os.path.join(pdbs, '{}.pdb'.format(p2)))
            [t1, non_res1], [t2, non_res2] = trim_non_residue(s1), trim_non_residue(s2)     # lists of residues
            a1, a2 = get_alignments(p1, p2, t1, t2)                 # alignments of residues
            processed1, processed2, middle1, middle2, front_back1, front_back2 = alignment_to_seqs(a1.seq, a2.seq)
            masks[p1] = [a1, non_res1, front_back1, middle1]
            masks[p2] = [a2, non_res2, front_back2, middle2]

            # Write to the fasta file
            f.writelines(['>{}'.format(p1), '\n',
                           processed1, '\n'])
            f.writelines(['>{}'.format(p2), '\n',
                           processed2, '\n'])

    return masks


def trim_non_residue(structure):
    """
    :arg structure: PDBParser structure
    returns list of nucleotide Residue objects in structure, trims off other stuff (e.g. HOH), returns the mask n list
    """
    l = []
    m = []
    for res in structure.get_residues():
        if res.get_resname() in residue_dict.keys():
            l.append(res)
            m.append(True)
        else:
            m.append(False)

    return l, m



def get_alignments(p1, p2, t1, t2):
    """
    runs clustal alignment on the two lists of residues
    :arg p1 p2: protein 1 and protein 2 names
    :arg t1: trimmed list of residues for p1
    :arg t2: trimmed list of residues for p2
    Returns clustal alignments for p1, p2 from lists t1, t2
    """

    tmp = rel_to_abs('tmp/pair.fasta', os.path.dirname(__file__))
    tmp_align = rel_to_abs('tmp/align.al', os.path.dirname(__file__))
    clustal = rel_to_abs('scripts/clustalo', os.path.dirname(__file__))

    res1 = ''.join(map(lambda k: residue_dict[k.get_resname()], t1))
    res2 = ''.join(map(lambda k: residue_dict[k.get_resname()], t2))

    with open(tmp, 'w') as g:
        # Write to the temp file
        g.writelines(['>{}'.format(p1), '\n',
                      res1, '\n'])
        g.writelines(['>{}'.format(p2), '\n',
                      res2, '\n'])

    system_call("{} -i {} -o {} --outfmt clu --force".format(clustal, tmp, tmp_align))
    from Bio import AlignIO
    alignment = AlignIO.read(tmp_align, "clustal")

    os.remove(tmp)
    os.remove(tmp_align)

    return list(alignment)


def alignment_to_seqs(s1, s2):
    # Returns a list of seqs to go to FASTA for alphafold and the masks used to get there
    # Trims off front/back and fills in middle


    ret1, ret2 = "", "" # Processed seqs

    # Find the front and the back of each
    front1, front2 = [], []
    back1, back2 = [], []
    i = 0
    seen_res1, seen_res2 = False, False
    for r1, r2 in zip(s1, s2):
        assert not (r1 == '-' and r2 == '-')

        if r1 == '-':
            if seen_res1:
                back1.append(i)
            else:
                front1.append(i)
                seen_res1 = False
                seen_res2 = True
        elif r2 == '-':
            if seen_res2:
                back2.append(i)
            else:
                front2.append(i)
                seen_res2 = False
                seen_res1 = True
        else:
            back1, back2 = [], []
            seen_res1 = True
            seen_res2 = True

        i += 1

    # Turn them into masks
    front_back1 = [False if i in front1 + back1 else True for i in range(len(s1))]
    front_back2 = [False if i in front2 + back2 else True for i in range(len(s1))]

    middle1 = []
    middle2 = []
    # Fill in the middle and create the mask (0 for filled, 1 for ok)
    for r1, r2, m in zip(s1, s2, np.logical_xor(front_back1, front_back2)):
        if m:
            middle1.append(True)
            middle2.append(True)
        else:
            if r1 == '-':
                ret1 += r2
                ret2 += r2
                middle1.append(False)
                middle2.append(True)
            elif r2 == '-':
                ret1 += r1
                ret2 += r1
                middle1.append(True)
                middle2.append(False)
            else:
                ret1 += r1
                ret2 += r2
                middle1.append(True)
                middle2.append(True)

    return [ret1, ret2, middle1, middle2, front_back1, front_back2]



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
    align(pairings)