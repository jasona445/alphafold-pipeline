import pandas as pd
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB import Superimposer
from utils.utils import *
import itertools
from core.align import *
from utils.lddt import lddt


def get_residues(proteins, folded=False, extension='pdb'):
    """
    Get residues from proteins

    :arg proteins: list of protein names in PDB format
    :arg folded: if True, looks for structures in data/folded-structures instead of data/experimental-structures
    :arg extension: pdb, ent, etc...

    Returns dictionary with protein name keys and lists of atoms values
    """

    structure_folder = rel_to_abs('data/experimental-structures', os.path.dirname(__file__)) if not folded else \
        rel_to_abs('data/folded-structures', os.path.dirname(__file__))
    parser = PDBParser(PERMISSIVE=True, QUIET=True)

    ret = {}
    for protein in proteins:
        if folded:
            path = os.path.join(structure_folder, '{}/ranked_0.{}'.format(protein, extension))
        else:
            path = os.path.join(structure_folder, '{}.{}'.format(protein, extension))
        structure = parser.get_structure(protein, path)
        ret[protein] = list(structure.get_residues())

    return ret


def trim_atoms(atoms1, atoms2):
    masked_atoms1 = []
    masked_atoms2 = []
    for a_s, fa_s in zip(atoms1, atoms2):
        aids = [x.get_id() for x in a_s]
        faids = [x.get_id() for x in fa_s]
        for a, fa in zip(sorted(a_s, key=lambda k: k.get_id()), sorted(fa_s, key=lambda k: k.get_id())):
            if a.get_id() in aids and fa.get_id() in faids:
                masked_atoms1.append(a)
                masked_atoms2.append(fa)
    return masked_atoms1, masked_atoms2


def get_atom_positions(atoms):
    """
    :arg atoms a dict of atoms
    """
    ret = {}
    for k, v in zip(atoms.keys(), atoms.values()):
        ret[k] = np.concatenate(list(map(lambda k: k.get_coord(), list(itertools.chain.from_iterable(v))))).reshape(-1, 3)
    return ret


def to_lddt_format(atoms):
    return np.stack([atom.get_coord() for atom in atoms]).reshape(1, -1, 3)


def calc_exp_distances(proteins):
    """
    :arg proteins is in pairs
    """

    ret = {}
    for p1, p2 in proteins:
        alignment_dict = align(proteins)
        all_residues = get_residues(list(itertools.chain.from_iterable(proteins)))

        a1, non_res1, front_back1, middle1 = alignment_dict[p1]
        a2, non_res2, front_back2, middle2 = alignment_dict[p2]
        res1 = iter(np.array(all_residues[p1], dtype=object)[non_res1])
        res2 = iter(np.array(all_residues[p2], dtype=object)[non_res2])
        residues1 = []
        residues2 = []

        for r1, r2 in zip(a1.seq, a2.seq):
            if r1 == '-':
                residues1.append(None)
                residues2.append(next(res2))
            elif r2 == '-':
                residues1.append(next(res1))
                residues2.append(None)
            else:
                residues1.append(next(res1))
                residues2.append(next(res2))

        front_back = np.logical_and(front_back1, front_back2)
        middle = np.logical_and(middle1, middle2)
        mask = np.logical_and(front_back, middle)

        residues1 = np.array(residues1, dtype=object)[mask]
        residues2 = np.array(residues2, dtype=object)[mask]


        # all_folded_residues = get_residues(list(itertools.chain.from_iterable(proteins)), folded=True)
        all_folded_residues = get_residues(proteins[0], folded=True)
        fres1 = iter(all_folded_residues[p1])
        fres2 = iter(all_folded_residues[p2])

        fresidues1 = []
        fresidues2 = []

        for fb in front_back:
            if not fb:
                fresidues1.append(None)
                fresidues2.append(None)
            else:
                fresidues1.append(next(fres1))
                fresidues2.append(next(fres2))
        fresidues1 = np.array(fresidues1, dtype=object)[mask]
        fresidues2 = np.array(fresidues2, dtype=object)[mask]

        atoms1 = list(map(lambda k: list(k.get_atoms()), residues1))
        atoms2 = list(map(lambda k: list(k.get_atoms()), residues2))
        fatoms1 = list(map(lambda k: list(k.get_atoms()), fresidues1))
        fatoms2 = list(map(lambda k: list(k.get_atoms()), fresidues2))

        catoms1, cfatoms1 = trim_atoms(atoms1, fatoms1)
        catoms2, cfatoms2 = trim_atoms(atoms2, fatoms2)

        # RMS and LDDT = [Folded v Regular WT, Folded v Regular MUT, Exp WT v Exp MUT, Folded WT v Folded MUT]
        RMS = []
        LDDT = []

        super_imposer = Superimposer()
        super_imposer.set_atoms(catoms1, cfatoms1)
        # super_imposer.apply(cfatoms1)
        lddt_ = lddt(to_lddt_format(catoms1), to_lddt_format(cfatoms1), np.ones([1, len(catoms1), 1]))
        RMS.append(super_imposer.rms); LDDT.append(lddt_)

        print('Folded v Regular WT {}, {}: RMS = {} lddt = {}'.format(p1, p2, super_imposer.rms, lddt_))

        super_imposer.set_atoms(catoms2, cfatoms2)
        # super_imposer.apply(cfatoms2)
        lddt_ = lddt(to_lddt_format(catoms2), to_lddt_format(cfatoms2), np.ones([1, len(catoms2), 1]))
        RMS.append(super_imposer.rms); LDDT.append(lddt_)


        print('Folded v Regular MUT {}, {}: RMS = {} lddt = {}'.format(p1, p2, super_imposer.rms, lddt_))

        catoms1, catoms2 = trim_atoms(atoms1, atoms2)
        cfatoms1, cfatoms2 = trim_atoms(fatoms1, fatoms2)

        super_imposer.set_atoms(catoms1, catoms2)
        # super_imposer.apply(catoms2)
        lddt_ = lddt(to_lddt_format(catoms1), to_lddt_format(catoms2), np.ones([1, len(catoms1), 1]))
        RMS.append(super_imposer.rms); LDDT.append(lddt_)

        print('Exp WT v Exp MUT {}, {}: RMS = {} lddt = {}'.format(p1, p2, super_imposer.rms, lddt_))

        super_imposer.set_atoms(cfatoms1, cfatoms2)
        # super_imposer.apply(cfatoms2)
        lddt_ = lddt(to_lddt_format(cfatoms1), to_lddt_format(cfatoms2), np.ones([1, len(cfatoms1), 1]))
        RMS.append(super_imposer.rms); LDDT.append(lddt_)

        print('Folded WT v Folded MUT {}, {}: RMS = {} lddt = {}'.format(p1, p2, super_imposer.rms, lddt_))
        ret[p1] = {"RMS": RMS, "LDDT": LDDT}
        ret[p2] = {"RMS": RMS, "LDDT": LDDT}

    return ret


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
    distance_dict = calc_exp_distances(pairings)
