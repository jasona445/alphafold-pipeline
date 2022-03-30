from utils.utils import *
import itertools


def download_files(proteins):
    """
    Takes in a list of proteins in PDB entry format (e.g. 1bz4)

    Downloads .ent's from PDB via curl
    """

    script = rel_to_abs('scripts/batch_download.sh', os.path.dirname(__file__))
    tmp = rel_to_abs('tmp/proteins.txt', os.path.dirname(__file__))
    output = rel_to_abs('data/experimental-structures', os.path.dirname(__file__))

    with open(tmp, 'w') as f:
        f.write(','.join(proteins))

    system_call("source \"{}\" -f \"{}\" -o \"{}\" -p".format(script, tmp, output))
    system_call("gunzip -r \"{}\"".format(output))

    # Clean up zips if gunzip fails
    contents = os.listdir(output)
    for item in contents:
        if item.endswith(".pdb.gz"):
            os.remove(os.path.join(output, item))
    os.remove(tmp)


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
    download_files(itertools.chain.from_iterable(pairings))







