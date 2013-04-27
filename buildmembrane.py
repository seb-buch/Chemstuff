#!/usr/bin/python
from __future__ import division
__author__ = 'Sébastien Buchoux <sebastien.buchoux@gmail.com'
__version__ = '0.1.0'

import argparse
import glob
import linecache
import numpy
import os
import shutil
import subprocess
import sys
import tempfile

from util import print_error

parser = argparse.ArgumentParser()

# Populate parser
parser.add_argument('--version', action='version', version=('%(prog)s ' + __version__))
parser.add_argument("--dim", nargs=2, default=[8, 8], type=int,
                    help="Size of the membrane (default is 8 by 8)")

parser.add_argument("--lipids", nargs="*", default=["dmpc"], type=str,
                    help="lipids to use to build membrane")

parser.add_argument("--lipcomp", nargs="*", default=[1], type=float,
                    help="lipid composition (must be same length as the --lipids argument")

parser.add_argument("--apl", default=0.65, type=float,
                    help="area per lipids in nm²")

parser.add_argument("--gap", default=0.5, type=float,
                    help="gap between box")

parser.add_argument("--list", "-l", help="show list of lipids available",
                    action="store_true",
                    default=False)

parser.add_argument("--output", "-o", help="output structure file",
                    default="confout.gro",
                    type=str)

parser.add_argument("--keeptemp", help="keep temporary files",
                    action="store_true",
                    default=False)
# Parse arguments
args = parser.parse_args()

# shortcuts for arguments
xdim, ydim = args.dim
libdir = "%s/lib" % os.path.dirname(os.path.realpath(__file__))
lipids = args.lipids
lipcomp = args.lipcomp
confout = os.path.relpath(args.output)
keeptemp = args.keeptemp
apl = args.apl
gap = args.gap
show_list = args.list
inc = numpy.sqrt(apl)

# Useful functions
def load_lipid_lib(path=libdir):
    """
    Load and check available lipid library

    :param path: File path where lipid topologies are
    """
    global lipid_lib, libdir

    path = os.path.realpath(path)
    sys.stdout.write("Loading library from '%s'... " % path)

    for fname in glob.glob("%s/*.gro" % path):
        name = os.path.basename(fname)[:-4]
        resname = linecache.getline(fname, 3)[5:9].strip()
        lipid_lib[name] = (fname, resname)

    print("%i lipids found" % len(lipid_lib))

    libdir = path


def run_cmd(cmd):
    popen = subprocess.Popen(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    output = popen.communicate()

    if popen.returncode != 0:
        print("FATAL: subprocess command ended badly (returned %i)" % popen.returncode)
        print("-> command: %s" % cmd)
        print("-> output:")
        print(output[0].decode())

        try:
            print("\nTemporary files are available here: %s" % tmpdir)
        except NameError:
            pass
        sys.exit(1)


def progress_bar(step=0, max_step=100, title="Work in progress", size=20):
    nb_sharp = int(round(step * size / max_step))
    nb_flat = size - nb_sharp

    sys.stdout.write("\r%s [%s%s] " % (title, "#" * nb_sharp,
                                       "-" * nb_flat))

    if step < max_step:
        num_size = len("%d" % max_step)
        progress = "%%%dd / %%%dd" % (num_size, num_size)
        sys.stdout.write(progress % (step, max_step))
    else:
        sys.stdout.write("Done!           \n")

    sys.stdout.flush()

# Load lipid lib
lipid_lib = {}
load_lipid_lib()

# List available lipids and exit
if show_list:
    print("Available lipids:")
    for key in lipid_lib.keys():
        print("  %s" % key)
    sys.exit(0)

# Check if number of lipids and composition are consistent
if len(lipids) != len(lipcomp):
    print_error("number of lipids and composition are not consistent")
    sys.exit(1)

for lipid in lipids:
    if lipid not in lipid_lib:
        print("ERROR: %s lipid not found in library" % lipid)
        print("Available lipids:")
        for key in lipid_lib.keys():
            print("  %s" % key)
        sys.exit(1)

num_comp = len(lipids)
composition = []
for i in range(num_comp):
    composition.append([])


num_lipid = xdim * ydim * 2


print("Creating a %d x %d lipid bilayers (total: %d lipids):" % (xdim, ydim, num_lipid))
print("-> Estimated membrane area: %.3f nm²" % (xdim * ydim * apl))


# Ajust composition
lipcomp = numpy.array(lipcomp)
lipcomp /= lipcomp.sum()

lipcomp_pprint = ""
for i, comp in enumerate(lipcomp):
    lipcomp_pprint += "%s(%.2f%%)" % (lipids[i], comp * 100)
    if i < (len(lipcomp) - 1):
        lipcomp_pprint += ", "

print("-> Desired lipid composition: %s" % lipcomp_pprint)

lipcomp_ajusted = numpy.rint(lipcomp * num_lipid).astype(int)
lipcomp_ajusted[0] -= lipcomp_ajusted.sum() - num_lipid # make sure the number of lipids is respected

lipcomp_pprint = ""
for i, comp in enumerate(lipcomp_ajusted):
    lipcomp_pprint += "%d %s(%.2f%%)" % (comp, lipids[i], comp / num_lipid * 100)
    if i < (len(lipcomp) - 1):
        lipcomp_pprint += ", "

print("-> Ajusted lipid composition: %s" % lipcomp_pprint)

# Randomize lipid composition
lipid_indices = numpy.arange(num_lipid)
numpy.random.shuffle(lipid_indices)


# Creating temp dir
tmpdir = tempfile.mkdtemp()


begin = 0
end = 0
resnr = 1
cmd_dict = {"fin": "",
            "fout": "",
            "resnr": resnr,
            "rotangle": 0,
            "trans_x": 0,
            "trans_y": 0,
            "trans_z": 0}


tmpout = "%s/tmpconf.gro" % tmpdir
tmpoutf = open(tmpout, "w")
tmpoutf.write("Gro file generate by %s v%s\n" % (parser.prog, __version__))
num_atoms_pos = tmpoutf.tell()
tmpoutf.write("              \n")
num_atoms = 0
try:
    for i, lipid in enumerate(lipids):
        size = lipcomp_ajusted[i]
        end += size

        # create copies of the lipids
        run_cmd("editconf -quiet -f %s -o %s -d 0" % (lipid_lib[lipid][0], "%s/lip%i_up.gro" % (tmpdir, i)))
        run_cmd("editconf -quiet -f %s -o %s -rotate 180 0 0 -d 0" % (lipid_lib[lipid][0], "%s/lip%i_down.gro" %
                                                                                           (tmpdir, i)))
        # get box height
        for step, j in enumerate(lipid_indices[begin:end]):
            progress_bar(step, end-begin, "-> Putting %i %s on randomized grid" % (size, lipid))
            z = -1 if j < 64 else 1
            if z > 0:
                j -= 64
                cmd_dict["fin"] = "%s/lip%i_up.gro" % (tmpdir, i)
            else:
                cmd_dict["fin"] = "%s/lip%i_down.gro" % (tmpdir, i)

            x = j % xdim
            y = j // xdim

            cmd_dict["fout"] = "%s/%i.gro" % (tmpdir, resnr)
            cmd_dict["resnr"] = resnr
            cmd_dict["trans_z"] = z * gap / 2
            cmd_dict["trans_x"] = x * inc
            cmd_dict["trans_y"] = y * inc

            cmd = "editconf -quiet -f %(fin)s -o %(fout)s -resnr %(resnr)i \
            -translate %(trans_x)f %(trans_y)f %(trans_z)f" % cmd_dict

            # Run editconf
            run_cmd(cmd)

            # Concatenate to tmp structure
            fo = open(cmd_dict["fout"])
            num_atoms_tmp = 0
            atom_count = 0
            for lino, line in enumerate(fo):
                if lino == 1:
                    num_atoms_tmp = int(line.strip())
                elif lino > 1:
                    if atom_count >= num_atoms_tmp:
                        break
                    else:
                        atom_count += 1
                        tmpoutf.write(line)
            fo.close()
            num_atoms += num_atoms_tmp

            # Update progress status
            progress_bar(step+1, end-begin, "-> Putting %i %s on randomized grid" % (size, lipid))




            # increment residue number
            resnr += 1
        begin += end

    tmpoutf.write("\n")
    tmpoutf.seek(num_atoms_pos)
    tmpoutf.write("%i" % num_atoms)
    tmpoutf.close()
    run_cmd("editconf -quiet -d 0 -f %s -o %s" % (tmpout, confout))

    print("\nMembrane saved to %s" % confout)
except KeyboardInterrupt:
    print("\n\nKeyboard interrupt!")

if keeptemp:
    print("Temporary files are still available in dir: " % tmpdir)
else:
    shutil.rmtree(tmpdir)