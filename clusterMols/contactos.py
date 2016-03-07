#!/usr/bin/env python

# Contactos.py - calculate similarities of molecular docking results

# Version 1.1

# Copyright 2007, 2008 Mikko Huhtala and Santeri Puranen

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.


import sys, math, itertools, operator, copy, os
from optparse import OptionParser

class Atom:

    def __init__(self, x = 0.0, y = 0.0, z = 0.0, index=False, type = ''):

        self.x = x
        self.y = y
        self.z = z
        self.index = index
        self.type = type

    def distance(self, other):

        return math.sqrt((self.x - other.x) ** 2 + \
                         (self.y - other.y) ** 2 + \
                         (self.z - other.z) ** 2)

    def squaredDistance(self, other):

        return (self.x - other.x) ** 2 + \
               (self.y - other.y) ** 2 + \
               (self.z - other.z) ** 2


class Conformation:

    def __init__(self):

        self.name = ''
        self.filename = ''
        self.protein = []
        self.xmin = 0
        self.ymin = 0
        self.zmin = 0
        self.xmax = 0
        self.ymax = 0
        self.zmax = 0
        self.ligand = []


class ContactDescriptor:

    def __init__(self, name = '', filename = ''):

        self.contacts = []
        self.name = name
        self.filename = filename

    def check(self, other):

        if len(self.contacts) != len(other.contacts):
            sys.stderr.write('\nError: descriptors of different length for %s and %s.\n\n' % \
                             (self.name, other.name))
            sys.exit(1)


    def set_brute_contacts(self, C, cutoff):

        self.contacts = []
        squaredCutoff = cutoff**2
        for i in C.ligand:
            for j in C.protein:
                if i.squaredDistance(j) <= squaredCutoff:
                    self.contacts.append(1)
                else:
                    self.contacts.append(0)


    def set_contacts(self, C, cutoff):

        self.contacts = [ 0 for i in xrange(len(C.protein) * len(C.ligand)) ]
        grid = GridSearch(C.protein, cutoff, [C.xmin, C.xmax, \
                                              C.ymin, C.ymax, \
                                              C.zmin, C.zmax])
        k = len(C.protein)
        base = 0
        for i in C.ligand:
            for j in grid.getWithinCutoff(i, cutoff):
                self.contacts[base+j.index] = 1
            base += k


    def similarity(self, other):

        self.check(other)
        shared = float(sum(itertools.imap(operator.and_, self.contacts, other.contacts)))
        maxcontacts = float(sum(itertools.imap(operator.or_, self.contacts, other.contacts)))
        return shared / maxcontacts


class TypeContactDescriptor(ContactDescriptor):

    def set_contacts(self, C, cutoff):

        self.contacts = [ [] for i in xrange(len(C.protein)) ]
        grid = GridSearch(C.protein, cutoff)
        for i in C.ligand:
            for j in grid.getWithinCutoff(i, cutoff):
                self.contacts[j.index].append(i.type)

    def atom_type_score(self, type):
        if type == 'H':
            return 0.5
        elif type[0] != 'C':
            return 3.0
        else:
            return 1.0

    def similarity(self, other):

        self.check(other)
        score = 0.0
        maxscore = 0.0
        for i in xrange(len(self.contacts)):
            contacts_in_other = copy.copy(other.contacts[i])
            union = []
            for j in self.contacts[i]:
                union.append(j)
                if j in contacts_in_other:
                    score += self.atom_type_score(j)
                    contacts_in_other.remove(j)
            union += contacts_in_other
            for j in union:
                maxscore += self.atom_type_score(j)
        if maxscore > 0:
            return score / maxscore
        else:
            return 0


class GridSearch:
    # A simple gridsearch algorithm for nearest neighbor search.

    def __init__( self, atomList, cutoff, bounds=[] ):

        # 'bounds' should contain the coordinate bounds (the
        # minimum and maximum values along each coordinate axis
        # separately) of the molecules being put into the search grid:
        #
        # bounds = [xmin,xmax,ymin,ymax,zmin,zmax]
        #
        # The default, empty list, causes the routine to find the bounds
        # by itself, at the cost of looping through all coordinates.
        self.spacing = cutoff
        if not bounds:
            bounds = self.findBounds( atomList )
        self.buildGrid( atomList, cutoff, bounds )

    def findBounds( self, atomList ):
        atom = atomList[0]
        xmin = xmax = atom.x
        ymin = ymax = atom.y
        zmin = zmax = atom.z
        for atom in atomList[1:]:
            xmin = min( xmin, atom.x )
            xmax = max( xmax, atom.x )
            ymin = min( ymin, atom.y )
            ymax = max( ymax, atom.y )
            zmin = min( zmin, atom.z )
            zmax = max( zmax, atom.z )

        return [xmin,xmax,ymin,ymax,zmin,zmax]

    def buildGrid( self, atomList, cutoff, bounds ):
        # Grid dimensions. Make large enough to fit the whole search structure.
        # Beware of overflow, though (x*y*z larger than max. for numeric type).
        self.dimensions = [ int( math.ceil( abs(bounds[1]-bounds[0]) / self.spacing ) ), \
                            int( math.ceil( abs(bounds[3]-bounds[2]) / self.spacing ) ), \
                            int( math.ceil( abs(bounds[5]-bounds[4]) / self.spacing ) ) ]

        # Coordinates of the grid center
        self.xcenter = (bounds[0]+bounds[1])*.5
        self.ycenter = (bounds[2]+bounds[3])*.5
        self.zcenter = (bounds[4]+bounds[5])*.5

        # constants needed for determining gridcell index
        #self.xfac = 1.0
        self.yfac = self.dimensions[0]
        self.zfac = self.dimensions[0]*self.dimensions[1]

        self.searchGrid = {}

        # calculate gridcell index for each 'atom' and store in dictionary (with
        # index as key and a list of corresponding 'atom' instances as value)
        for i in atomList:
            self.searchGrid.setdefault( self.getIndex( i ), [] ).append( i )

    def getIndex( self, atom ):
        # return the index of the gridcell that corresponds
        # to the coordinates of 'atom'
        return int( ( (atom.x-self.xcenter) // self.spacing ) + \
                    ( (atom.y-self.ycenter) // self.spacing ) * self.yfac + \
                    ( (atom.z-self.zcenter) // self.spacing ) * self.zfac )


    def getWithinCutoff( self, atom, cutoff ):
        # return all atom instances that are within 'cutoff' from 'atom'

        # get gridcell index of 'atom'
        i = self.getIndex( atom )
        y = self.yfac
        z = self.zfac

        atoms = []
        # loop over indices of gridcells that surround 'atom' (including
        # gridcell index of 'atom' itself; 3^3 = 27 indices in total)
        for index in [ i-1-y-z, i-y-z, i+1-y-z, i-1-z, i-z, i+1-z, i-1+y-z, i+y-z, i+1+y-z, \
                       i-1-y, i-y, i+1-y, i-1, i, i+1, i-1+y, i+y, i+1+y, \
                       i-1-y+z, i-y+z, i+1-y+z, i-1+z, i+z, i+1+z, i-1+y+z, i+y+z, i+1+y+z ]:
            atoms.extend( self.searchGrid.get( index, [] ) )

        hits = []
        # square cutoff distance in order to save call
        # to math.sqrt() for each distance calculation
        cutoffSquared = cutoff**2
        for a in atoms:
            if atom.squaredDistance( a ) <= cutoffSquared:
                hits.append( a )

        return hits


def remove_noninformative_pos(descriptors):

    size = len(descriptors[0].contacts)
    for i in descriptors[1:]:
        if len(i.contacts) != size:
            sys.stderr.write('Error: descriptors are of different sizes. Use option \'-t\' if the input includes different ligands.\n')
            sys.exit(1)

    # are these lists of ints or lists of lists?
    type_descriptor = 0
    if type(descriptors[0].contacts[0]) != type(1):
        type_descriptor = 1


    # create a list of informative positions
    if type_descriptor:
        pseudo_desc = [ [] for i in xrange(len(descriptors[0].contacts)) ]
    else:
        pseudo_desc = [ 0 for i in xrange(len(descriptors[0].contacts)) ]

    desc_count = len(descriptors)
    for d in descriptors:
        desc_count -= 1
        pseudo_desc = map(operator.add, pseudo_desc, d.contacts)
        sys.stderr.write("\rRemoving non-informative positions from descriptors: %d " % desc_count)

    to_keep = []
    for index,val in enumerate(pseudo_desc):
        if val:
            to_keep.append(index)

    # create a new list of descriptors and copy only informative
    # positions into it
    compacted = []
    for i in descriptors:
        if type_descriptor:
            D = TypeContactDescriptor(name = i.name, filename = i.filename)
        else:
            D = ContactDescriptor(name = i.name, filename = i.filename)
        compacted.append(D)

    for cj, dj in itertools.izip( compacted, descriptors ):
        for i in to_keep:
            cj.contacts.append(dj.contacts[i])

    sys.stderr.write("\rRemoving non-informative positions from descriptors: done")

    return compacted


def get_next_from_mol2(F, receptor = False, ligresid = 'UNK'):

    if not receptor:
        new_conformation = Conformation()
    else:
        new_conformation = copy.deepcopy(receptor)
    new_conformation.filename = F.name
    I = F.readline()
    index = itertools.count()
    read = 0
    name = 0
    first_atom = 1
    while I:

        if not I.strip():
            I = F.readline()
            continue

        if I[:13] == "@<TRIPOS>ATOM":
            read = 1
            I = F.readline()
            continue

        if read and I[:9] == "@<TRIPOS>":
            break

        if I[:17] == "@<TRIPOS>MOLECULE":
            name = 1
            I = F.readline()
            continue

        if name:
            new_conformation.name = I.strip().replace(' ', '_')
            name = 0

        if read:
            if I.strip():
                fields = I.split()
                x = float(fields[2])
                y = float(fields[3])
                z = float(fields[4])

                if first_atom:
                    first_atom = 0
                    new_conformation.xmin = x
                    new_conformation.ymin = y
                    new_conformation.zmin = z
                    new_conformation.xmax = x
                    new_conformation.ymax = y
                    new_conformation.zmax = z

                new_conformation.xmin = min(x, new_conformation.xmin)
                new_conformation.ymin = min(y, new_conformation.ymin)
                new_conformation.zmin = min(z, new_conformation.zmin)
                new_conformation.xmax = max(x, new_conformation.xmax)
                new_conformation.ymax = max(y, new_conformation.ymax)
                new_conformation.zmax = max(z, new_conformation.zmax)

                type = fields[5]
                res_type = fields[7]
                if receptor:
                    new_conformation.ligand.append(Atom(x, y, z, 0, type))
                elif res_type[:3] == ligresid and not receptor:
                    new_conformation.ligand.append(Atom(x, y, z, 0, type))
                else:
                    new_conformation.protein.append(Atom(x, y, z, index.next(), type))

        I = F.readline()

    return new_conformation



def parse_commandline():

    usage = "usage: %prog [options] mol2file1 [mol2file2 ...] \n\nOne or more MOL2 files can be given. Each file may contain one or\nmore models. Residue UNK (default, see option -l)  is the ligand\nand everything else is the protein. The order of lines must be consistent\nacross models, except for ligands when using -t (clustering different\nligands)."

    parser = OptionParser(usage=usage)

    parser.add_option("-c", "--cut-off", type="float", action = "store", dest = "cutoff", default = 3.0, \
                      help = "Cut-off distance for atom-atom contacts in Angstroms. The default is 3.0.")

    parser.add_option("-o", "--output", action = "store", dest = "output", default = "contactos_out", \
                      help = "Prefix for output file names. Will create '.log', '.mcl_pairs', '.similarity_matrix' and optionally '.mcl_out' and the directory '.mcl_clusters'. Default prefix is 'contactos_out'.")

    parser.add_option("-t", "--types", action = "store_true", dest = "types", default = False, \
                      help = "Generate descriptors based only on ligand atom types in contact with protein. This is a lot less sensitive than the default algorithm, but it allows a set of different ligand molecules in the input. Only the receptor needs to be the same in all of the docking results. ")

    parser.add_option("-l", "--ligand-resid", action = "store", dest = "ligresid", default = "UNK", \
                      help = """Residue identifier of the ligand. Default is 'UNK'. This is used in the default situation, where each receptor conformation is accompanied by a ligand in the same file. This option is meaningless if option -r is used.""")

    parser.add_option("-r", "--receptor", action = "store", dest = "receptor", default = "", \
                      help = "Specify a receptor structure file. By default, it is assumed that each input structure file contains a receptor conformation and a ligand. If --receptor is given, only one single receptor conformation is read from the given file instead, and the rest of the mol2 input files are assumed to be ligand poses.")

    parser.add_option("-m", "--run-mcl", action = "store_true", dest = "mcl", default = False, \
                      help = "Run MCL clustering on the generated similarity matrix. The executable 'mcl' must be in the search path. The input structure files are automatically symlinked into directories corresponding to the clusters found by MCL. See also --mcl-inflation.")

    parser.add_option("-I", "--mcl-inflation", type="float", action = "store", dest = "inflation", default = 2.0, \
                      help = "The inflation parameter passed to MCL. Affects clustering granularity. Must be between 0.1 and 30.0. The default is 2.0. Useful only with --run-mcl. Experimentation is likely needed to find a good value. Too coarse a granularity setting will result in the largest cluster containing unrelated ligand poses and too fine granularity will result in many clusters containing only one member. Larger values result in more clusters.")

    return parser.parse_args()


if __name__ == "__main__":


    (options, args) = parse_commandline()
    descriptors = []
    filename_mapping = {}
    filenames = args

    if not filenames:
        sys.stderr.write('Get help with -h.\n')
        sys.exit(1)

    out_log = file(options.output + '.log', 'w')
    out_pairs = file(options.output + '.mcl_pairs', 'w')
    out_matrix = file(options.output + '.similarity_matrix', 'w')

    out_log.write('\n:: contactos.py\n\nAtom-atom distance cut-off: %5.2f\n' % options.cutoff)

    if options.types:
        out_log.write( "\nContact descriptor: protein atom - ligand atom type (can cluster different ligands)\n" )
    else:
        out_log.write( "\nContact descriptor: protein atom - ligand atom (all ligands must be the same)\n" )


    # if receptor is in a separate file, get it now
    receptor = False
    if options.receptor:
        F = file(options.receptor, 'r')
        receptor = get_next_from_mol2(F, ligresid = options.ligresid)
        out_log.write( "\nUsing a single protein conformation from %s\n" % options.receptor )
    else:
        out_log.write( "\nUsing multiple protein conformations, ligand residue is %s\n" % options.ligresid )

    # generate descriptors for all conformations
    count = 0
    for i in filenames:
        F = file(i, 'r')
        C = get_next_from_mol2(F, receptor, ligresid = options.ligresid)
        while C.protein and C.ligand:
            count += 1
            sys.stderr.write("\rGenerating descriptors. Conformation no %d : %s                    " % (count, C.name))
            if options.types:
                D = TypeContactDescriptor(C.name, C.filename)
            else:
                D = ContactDescriptor(C.name, C.filename)
            D.set_contacts(C, options.cutoff)
            descriptors.append(D)
            filename_mapping[C.name] = C.filename
            C = get_next_from_mol2(F, ligresid = options.ligresid)
    sys.stderr.write('\n')

    if not options.types:
        out_log.write('\nLigand atom - protein atom pairs: %d\n' % len(descriptors[0].contacts))
    else:
        out_log.write('\nNumber of protein atoms: %d\n' % len(descriptors[0].contacts))

    # get rid of non-informative positions in descriptors
    descriptors = remove_noninformative_pos(descriptors)
    sys.stderr.write('\n')
    out_log.write('\nInformative positions: %d\n' % len(descriptors[0].contacts))


    # calculate similarities
    out_log.write('\nNumber of ligand-protein contacts found in models:\n\n')
    size = len(descriptors)
    matrix = []
    count = (size * (size + 1)) / 2
    k = 0
    for i in xrange(size):
        matrix.append([])
        k += 1
        for j in range(k):
            sys.stderr.write("\rCalculating similarities: %d " % count)
            matrix[i].append(descriptors[i].similarity(descriptors[j]))
            count -= 1
            if i == j:
                if not options.types:
                    out_log.write('%s %d\n' % (descriptors[i].name, sum(descriptors[i].contacts)))
                else:
                    total = 0
                    for m in descriptors[i].contacts:
                        total += len(m)
                    out_log.write('%s %d\n' % (descriptors[i].name, total))


    sys.stderr.write("\rCalculating similarities: done\n")


    # write similarity matrix
    out_matrix.write("%d\n" % size)
    for i in descriptors:
        out_matrix.write("%s " % i.name)
    out_matrix.write('\n')
    for i in matrix:
        for j in i:
            out_matrix.write("%6.5f " % j)
        out_matrix.write('\n')
    out_matrix.close()


    # write list of pairs, in both directions for MCL
    k = 0
    for i in xrange(size):
        k += 1
        for j in xrange(k):
            if i == j:
                continue
            out_pairs.write("%s %s %6.5f\n" % (descriptors[i].name,\
                                               descriptors[j].name,\
                                               matrix[i][j]))
            out_pairs.write("%s %s %6.5f\n" % (descriptors[j].name,\
                                               descriptors[i].name,\
                                               matrix[i][j]))
    out_pairs.close()

    # if requested, run MCL and parse its output
    if options.mcl:
        sys.stderr.write("Running MCL clustering.\n\n")
        sys.stderr.write("----------------------------------------------------------------------\n\n")

        mcl_out_file = options.output + '.mcl_out'

        if os.path.isfile(mcl_out_file):
            os.remove(mcl_out_file)

        command = 'mcl %s --abc -o %s -I %f' % (options.output + '.mcl_pairs', mcl_out_file, options.inflation)
        sys.stderr.write(command + '\n')
        out_log.write('\nRunning MCL with the following command\n\n' + command + '\n')

        os.system(command)

        if not os.path.isfile(mcl_out_file):
            sys.stderr.write("\n\nNo MCL output file found. Running MCL seems to have failed.\n")
            sys.exit()

        sys.stderr.write("----------------------------------------------------------------------\n\n")

        mcl_lines = file(mcl_out_file, 'r').readlines()
        mcl_clusters = []
        for i in mcl_lines:
            mcl_clusters.append(i.split())

        os.mkdir(options.output + '.mcl_clusters')
        os.chdir(options.output + '.mcl_clusters')

        for i,j in enumerate(mcl_clusters):
            dirname = str( i + 1 )
            os.mkdir(dirname)
            os.chdir(dirname)
            for k in j:
                target = os.path.basename(filename_mapping[k])
                if filename_mapping[k][0] != '/':
                    source = '../../' + filename_mapping[k]
                else:
                    source = filename_mapping[k]
                if not os.path.exists(target):
                    os.symlink(source, target)
            os.chdir('..')

        os.chdir('..')
