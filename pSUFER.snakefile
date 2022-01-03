#!/usr/bin/env python3
import logging
from collections import OrderedDict, Counter
import re
import numpy as np
from copy import deepcopy
import os
from functools import reduce
from Bio.SeqUtils import seq1, seq3
from Bio.PDB.PDBParser import PDBParser
import Bio.PDB.Structure
import gzip



__author__ = 'Rosalie Lipsh'

LGR = logging.getLogger(__name__)

def parse_pdb_file(pdb_path):
     """Reads a pdb file and return a structure object"""
     if isinstance(pdb_path, Bio.PDB.Structure.Structure):
         return pdb_path
     if pdb_path.endswith('.gz'):
         with gzip.open(pdb_path, 'rb') as f:
             file_content = f.read().decode('utf-8')
             tmp = 'tmp_'+pdb_path.split('/')[-1]
             open(tmp, 'w').write(file_content)
             pdb = PDBParser().get_structure('', tmp)
             os.remove(tmp)
     else:
         pdb = PDBParser().get_structure('',pdb_path)
     return pdb

def check_pdb(pdb):
    """If the pdb is a path, parses it into a pdb object"""
    if not isinstance(pdb, Bio.PDB.Structure.Structure):
        pdb = parse_pdb_file(pdb)
    return pdb

def is_hetatm(residue):
    assert type(residue) == Bio.PDB.Residue.Residue
    return residue.id[0].replace(' ', '')

def get_sequence(pdb, chain=''):
  """
  :param chain: If nor specified, reurns the sequence of all chains
  :param pdb:
  :return:
  """
  pdb = check_pdb(pdb)
  chain = chain.upper()
  seqs = OrderedDict([(c.id, ''.join([seq1(r.resname) for r in c if not is_hetatm(r)])) for c in pdb.get_chains()])
  if chain:
      seq = seqs[chain]
  else:
      seq = ''.join([s for c, s in seqs.items()])
  return seq


class ResFile:
    """Currently works only on resfiles from type PIKAA"""
    def __init__(self, path='', pdb='', pssm=''):
        """
        :param pdb: path to pdb the resfile belonges to (position nubering is
        by it)
        :param pssm: WARNING use a pdb renumbered as the pssm
        """
        self.path = path
        self.resfile = OrderedDict()
        if self.path:
            self.parse_res_file(self.path)
        self.rearranged = False
        self.set_pdb(pdb)
        self.set_pssm(pssm)

    def __update(self, pos, allowed, union=True, force=False):
        """
        Upadtes the allowed residues in pos
        :param force: if True, will copy the residues as is, no checks
        """
        if not hasattr(self, 'resfile'):
            self.resfile = OrderedDict()
        if force:
            self.resfile[pos] = allowed
            return
        allowed = allowed.upper()
        org = allowed
        allowed = re.sub('[0-9BJOUXZ]', '', allowed)
        removed = set(org) - set(allowed)
        if removed:
            LGR.warning(('The following characters were removed from '
                            'position {}: {}'.format(pos, ','.join(removed))))
        if pos in self.resfile.keys() and union:
            LGR.warning(('{} is in both resfiles. taking the union of '
                           'the residues').format(pos))
            allowed = ''.join(set(self.resfile[pos]) | set(allowed))
            self.rearranged = False
        self.resfile[pos] = allowed
        LGR.debug('Updated resfile: {}={}'.format(pos, allowed))

    def __sort(self):
        """
        Sorts the resfile by chain and position
        """
        self.resfile = OrderedDict(sorted(self.resfile.items(),
                                          key=lambda x: (x[0][-1],
                                                         int(x[0][:-1]))))

    def parse_res_file(self, resfile_path):
        """
        :param resfile_path:
        :return: list of tuples (residueCHAIN, the allowed residues as a string)
        """
        LGR.info('Parsing the resfile ' + resfile_path)
        resfile = open(resfile_path, 'r').readlines()
        resfile = OrderedDict((''.join(line.split()[:2]),
                               line.split()[-1].upper())
                              for line in resfile if 'PIKAA' in line)
        for pos, allowed in resfile.items():
            self.__update(pos, allowed)
        self.__sort()

    def init_from_dict(self, res_dict):
        """
        :param res_dict: [position] = string of AAs
        """
        for pos, allowed in res_dict.items():
            pos = str(pos)
            self.__update(pos, allowed)
        self.__sort()
        self.rearranged = False
        if self.pdb:
            self.rearrange_resfile()

    def set_pdb(self, pdb_path):
        self.pdb = parse_pdb_file(pdb_path) if pdb_path else ''

    def set_pssm(self, pssm_path):
        # self.pssm = PSSM(name='', pssm_file=pssm_path) if pssm_path else None
        self.pssm = PSSM(pssm_path) if pssm_path else None

    def __len__(self):
        return len(self.resfile)

    def __eq__(self, other):
        if sorted(self.resfile.keys()) != sorted(other.resfile.keys()):
            return False
        for pos in self.resfile.keys():
            if not sorted(self.resfile[pos]) == sorted(other.resfile[pos]):
                return False
        return True

    def __str__(self):
        return '\n'.join(['{:<6} {}'.format(k, v) for k, v in
                          self.resfile.items()])

    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        new_rf = deepcopy(self)
        if self.pdb and other.pdb:
            if get_sequence(self.pdb) != get_sequence(other.pdb):
                raise ValueError(('concatenated resfiles should come from the '
                                  'same protein:\n {}\n{}').format(self.pdb,
                                                                   other.pdb))

        new_rf.init_from_dict(other.resfile)
        return new_rf

    def __sub__(self, other):
        """Takes the difference between the allowed residues BUT leaves the WT
        identity
        """
        if not self.pdb:
            raise ValueError('Cannot substract without the WT')
        new_rf = deepcopy(self)
        new_rf.rearrange_resfile()
        for pos, allowed in new_rf:
            new_rf.resfile[pos] = ''.join(set(allowed) - set(other[pos]))
            if allowed[0] not in new_rf.resfile[pos]:
                new_rf.resfile[pos] += allowed[0]
        new_rf.rearrange_resfile()
        return new_rf

    def __iter__(self):
        for residue, allowed in self.resfile.items():
            yield residue, allowed

    def __getitem__(self, position):
        position = str(position)
        return self.resfile[position] if position in self.resfile else ''

    def __setitem__(self, key, value):
        self.__update(pos=key, allowed=value, union=False)

    def write_resfile(self, path, sort=True):
        """Writes the resfile to a new file located in path"""
        if sort:
            self.__sort()
        lines = ['nataa\nstart\n']
        lines += ['{}\t{}\tPIKAA\t{}\n'.format(k[:-1], k[-1], allowed)
                  for k, allowed in self.resfile.items() if allowed]
        open(path, 'w').writelines(lines)

    def write_reduces_resfile(self, path, sort=True):
        if sort:
            self.__sort()
        lines = ['{}{}\t{}\r\n'.format(k[:-1], k[-1], allowed)
                 for k, allowed in self.resfile.items()]
        open(path, 'w').writelines(lines)


    def rearrange_resfile(self, remove_only_wt=False):
        """
        Rearranges the allowed residues to place thw WT residues first.
        :param resfile: a resfile parsed to tuples (from parse_res_file).
        :param pdb_path: path to a PDB file corresponding to the resfile
        :param remove_only_wt: if True, removes records with only the WT
        allowed
        """
        if not self.pdb:
            LGR.warning('pdb file was not set, cannot rearrange resfile')
            return
        resfile_new = OrderedDict()
        for residue, allowed_res in self.resfile.items():
            res = int(residue[:-1])
            chain = residue[-1]
            wt = seq1(self.pdb[0][chain][res].resname)
            if 1 == len(allowed_res) and allowed_res == wt and remove_only_wt:
                continue
            try:
                wt_i = allowed_res.index(wt)
                allowed_res = wt + allowed_res[:wt_i] + allowed_res[wt_i + 1:]
            except ValueError: # if wt is not in the residues
                LGR.warning(('WT residue {} was not in position {} in the res'
                             ' file').format(wt, residue))
                allowed_res = wt + allowed_res
            resfile_new[residue] = allowed_res
        self.resfile  = resfile_new
        self.rearranged = True

    def get_allowed(self, position):
        position = str(position)
        return self.resfile[position] if position in self.resfile else ''

    def get_wt(self, position):
        return self.resfile[position][0]

    def get_positions(self):
        return list(self.resfile.keys())

    def chains(self):
        return list(set([p[-1] for p in self.resfile.keys()]))

    def possible_permutations(self):
        """
        If too big, will return 10**20
        """
        items = [len(v) for k, v in self.resfile.items()]
        permuts = reduce(lambda x, y: x * y, items, 1)
        return permuts if permuts > 0 else 10**20
        #return np.prod([len(v) for k, v in self.resfile.items()])

    def renumber(self, newk):
        """
        BE CAREFUL when using!!
        :param newk: Dictionary with [old key] = new key for example
        query.structure.get_renumbered(resfile.get_positions(), org2new=True)
        """
        self.resfile = OrderedDict([(newk[res], allowed)
                                    for res, allowed in self.resfile.items()])

    def remove_position(self, pos):
        """"""
        if pos in self.get_positions():
            del self.resfile[pos]
        else:
            LGR.warning(f'position {pos} is not in resfile')



def unify_res_file(paths, pdb=None, out_path=os.getcwd()):
    """Unifies resfiles, each threshold into a different file.
    The files shoukd obviously come from the same protein
    :param paths: paths to resfiles. should not have points in the file name
    except from the ones indicating the threshold
    :param pdb: path to a pdb file. will change the first aa in each row to the
    wt
    :param out_path: path to directory where to print the unified resfiles
    :return: a dictionary [threshold] = unified resfile
            and the paths to all the new unified resfilles
    """
    d = dict()
    for f in paths:
        f = os.path.abspath(f)
        threshold = float(f[f.find('.') + 1:])
        rf = ResFile(f)
        if pdb:
            rf.set_pdb(pdb)
            rf.rearrange_resfile()
        if threshold in d.keys():
            d[threshold] += rf
        else:
            d[threshold] = rf
    out_paths = list()
    for threshold, resfile in d.items():
        path = os.path.join(out_path, 'designable_aa_resfile.' + str(threshold))
        out_paths.append(path)
        resfile.write_resfile(path)
    return d, out_paths
