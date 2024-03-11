import re
import os
import socket

import numpy as np

import handygenome.sv.breakends as breakends
#import handygenome.variant.vcfspec as libvcfspec
from handygenome.variant.vcfspec import Vcfspec


class IGVHandle:
    def __init__(self, port, pathsub=None):
        self.port = port
        self.pathsub = pathsub

    def cmd(self, msg):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.connect(('localhost', self.port))
            s.send(msg.encode('utf-8'))
            s.shutdown(socket.SHUT_WR)
            data = s.recv(1024).decode()
            print(data, end='')

    def new(self):
        self.cmd('new')
        
    def goto(self, *args, width=200):
        cmd_src = list()
        for locus in args:
            if isinstance(locus, tuple):
                if len(locus) == 3:
                    cmd_src.append(f'{locus[0]}:{locus[1] - width}-{locus[2] + width}')
                elif len(locus) == 2:
                    cmd_src.append(f'{locus[0]}:{locus[1] - width}-{locus[1] + width}')
                else:
                    raise Exception(f'Argument of type "tuple" must be of length 2 or 3.')
            elif isinstance(locus, str):
                locus_sp = locus.split(':')
                chrom = locus_sp[0]
                start1, end1 = locus_sp[1].split('-')
                start1 = int(start1)
                end1 = int(end1)
                cmd_src.append(f'{chrom}:{start1 - width}-{end1 + width}')
            elif isinstance(locus, Vcfspec):
                chrom = locus.chrom
                if locus.get_mutation_type(alt_index=0) == 'del':
                    start1 = locus.pos + 1 - width
                    end1 = locus.pos + 1 + width
                else:
                    start1 = locus.pos - width
                    end1 = locus.pos + width

                cmd_src.append(f'{chrom}:{start1}-{end1}')
            elif isinstance(locus, breakends.Breakends):
                chrom1 = locus.chrom_bnd1
                start1 = locus.pos_bnd1 - width
                end1 = locus.pos_bnd1 + width
                cmd_src.append(f'{chrom1}:{start1}-{end1}')

                chrom2 = locus.chrom_bnd2
                start2 = locus.pos_bnd2 - width
                end2 = locus.pos_bnd2 + width
                cmd_src.append(f'{chrom2}:{start2}-{end2}')
            else:
                raise Exception(
                    f'Arguments must be a str, tuple, Vcfspec, or Breakends.'
                )

        self.cmd(f'goto {" ".join(cmd_src)}')

    def load(self, filepaths):
        # make into absolute paths
        filepaths = [os.path.abspath(x) for x in np.atleast_1d(filepaths)]
        # do substition
        if self.pathsub is not None:
            filepaths = [re.sub(self.pathsub[0], self.pathsub[1], x) for x in filepaths]

        for path in filepaths:
            self.cmd(f'load {path}')
    
    def zoomin(self, n=1):
        for _ in range(n):
            self.cmd('zoomin')

    def zoomout(self, n=1):
        for _ in range(n):
            self.cmd('zoomout')

    def viewaspairs(self):
        self.cmd('viewaspairs')

    def viewaspairs_off(self):
        self.cmd('viewaspairs false')

    def snapshot(self, path, maxPanelHeight=1000):
        basename = os.path.basename(path)
        dirname = os.path.dirname(path)

        if not os.path.exists(dirname):
            os.makedirs(dirname)
        self.cmd(f'snapshotDirectory {dirname}')

        self.cmd(f'maxPanelHeight {maxPanelHeight}')
        self.cmd(f'snapshot {basename}')

    def group(self, option, tagname=None):
        if tagname is not None:
            assert option == 'TAG'
            cmd = f'group {option} {tagname}'
        else:
            cmd = f'group {option}'

        self.cmd(cmd)

    def group_pairorient(self):
        self.group('PAIR_ORIENTATION')
        
    def sort(self, key):
        self.cmd(f'sort {key}')
