import socket

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
breakends = importlib.import_module('.'.join([top_package_name, 'sv', 'breakends']))


class IGVHandle:
    def __init__(self, port):
        self.port = port

    def cmd(self, msg):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.connect(('localhost', self.port))
            s.send(msg.encode('utf-8'))
            s.shutdown(socket.SHUT_WR)
            data = s.recv(1024).decode().replace('\n', '')
            print(data)
        
    def goto(self, loci, width=500):
        assert isinstance(loci, (list, tuple)), (
            f'"loci" argument must be a list or a tuple.')

        cmd_src = list()
        for locus in loci:
            if isinstance(locus, str):
                cmd_src.append(locus)
            elif isinstance(locus, common.Vcfspec):
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
                raise Exception(f'Entries of "loci" must be a string, '
                                f'Vcfspec, or Breakends.')

        self.cmd(f'goto {" ".join(cmd_src)}')

    def load(self, filepaths):
        assert isinstance(filepaths, (list, tuple)), (
            f'"filepaths" argument must be a list or a tuple.')

        for filepath in filepaths:
            self.cmd(f'load {filepath}')
    
    def zoomin(self, n=1):
        for _ in range(n):
            self.cmd('zoomin')

    def zoomout(self, n=1):
        for _ in range(n):
            self.cmd('zoomout')

    def viewaspairs(self, turnoff=False):
        self.cmd('viewaspairs')

    def viewaspairs_off(self):
        self.cmd('viewaspairs false')
        
