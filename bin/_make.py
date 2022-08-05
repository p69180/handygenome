#!/home/users/pjh/tools/miniconda/210821/miniconda3/envs/genome_v5/bin/python
import os
import stat
import re
import textwrap

PYTHON = '/home/users/pjh/tools/miniconda/210821/miniconda3/envs/genome_v5/bin/python'
PKGDIR = '/home/users/pjh/scripts/python_genome_packages'
PKGPATH = PKGDIR + '/handygenome'
TARGETDIR = '.'


# delete existing ones
for fname in os.listdir('.'):
    if not fname.startswith('_'):
        os.remove(fname)

# create new ones
for top, dirs, files in os.walk(os.path.join(PKGPATH, 'tools')):
    for fname in files:
        if fname.endswith('.py') and not fname.startswith('_'):
            target_name = os.path.join(TARGETDIR, re.sub('.py$', '', fname))

            abspath = os.path.join(top, fname)
            nopy = re.sub('.py$', '', abspath)
            strip_leading = re.sub(PKGDIR + '/', '', nopy)
            modulepath = '.'.join(strip_leading.split('/'))

            text = textwrap.dedent(f"""\
                #!{PYTHON}
                import sys
                sys.path.append('{os.path.dirname(PKGPATH)}')
                from {modulepath} import main
                main(sys.argv[1:])""")
            with open(target_name, 'wt') as outfile:
                outfile.write(text)

            os.chmod(target_name, stat.S_IRUSR 
                                  | stat.S_IWUSR
                                  | stat.S_IXUSR
                                  | stat.S_IRGRP
                                  | stat.S_IXGRP)
