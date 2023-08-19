#!/home/users/pjh/tools/miniconda/221104/miniconda3/envs/genome_v7/bin/python
import sys
import os
import stat
import re
import textwrap
import shutil

PYTHON = sys.executable
PROJECT_PATH = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
PACKAGE_PATH = os.path.join(PROJECT_PATH, 'handygenome')
TARGETDIR = os.path.join(PROJECT_PATH, 'bin')


def main():
    # delete existing ones
    if os.path.exists(TARGETDIR):
        shutil.rmtree(TARGETDIR)
    os.mkdir(TARGETDIR)

    # create new ones
    for top, dirs, files in os.walk(os.path.join(PACKAGE_PATH, 'bin')):
        for fname in files:
            if fname.endswith('.py') and not fname.startswith('_'):
                target_name = os.path.join(TARGETDIR, re.sub('.py$', '', fname))

                abspath = os.path.join(top, fname)
                nopy = re.sub('.py$', '', abspath)
                strip_leading = re.sub(PROJECT_PATH + '/', '', nopy)
                modulepath = '.'.join(strip_leading.split('/'))

                text = textwrap.dedent(
                    f"""\
                    #!{PYTHON}
                    import sys
                    sys.path.append('{os.path.dirname(PACKAGE_PATH)}')
                    from {modulepath} import main
                    main(sys.argv[1:])"""
                )
                with open(target_name, 'wt') as outfile:
                    outfile.write(text)

                os.chmod(
                    target_name, 
                    (
                        stat.S_IRUSR 
                        | stat.S_IWUSR
                        | stat.S_IXUSR
                        | stat.S_IRGRP
                        | stat.S_IXGRP
                    )
                )


if __name__ == '__main__':
    main()
