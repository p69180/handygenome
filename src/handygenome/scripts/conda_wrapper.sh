#!/bin/bash
# 210913
	# init
	# set -eu leads to an error while activating conda R environment
# 220715
    # Does not depend on envbash_wrapper

#set -eu


set_constants()
{
	default_condadir=/home/users/pjh/tools/miniconda/210821/miniconda3
	default_envname=genome_v5
	#default_condadir=/opt/miniconda_210914/miniconda3
    env=/usr/bin/env
    bash=/usr/bin/bash
}

usage()
{
	(
		echo "Usage: $(basename $0) [--condadir <conda base direcory>] --envname <envname> arg1 [arg2 [...]]"
		echo "    arg1, arg2, ... are evaluated under the conda environment."
		echo "    --envname : (optional) conda environment name to activate. Default: $default_envname"
        echo "    --condadir : (optional) conda top directory (dirname of 'envs' directory). Default: $default_condadir"
	) > /dev/stderr
	exit 1
}


parse_arguments()
{
	declare -ga args
	posarg_flag=0

	while [[ $# -gt 0 ]] ; do
		if [[ $posarg_flag = 0 ]] ; then
			case "$1" in
				-*)
					case "$1" in 
						--envname)
							shift ; envname="$1" ; shift ;;
						--condadir)
							shift ; condadir="$1" ; shift ;;
						-h|--help)
							usage ;;
						--)
							posarg_flag=1
							shift
							;;
						*)
							echo Invalid argument: $1 > /dev/stderr
							usage
							;;
					esac
					;;
				*)
					posarg_flag=1
					;;
			esac
		elif [[ $posarg_flag = 1 ]] ; then
			while true ; do
				if [[ $# = 0 ]] ; then
					break
				else
					args+=( "$1" )
					shift
				fi
			done
		fi
	done
}


sanity_check()
{
	if [[ ${#args[@]} = 0 ]] ; then
		echo "There must be at least one argument after --args option."
		exit 1
	fi

	#if [[ -z ${envname:-} ]] ; then
	#	echo "--envname option is required." > /dev/stderr
	#	exit 1
	#fi
}


set_default_params()
{
	if [[ -z ${condadir:-} ]] ; then
		condadir=$default_condadir
	fi

	if [[ -z ${envname:-} ]] ; then
		envname=$default_envname
	fi
}


# MAIN
set_constants
parse_arguments "$@"
sanity_check
set_default_params

$env -i HOME=${HOME} $bash --noprofile --norc <(
cat <<EOF
# writing "set -eu" here results in error.

eval $(declare -p args)
eval $(declare -p condadir)
eval $(declare -p envname)

eval "\$(\${condadir}/bin/conda shell.bash hook)"
conda activate \$envname
export LD_LIBRARY_PATH=\${condadir}/envs/\${envname}/lib
eval "\${args[@]}"
EOF
)

