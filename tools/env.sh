# source this file: puts every module on PYTHONPATH so that `python -m landmarks ...` etc. work from anywhere
#   source tools/env.sh
_here="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"
_p="$_here/common"
for d in "$_here"/[0-9][0-9]_*; do [ -d "$d" ] && _p="$_p:$d"; done
export CICHLID_REPO="$_here"
export PYTHONPATH="$_p${PYTHONPATH:+:$PYTHONPATH}"
unset _here _p
