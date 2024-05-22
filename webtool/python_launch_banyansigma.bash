#!/bin/bash
#Keywords: NAME, RA, DEC, PMRA, EPMRA, PMDEC, EPMDEC, RV, ERV, PLX, EPLX, IP
export PATH="/usr/local/anaconda3/bin:$PATH"
export HOME="/home/ipm/banyan/banyansigma/answer"
wwwdir='/home/ipm/banyan'
python3 -c "import sys; sys.path.append('$wwwdir'); from banyan_sigma_wrapper import banyan_sigma_wrapper; void = banyan_sigma_wrapper(name='$1',ra='$2',dec='$3',pmra='$4',epmra='$5',pmdec='$6',epmdec='$7',rv='$8',erv='$9',plx='${10}',eplx='${11}',ip='${12}')"
echo "Success running python_launch_banyansigma.bash"
