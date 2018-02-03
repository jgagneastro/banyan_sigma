#!/bin/bash
export PATH="/usr/local/anaconda3/bin:$PATH"
export HOME="/home/ipm/banyan/banyansigma/answer"
wwwdir='/home/ipm/banyan'
python3 -c "import sys; sys.path.append('$wwwdir'); from name_resolver_webtool import name_resolver_webtool; void = name_resolver_webtool(name='$1')"
echo "Success running python_resolve_name.bash"
