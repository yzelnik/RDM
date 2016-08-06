#!/bin/tcsh
qstat -u kinast | grep qw | gawk '{print $1}' > Worklist
# qstat -u kinast | gawk '{print $1}' > Worklist
cat Worklist | xargs qdel
qstat -u kinast | grep r | gawk '{print $1}' > Worklist
cat Worklist | xargs qdel
