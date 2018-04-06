#!/bin/csh -f

set nodes = `cat $PBS_NODEFILE`
set old_node = ''
@ n = 1
foreach node ( $nodes )
  if ($node != $old_node) then
    ssh -q $node killall $1 &
    @ n = $n + 1
    set old_node = $node
    echo $n $node
  endif
end
wait
echo Finished

exit 0
