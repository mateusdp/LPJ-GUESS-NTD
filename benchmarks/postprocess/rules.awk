function o(x) {
  print $1, $2, x
}

BEGIN {
  pft = 21 # why is this thing here?
  lat = 2
  trees = 20
  grass = 19
  tr = 17
  te = 18
  bt = 16
  total = 15
  bne = 3
  bine = 4
  bns = 5
  tene = 6
  tebs = 7
  ibs = 8
  tebe = 9
  trbe = 10
  tribe = 11
  trbr = 12
  bnee = 21
  trbee = 22
}

{
if ($trees > 2.5 && $trbe > .6*$trees) o(9)
else if ($trees > 2.5 && $trbr > .6*$trees) o(10)
else if ($trees > 2.5 && $tr > .5*$trees && (($trbe > $tebe && $trbe > $tebs) || ($trbr > $tebe && $trbr > $tebs))) o(8)
else if ($trees > 2.5 && $bt > .8*$trees && ($bnee > $bns || $ibs > $bns)) o(2)
else if ($trees > 2.5 && $bt > .8*$trees && ($bns > $bnee && $bns > $ibs)) o(1)
else if ($trees > 2.5 && $te > .8*$tress && $tebe > .5*$trees) o(6)
else if ($trees > 2.5 && $te > .8*$tress && $tebs > .5*$trees) o(5)
else if ($trees > 2.5 && $te > .8*$tress && $tene > .5*$trees) o(4)
else if ($trees > 2.5 && $bt > .2*$trees) o(3)
else if ($trees > 2.5) o(7)
else if ($trees > .5 && $trees < 2.5 && $bt > .8*$trees && ($bnee > $bns || $ibs > $bns)) o(2)
else if ($trees > .5 && $trees < 2.5 && $bt > .8*$trees && ($bns > $bnee && $bns > $ibs)) o(1)
else if ($trees > .5 && $trees < 2.5 && $trees > .8*$total) o(15)
else if ($trees > .5 && $trees < 2.5 && $total > 2.5) o(11)
else if ($trees > .5 && $trees < 2.5) o(12)
else if ($trees < .5 && $grass > .2 && $lat > 54) o(18)
else if ($grass > 2.0) o(13)
else if ($trees > .2 && $grass < 1.0) o(16)
else if ($grass > .2) o(14)
else if ($total > .2) o(16)
else if ($total <= .2) o(17)
}
