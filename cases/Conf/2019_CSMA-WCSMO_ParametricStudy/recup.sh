#!/bin/bash

classicPara='classic/parametricStudy'
xfemPara='xfem/parametricStudy'
set -x

rsync -avu Gypse:/home/Work/SILEX3.4/calculs/cavity_and_plate/$classicPara/*.mat $classicPara/RAW/.
rsync -avu Nancite:~/Documents/Recherche/OptimAcoustic/SILEXfull/calculs/cavity_and_plate/$classicPara/*.mat $classicPara/RAW/.

rsync -avu Gypse:/home/Work/SILEX3.4/calculs/cavity_and_plate/$classicPara/*.log $classicPara/.
rsync -avu Nancite:~/Documents/Recherche/OptimAcoustic/SILEXfull/calculs/cavity_and_plate/$classicPara/*.log $classicPara/RAW/.

rsync -avu Gypse:/home/Work/SILEX3.4/calculs/cavity_and_plate/$xfemPara/*.mat $xfemPara/.
rsync -avu Nancite:~/Documents/Recherche/OptimAcoustic/SILEXfull/calculs/cavity_and_plate/$xfemPara/*.mat $xfemPara/RAW/.

rsync -avu Gypse:/home/Work/SILEX3.4/calculs/cavity_and_plate/$xfemPara/*.log $xfemPara/.
rsync -avu Nancite:~/Documents/Recherche/OptimAcoustic/SILEXfull/calculs/cavity_and_plate/$xfemPara/*.log $xfemPara/RAW/.

