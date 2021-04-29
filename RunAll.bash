#!/bin/bash
python 1_ProteinCellCycleClusters.py &> output/1_ProteinCellCycleClusters.log
python 2_ProteinFucciPsuedotime.py &> output/2_ProteinFucciPsuedotime.log
python 3_RNAFucciPseudotime.py &> output/3_RNAFucciPseudotime.log
python 4_TemporalDelay.py &> output/4_TemporalDelay.log
python 5_ProteinProperties.py &> output/5_ProteinProperties.log