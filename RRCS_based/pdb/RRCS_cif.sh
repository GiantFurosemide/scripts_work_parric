#!/bin/bash
# usage:RRCS_cif.sh 5kvm.pdb
#
# Please set scripts lib which contains RRSC.py cif2pdb.py to SCRIPTS_DIR
SCRIPTS_DIR=$PWD

RRSC_cif(){
	
	cif_f=$1
	pdb_f="$cif_f".pdb
	python3 $SCRIPTS_DIR/cif2pdb.py $cif_f
	python3 $SCRIPTS_DIR/RRCS.py $pdb_f
}


RRSC_cif $1
