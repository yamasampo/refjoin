# Copyright 2019 by Kent Kawashima.  All rights reserved.

CODON_INTCODE = {cod: i for i, cod in enumerate(a+b+c for a in 'ATCGNX' for b in 'ATCGNX' for c in 'ATCGNX')}
CODON_INTCODE['---'] = 196

INTCODE_CODON = {i: codon for codon, i in CODON_INTCODE.items()}
