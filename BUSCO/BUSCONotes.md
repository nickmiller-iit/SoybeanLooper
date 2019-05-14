# BUSCO analysis

We are planning a USDA-NIFA foundational proposal. As part of this, we will improve the genome assembly. It would be good to know how good it is already. To that end, we will get some BUSCO stats.

Ran BUSCO in genome mode, with the endopterygote lineage data. The initial Augustus run was in "fly" mode because there are no trained params for lepidoptera. Once the config.ini file was set up correctly for the conda install of BUSCO, things went pretty smoothly.

## Results

	C:63.6%[S:63.1%,D:0.5%],F:27.5%,M:8.9%,n:2442

	1551	Complete BUSCOs (C)
	1540	Complete and single-copy BUSCOs (S)
	11	Complete and duplicated BUSCOs (D)
	671	Fragmented BUSCOs (F)
	220	Missing BUSCOs (M)
	2442	Total BUSCO groups searched

All told, this isn't bad. Granted it's not what we would want to see for a draft genome assembly, but there's an impressive number of complete, or at least fragmented BUSCOs for a quick and dirty assembly.
