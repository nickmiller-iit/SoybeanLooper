# Acetylcholinesterase gene

Getting the ryanodine receptor gene proved tricky/messy. Acetylcholinesterase is a much smaller protein: ~640 - 690 aa. Lest see if we get better results here.

Used the same basic approach as for the ryanodine receptor. Used tblastn with three noctuid achE full aa sequences as a qury to get candidate scaffolds, then used exonerate to align the *H. armigera* aa sequence to candidate scaffolds.

In this case it was pretty easy to fish out the scaffolds. There were a bunch of hits with relativly low % identity and 4 hits with % identity > 90.

scaffold43722: 1 - 400
scaffold1417: 395 - 501
scaffold 23910: 494 - 549
scaffold: 54604: 550 - 638

The *H. armigera* protein is 638 aas, so it looks like we managed to recover the entire gene, albeit in a fragmented state.
