# Identifying MIPs probes with MIPGEN

The MIPGEN tool from the Shendure lab is designed to identify MIPs probes for specified targets. It also checks that probes won't hybridize elsewhere in the genome and calculates a score reflecting the probability that the probe will produce a viable assay. All of this is very helpful. However MIPGEN is intended for designing tiled sets of probes that cover larger features. We don't want to do this, because we are interested in surveying sites scattered throughout the genome. This isn't a big deal as we can just pick one probe from a tiled set.

## MIPGEN is a little fussy

After playing around with MIPgen, I noticed some things that cause it to fail. This was figured out by trial and error, since the MIPGEN documentation and paper didn't really give me any clues as to what was going on. Known issues include:

 * MIPGen chokes on fasta files that contain pipe characters in the identifier lines. This appears to be because it passes identifiers directly to samtools faidx without quoting the identifier.
 * MIPGen fails on features at the start of a contig. No idea why, but avoiding features < 300 bases from the start seems to serve as a workaround.
 * MIPGen fails at the tiling stage on large features. This seems to be because if it cannot tile across the feature it fails.

To work around these issues the basic approach is

 1. Update the redundans assembly to remove pipe chars from identifiers.
 2. Regenerate BED files to match identifiers.
 3. Generate a new BED file with 1bp targets, evenly spaced through the non-repetitive regions at 500bp intervals. targets within 200 bases of the start or end of the scaffold are dropped. This has the effect that, at the tiling stage, MIPgen only needs to "tile" over a single base, so it just picks the best probe for each 1bp target. This also reduces the sizes of the output files considerably.

These workarounds were figured out in discussion with Evan Boyle and Jay Shendure. Be sure to acknowledge them in the paper.

Addistional conversations with Evan established that the short scaffolds we have were causing a problem. Here's Evan's explanation in full:

*"I was able to recognize what is happening. I didn’t anticipate MIPgen being used on such small scaffolds. By design it acquires sequence from +/- 1kb of the area of interest. I found that this information slightly improved performance.

Right not that sequence doesn’t exist because the scaffold is too small. Fortunately, with the default logistic scoring, this sequence is not needed anyway and you can basically turn this off.

I don’t think I have time to add an option to the code (mostly because I am always afraid of breaking it incidentally) in the next few days, but if you are savvy enough, if you go to lines 1126 and 1127 of the main script (mipgen.cpp) and replace the lines as such (simply deleting the extra 1000 nucleotides):"*
 
```
+ boost::lexical_cast<string>(feature->start_position_flanked - max_capture_size) + "-“ \ + boost::lexical_cast<string>(feature->stop_position_flanked + max_capture_size
 + 14) +
 "
 " \
 ```
 
 *"Then as long as you keep -score_method logistic, I think it should run as you expect and tile those regions"*


## Running MIPgen on a reduced set of targets

After attempting to run MIPgen with all possible targets, it became apparent that doing so would
 
  1. Take several days, possibly more than a week.
  2. Generate huge output files, totalling sevaeral Tb

The solution is to sample from the total set of potential targets. Generating a sample of 10000 targets should still give us plenty of scope to select 50 high-quality MIPs that a re scattered throughout the genome. The R script that does the subsampling accepts a random seed as an argument. This is used to ensure that the same sample of targets can be selected every time the script is run.

Made a subsample of 10000 targets and ran MIPgen. The picked_mips file contains 3224 MIPs with logistic scores >= 0.98 (ie very high scoring), in 3097 different scaffolds. 

## Choosing MIPs to test

Details are documented in the R markdown notebook. Basic approach was

 1. Filter picked MIPs to retain the MIPs with the highest logistic score in each scaffold.
 2. Sort the set of MIPs by logistic score
 3. Select the 50 MIPs with the highest logistic scores
