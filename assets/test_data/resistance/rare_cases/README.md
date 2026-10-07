# Rare-subtype comparison cases

Eleven supplied 15% IUPAC genomes cover 2c, 2i, 2k, 3h, 4c, 4d, 4k, 4m, 4r, 6e, and genotype 1 with no assigned subtype. Headers were replaced with randomly assigned `G2P_RARE_###` identifiers. Original identifiers and source files must remain outside Git.

`expected.json` preserves the matched geno2pheno CSV outcomes for all 16 drugs, including scored amino acids, raw fold/evidence fields, per-gene subtype labels and selected references, coverage, and subtype-relative substitutions. The genotype-1 case retains raw service label `1*2` while storing genotype `1` and an unassigned subtype. Sequence checksums ensure that each expectation stays associated with the submitted genome.

These are target results for full caller parity. The integrity tests validate the inputs and expectations; passing them does not establish that the current caller reproduces the online predictions or scored-mutation lists. Raw evidence codes and fold values are retained without adding a clinical interpretation.

## Recovery and remaining validation

The temporary re-identification key and source copies were lost when WSL restarted on 2026-10-06. On 2026-10-07, the restored `rare_subtypes.fa` and `g2p_hcv_data.csv` allowed all eleven associations to be recovered by sequence checksum and checked against regional subtype/reference results. The recovered key is stored in a persistent, private directory outside Git.

For recovery, match sequences by SHA-256 of uppercase bases without whitespace, using `sequence_sha256` in `expected.json`; do not infer identifiers from fixture order. The supplied CSV identifiers have a `(No identifier given) - ` prefix before the FASTA identifier. Store any reconstructed key outside the repository, not in `/tmp`.

The replacement local VADR run completed on 2026-10-07 using the existing pipeline wrapper, processing one genome at a time. The `.gff3` fixtures retain its three resistance-gene annotations without correcting the input sequences or boundaries. They are reproducible caller inputs, not a declaration that each whole genome is biologically intact. The candidate reference panel is not active in the production workflow.

The user reports mixed signals in `G2P_RARE_001` (2i), `G2P_RARE_005` (3h), `G2P_RARE_007` (2c), `G2P_RARE_009` (6e), and `G2P_RARE_010` (4r). These remain valid comparisons of the submitted IUPAC sequences, not ground truth for strain-specific haplotypes or linked compound mutations. The other cases have not been independently certified single-strain.

Boundary review found that VADR trims an ambiguous final base from 2i NS5A (7669 rather than reference-supported 7670) and 4c NS3 (5313 rather than 5314). The 3h NS5A interval, 6269–7622, matches subtype-reference alignment but contains internal indels and downstream stop codons. Read review supports two bases missing from that consensus near position 7315, distinct from an approximately equal three-base length mixture near the gene end. All 35 catalogue resistance-site codons across these three genes independently matched nucleotide projections through their subtype references. No boundaries or bases have been silently repaired in these fixtures.

The regression check requires all 16 drug rows and compares resistant/reduced-susceptibility, susceptible, and licensing outcomes. It deliberately collapses susceptible* into susceptible for that check: scored-substitution presentation still needs its own subtype-dependent validation.
