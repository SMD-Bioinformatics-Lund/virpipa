# Resistance numbering reference

`AF011751.1.fasta` contains the 9,599-base H77 sequence named as the amino-acid numbering reference in [the geno2pheno HCV methods](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0155869), section 2.1.2.3. Accession: [AF011751.1](https://www.ncbi.nlm.nih.gov/nuccore/AF011751.1).

The sequence was extracted from the supplied two-reference FASTA on 2026-10-06 and checked against the supplied MAFFT alignment. SHA-256 of the uppercase sequence without header or whitespace: `72fb77d5ba74c2d865e542b48e27f8662d3173752dbf217a6e5f08dc1537aa73`.

NS3 occupies nucleotides 3420–5312, NS5A 6258–7601, and NS5B 7602–9374 (1-based inclusive). These coordinates match the existing AF009606.1 reference; the 47-base length difference occurs after the coding region. In these three genes, the only amino-acid difference is NS3 E176G (AF009606.1 → AF011751.1), outside the positions scored by the committed rules snapshot. NS5A nucleotide differences are synonymous and NS5B is identical.

The pipeline, module test, and standalone resistance callers use this file. Its subdirectory keeps it outside the top-level assembly-reference glob. Resistance JSON records the reference identifier and file checksum.

## Verified subtype alignment bridges

`panel.json` defines one reporting genome for each of the ten additional subtypes observed in the supplied geno2pheno results, plus AF011751.1 for 1a and AJ851228.1 for the genotype-1 case with unassigned subtype. AJ851228.1 is the reference observed for that case, not a universal reference for every unassigned genotype-1 sequence; its subtype remains null. The caller uses exact assigned subtype matches as alignment bridges, without changing pipeline subtype assignment. Unsupported and unassigned subtypes use an explicit conservative H77 fallback.

The eleven added genomes in `genomes/` are original accession-version sequences retrieved from NCBI GenBank on 2026-10-06. Each was checked against its installed HCV-GLUE BLAST entry: in these eleven entries, removing the database's `N` padding exactly recovers the original sequence. This is a verified property of these records, not a general instruction to delete Ns from genomes.

Six common-subtype repository genomes have also been verified exactly against that alignment: 1b D90208, 2a D00944, 2b D10988, 3a D17763, 5a Y13184 and 6a Y12083. Their original source files and sequence checksums are recorded in the manifest; accession versions and equivalence to geno2pheno reporting-reference choice are not established. These copies are kept inside the panel directory so Nextflow can stage the complete, self-contained reference input. Repository assembly references lacking an exact verified aligned counterpart were not promoted.

`alignments/NS3.fasta`, `NS5A.fasta`, and `NS5B.fasta` preserve the installed HCV-GLUE alignment columns, representing verified padding as gaps. The H77 member uses AF011751.1's coding sequence in the same coordinate columns. Gene boundaries are projected through that alignment and checked against each GenBank polyprotein CDS reading frame. The manifest pins genome sequence checksums, gene coordinates, completeness, and the source database checksum. DQ155561.1 has an incomplete NS5B end, which is explicitly marked.

These files supply reference sequences and coordinate links; they do not define a rare-to-common subtype translation table or a susceptible* scoring baseline. Broadening bridge coverage and validating per-region classification remain future work.
