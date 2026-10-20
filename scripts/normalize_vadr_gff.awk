# Preserve VADR qualifiers and identifiers when projecting peptides to genes.
BEGIN { FS=OFS="\t"; split("core E1 E2 p7 NS2 NS3 NS4A NS4B NS5A NS5B", genes, " "); print "##gff-version 3" }
/^#/ { next }
NF == 9 {
    sub(/_pilon$/, "", $1)
    if ($3 == "5'UTR" || $3 == "3'UTR") { print; next }
    if ($3 != "mat_peptide" && $3 != "misc_feature") next
    for (i=1; i<=10; i++) {
        gene=genes[i]
        if ($9 ~ "[= ]" gene "[ ;]") {
            if (seen[gene]++) { print "Duplicate VADR feature for " gene > "/dev/stderr"; exit 1 }
            original=$3; $3="gene"
            $9=$9 "gene=" gene ";annotation_type=" original ";"
            if (original == "misc_feature") $9=$9 "annotation_inferred=true;"
            print
        }
    }
}
