#!/usr/bin/awk -f
# Coordinate-preserving consensus from a single-sample, all-sites VCF.
# Analytical modes mask missing observations; majority retains scaffold bases.
BEGIN {
    FS = "\t"; OFS = "\t"
    if (MIN_AF == "") MIN_AF = 0.10
    if (MIN_DP == "") MIN_DP = 5
    if (MODE == "") MODE = (MIN_AF == 1 ? "majority" : "iupac")
    if (MODE != "iupac" && MODE != "majority" && MODE != "analytical_majority") {
        print "Invalid consensus MODE" > "/dev/stderr"; exit 1
    }
    codes[1]="A"; codes[2]="C"; codes[4]="G"; codes[8]="T"
    codes[3]="M"; codes[5]="R"; codes[9]="W"; codes[6]="S"
    codes[10]="Y"; codes[12]="K"; codes[7]="V"; codes[11]="H"
    codes[13]="D"; codes[14]="B"; codes[15]="N"
    bits["A"]=1; bits["C"]=2; bits["G"]=4; bits["T"]=8
}

FILENAME == ARGV[1] {
    if ($0 ~ /^#/) next
    key = $1 SUBSEP $2
    vcf_contigs[$1]=1
    reason[key]="invalid_evidence"
    ref = toupper($4)
    nalt = split(toupper($5), alt, ",")
    if ($5 == ".") nalt=0
    nf = split($9, fmt, ":"); split($10, values, ":")
    ad = ""; depth = -1
    for (j=1; j<=nf; j++) {
        if (fmt[j] == "AD") ad = values[j]
        if (fmt[j] == "DP" && values[j] ~ /^[0-9]+$/) depth = values[j]+0
    }
    if (ad == "" || ad == ".") next
    nad = split(ad, counts, ",")
    if (nad != nalt+1) next
    delete totals
    valid = 1; total = 0; indel = 0
    for (j=1; j<=nad; j++) {
        if (counts[j] !~ /^[0-9]+$/) { valid=0; break }
        base = (j == 1 ? ref : alt[j-1])
        if (base in bits) { totals[base] += counts[j]; total += counts[j] }
        else if (base != "." && base != "<*>" && base != "<NON_REF>") indel=1
    }
    if (indel) {
        indel_records++
        for (j=0; j<length(ref); j++) {
            blocked[$1 SUBSEP ($2+j)]=1
            reason[$1 SUBSEP ($2+j)]="indel_ambiguity"
        }
    }
    if (valid && !indel) usable_depth[key]=total
    # Mixed SNV/indel records cannot establish a nucleotide frequency safely.
    if (!valid || indel || total <= 0 || total < MIN_DP || (depth >= 0 && depth < MIN_DP)) next
    if (key in replacement) { blocked[key]=1; reason[key]="duplicate_record"; next }
    mask=0; maximum=-1; chosen=""; tied=0
    for (j=1; j<=4; j++) {
        base=substr("ACGT", j, 1)
        if (totals[base] > 0 && totals[base]/total >= MIN_AF) mask += bits[base]
        if (totals[base] > maximum) { maximum=totals[base]; chosen=base; tied=1 }
        else if (totals[base] == maximum) tied++
    }
    if (tied > 1) ties++
    reason[key]=""
    if (MODE == "majority" || MODE == "analytical_majority") {
        if (ref in bits && totals[ref] == maximum) chosen=ref
        replacement[key]=chosen
    } else {
        replacement[key]=(mask in codes ? codes[mask] : "N")
        if (mask == 15) reason[key]="all_four_alleles"
    }
    next
}

FILENAME == ARGV[2] {
    if ($0 ~ /^>/) {
        nh=split(substr($0,2), header, /[ \t]/); chrom=header[1]; pos=0
        for (vcf_contig in vcf_contigs) if (vcf_contig != chrom) {
            print "VCF contig does not match consensus reference" > "/dev/stderr"; exit 1
        }
        print ">" chrom
        next
    }
    output=""
    for (j=1; j<=length($0); j++) {
        pos++; key=chrom SUBSEP pos
        base=(key in replacement && !(key in blocked) ? replacement[key] : (MODE == "majority" ? toupper(substr($0,j,1)) : "N"))
        output=output base
        if (EVIDENCE_FILE != "") {
            why=(key in blocked ? reason[key] : (key in replacement ? reason[key] : (reason[key] == "invalid_evidence" && !(key in usable_depth) ? "invalid_evidence" : (usable_depth[key] > 0 ? "insufficient_or_invalid_evidence" : "missing_coverage"))))
            print pos, (usable_depth[key]+0), why > EVIDENCE_FILE
        }
    }
    print output
}

END {
    if (STATS_FILE != "") print (indel_records+0), (ties+0) > STATS_FILE
}
