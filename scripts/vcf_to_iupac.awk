#!/usr/bin/awk -f
# vcf_to_iupac.awk
# Usage:
#   awk -v MIN_AF=0.10 -v MIN_DP=5 -f vcf_to_iupac.awk sample.vcf polished.fa > polished_iupac.fa

BEGIN {
    FS = "\t"; OFS = "\t";

    # Defaults if user didn't provide -v MIN_AF / -v MIN_DP
    if (MIN_AF == "") MIN_AF = 0.10;
    if (MIN_DP == "") MIN_DP = 5;

    # IUPAC table for REF,ALT (biallelic top-alt only)
    iupac["A,G"] = "R"; iupac["G,A"] = "R";
    iupac["C,T"] = "Y"; iupac["T,C"] = "Y";
    iupac["G,C"] = "S"; iupac["C,G"] = "S";
    iupac["A,T"] = "W"; iupac["T,A"] = "W";
    iupac["G,T"] = "K"; iupac["T,G"] = "K";
    iupac["A,C"] = "M"; iupac["C,A"] = "M";

    # Which ARGV index is the VCF? We expect user to supply vcf first, fasta second
    vcf_index = 1;
    fasta_index = 2;
}

# -----------------------
# Phase 1: read VCF (first file)
FILENAME == ARGV[vcf_index] {
    # skip meta headers
    if ($0 ~ /^##/) next;

    # header line with column names
    if ($0 ~ /^#CHROM/) {
        # sample columns start at column 10 (VCF standard)
        sample_col = 10;
        next;
    }

    # Parse VCF record
    chrom = $1;
    pos = $2 + 0;
    ref = $4;
    split($5, alts, ","); n_alt = length(alts);

    # Parse FORMAT and sample fields
    split($9, fmtarr, ":");
    # sample field (single-sample VCF assumed) is column sample_col
    sample_field = $(sample_col);

    # If sample field is missing or '.', skip
    if (sample_field == "" || sample_field == ".") next;
    split(sample_field, valarr, ":");

    # find AD and DP indices (in FORMAT)
    ad_idx = -1; dp_idx = -1;
    for (i = 1; i <= length(fmtarr); i++) {
        if (fmtarr[i] == "AD") ad_idx = i;
        if (fmtarr[i] == "DP") dp_idx = i;
    }

    # Get DP (prefer FORMAT/DP; fallback to INFO DP)
    dp = -1;
    if (dp_idx > 0 && valarr[dp_idx] != ".") {
        dp = valarr[dp_idx] + 0;
    } else {
        # look in INFO column ($8) for DP=number
        if (match($8, /(^|;)DP=([0-9]+)/, m)) {
            dp = m[2] + 0;
        }
    }
    if (dp < MIN_DP) next;   # not enough coverage

    # Get AD (allelic depths) - required
    if (ad_idx < 0 || valarr[ad_idx] == ".") next;
    split(valarr[ad_idx], adarr, ",");

    # sanity: adarr[1] is ref count; adarr[2..] are alt counts
    ref_count = (adarr[1] + 0);
    total = ref_count;
    max_alt_count = 0;
    max_alt_idx = 0;
    for (i = 1; i <= n_alt; i++) {
        alt_count = (adarr[i+1] + 0);
        alt_counts[i] = alt_count;
        total += alt_count;
        if (alt_count > max_alt_count) { max_alt_count = alt_count; max_alt_idx = i; }
    }
    if (total == 0) next;

    # Compute frequency for top ALT (relative to total depth)
    af_top_alt = max_alt_count / total;

    if (af_top_alt >= MIN_AF && max_alt_idx > 0) {
        altbase = alts[max_alt_idx];
        key = chrom SUBSEP pos;              # per-contig key
        iupackey = ref "," altbase;
        if (iupac[iupackey] != "") {
            replacement[key] = iupac[iupackey];
        }
        # if no IUPAC mapping (non-ACGT), we skip replacement
    }
    next;
}

# -----------------------
# Phase 2: read FASTA (second file) and apply replacements
FILENAME == ARGV[fasta_index] {
    # header line
    if ($0 ~ /^>/) {
        # extract contig name (first token after '>')
        header = $0;
        split(header, parts, /[ \t]/);
        curchr = substr(parts[1], 2);  # remove '>'
        print $0;                      # print header unchanged
        pos = 0;                       # reset coordinate for this contig
        next;
    }

    # sequence line: iterate bases and replace if requested
    seq = $0;
    out = "";
    for (i = 1; i <= length(seq); i++) {
        pos++;
        base = substr(seq, i, 1);
        key = curchr SUBSEP pos;
        if (key in replacement) {
            out = out replacement[key];
        } else {
            out = out base;
        }
    }
    print out;
    next;
}

# If there are extra files, do nothing
