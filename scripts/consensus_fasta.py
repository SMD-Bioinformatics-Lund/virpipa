#!/usr/bin/env python3
import sys
from collections import Counter

# IUPAC ambiguity codes
IUPAC_CODES = {
    frozenset('A'): 'A',
    frozenset('C'): 'C',
    frozenset('G'): 'G',
    frozenset('T'): 'T',
    frozenset('AG'): 'R',
    frozenset('CT'): 'Y',
    frozenset('GC'): 'S',
    frozenset('AT'): 'W',
    frozenset('GT'): 'K',
    frozenset('AC'): 'M',
    frozenset('CGT'): 'B',
    frozenset('AGT'): 'D',
    frozenset('ACT'): 'H',
    frozenset('ACG'): 'V',
    frozenset('ACGT'): 'N'
}

def read_fasta(file_path):
    """
    Read a FASTA file and return a list of sequences.
    Returns reference sequence separately from other sequences.
    """
    sequences = []
    current_sequence = []
    reference = None
    reference_header = None

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
#            print(line)
#            print("end ofline")
            if line.startswith('>'):
#                print("startswith")
#                print(current_sequence)
                if current_sequence == []:
                    reference_header = line[1:]
#                    print("HEADERTIME")
#                    print(reference_header)
                if current_sequence:
                    if reference is None:
                        reference = ''.join(current_sequence)
                    else:
                        sequences.append(''.join(current_sequence))
                    current_sequence = []
            else:
                current_sequence.append(line.upper())  # Convert to uppercase

    # Handle last sequence
    if current_sequence:
        if reference is None:
            reference = ''.join(current_sequence)
        else:
            sequences.append(''.join(current_sequence))

    return reference, reference_header, sequences

def get_consensus_base(base_counts):
    """
    Determine the consensus base using IUPAC ambiguity codes when appropriate.
    Returns the consensus base and its score.
    """
    total_bases = sum(base_counts.values())
    if not total_bases:
        return '-', 0.0

    # Get the highest count
    max_count = max(base_counts.values())

    # Find all bases that appear with the highest frequency
    most_common_bases = sorted([base for base, count in base_counts.items()
                              if count == max_count])

    # always return a base, even if ambiguous
    return most_common_bases[0], (max_count / total_bases * 100)
    # # If only one base has the highest count, return that base
    # if len(most_common_bases) == 1:
    # return most_common_bases[0], (max_count / total_bases * 100)

    # # Otherwise, use IUPAC code for the ambiguous position
    # iupac_base = IUPAC_CODES.get(frozenset(most_common_bases), 'N')
    # return iupac_base, (max_count / total_bases * 100)

def get_consensus(reference, sequences, gapped):
    """
    Generate a consensus sequence from a list of DNA sequences.
    Uses reference sequence when all other sequences have gaps.
    Uses IUPAC ambiguity codes for positions with multiple equally frequent bases.
    Returns consensus sequence, scores, and information about reference usage.
    """
    if not sequences:
        return reference, [100.0] * len(reference), [True] * len(reference)

    # Ensure all sequences are the same length
    seq_length = len(reference)
    if not all(len(seq) == seq_length for seq in sequences):
        raise ValueError("All sequences must be the same length")

    consensus = []
    scores = []
    reference_used = []  # Track positions where reference was used

    # Look at each position across all sequences
    for i in range(seq_length):
        # Get all bases at this position
        position_bases = [seq[i] for seq in sequences]

        # Remove gaps from consideration for consensus
        valid_bases = [base for base in position_bases if base != '-']

        if not valid_bases:
            if gapped:
                consensus_base = "-"
                consensus_score = 100.0  # Reference base score
                used_reference = True
            else:
                # If no valid bases in other sequences, use reference
                consensus_base = reference[i]
                consensus_score = 100.0  # Reference base score
                used_reference = True
        else:
            # Count occurrences of each base, excluding gaps
            base_counts = Counter(valid_bases)
            consensus_base, consensus_score = get_consensus_base(base_counts)
            used_reference = False

        consensus.append(consensus_base)
        scores.append(consensus_score)
        reference_used.append(used_reference)

    return ''.join(consensus), scores, reference_used

def main():
    if len(sys.argv) == 1:
        print("Usage: python consensus_fasta.py <input_fasta_file> <gapped>", file=sys.stderr)
        sys.exit(1)

    if len(sys.argv) >2:
        gapped = True
    else:
        gapped = False

    input_file = sys.argv[1]

    try:
        # Read sequences from input file
        reference, reference_header, sequences = read_fasta(input_file)
#        print(sequences)
        if not reference:
            print("No sequences found in input file", file=sys.stderr)
            sys.exit(1)

        # Generate consensus sequence
        consensus, scores, reference_used = get_consensus(reference, sequences, gapped)

        # Write output
#        print(f">{reference_header} length={len(consensus)}")
        print(f">{reference_header}")
        print(consensus)

        # Print additional statistics
        print("\nConsensus analysis at each position:", file=sys.stderr)
        for pos, (base, score, ref_used) in enumerate(zip(consensus, scores, reference_used), 1):
            if ref_used:
                print(f"Position {pos}: {base} (reference base used)", file=sys.stderr)
            elif base in 'ACGT':
                print(f"Position {pos}: {base} ({score:.1f}%)", file=sys.stderr)
            else:
                print(f"Position {pos}: {base} (ambiguous position, highest frequency: {score:.1f}%)", file=sys.stderr)

        # Calculate statistics
        total_ref_positions = sum(reference_used)
        total_positions = len(consensus)

        print(f"\nTotal sequences analyzed: {len(sequences)} (excluding reference)", file=sys.stderr)
        print(f"Reference sequence used in {total_ref_positions} of {total_positions} positions ({(total_ref_positions/total_positions)*100:.1f}%)", file=sys.stderr)

        # Calculate average score excluding reference positions
        non_ref_scores = [score for score, ref_used in zip(scores, reference_used) if not ref_used]
        if non_ref_scores:
            avg_score = sum(non_ref_scores) / len(non_ref_scores)
            print(f"Average consensus score (excluding reference positions): {avg_score:.1f}%", file=sys.stderr)

        # Print IUPAC code reference
        print("\nIUPAC ambiguity codes used:", file=sys.stderr)
        print("R = A/G   Y = C/T   S = G/C   W = A/T   K = G/T   M = A/C", file=sys.stderr)
        print("B = C/G/T   D = A/G/T   H = A/C/T   V = A/C/G   N = A/C/G/T", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
