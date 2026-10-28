import numpy as np

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """

    # keeping only the positions where both sequences have a real base (skips gaps "-" and unknown "N")
    valid = set("ACGT")
    pairs = [(r, d) for r, d in zip(reference_sequence, distant_sequence) if r in valid and d in valid]
    length = len(pairs)
    mismatches = sum(1 for r, d in pairs if r != d)

    p = mismatches / length # proportion of sites that differ

    # jukes-cantor correction (per site), then scale back to the same unit as the raw count
    distance_per_site = -3 / 4 * np.log(1 - 4 / 3 * p)

    return distance_per_site * length
          


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    raise NotImplementedError()
