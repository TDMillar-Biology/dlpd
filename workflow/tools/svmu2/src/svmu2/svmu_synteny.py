"""
svmu2_synteny.py

Synteny-aware PAM orthologue assignment for SVMU2.
Uses SVMU2/delta_parser refined synteny blocks to project reference PAMs
into query genomes and determine one-to-one positional orthology.

Author: Trevor Millar (with ChatGPT)
"""

from dataclasses import dataclass
from typing import List, Optional


# ------------------------------------------------------------
# 1. Lightweight wrapper around delta_parser AlignmentBlock
# ------------------------------------------------------------

## A lighter weight representation of the delta parser blocks
@dataclass
class SyntenyBlock:
    ref_chr: str
    ref_start: int
    ref_end: int

    query_chr: str
    query_start: int
    query_end: int

    slope: float
    block_id: str
    original: object  # for debugging

    def __post_init__(self):
        self.orientation = 1 if self.slope >= 0 else -1

    def contains_ref(self, pos: int) -> bool:
        return self.ref_start <= pos <= self.ref_end

    def contains_query(self, pos: int, tol: int = 1000) -> bool:
        qmin = min(self.query_start, self.query_end) - tol
        qmax = max(self.query_start, self.query_end) + tol
        return qmin <= pos <= qmax

    def project_ref_to_query(self, ref_pos: int) -> float:
        offset = ref_pos - self.ref_start
        if self.orientation == 1:
            return self.query_start + offset * self.slope
        else:
            return self.query_start - offset * self.slope

# ------------------------------------------------------------
# 2. Build synteny block list from delta_parser refined_path
# ------------------------------------------------------------

## build list of synteny blocks from input list of delta parser blocks
def build_synteny_blocks(refined_path: List[object]) -> List[SyntenyBlock]:
    blocks = []
    
    for i, blk in enumerate(refined_path):
        b = SyntenyBlock(
            ref_chr    = blk.reference,
            ref_start  = blk.reference_start,
            ref_end    = blk.reference_end,

            query_chr  = blk.query,
            query_start = blk.query_start,
            query_end   = blk.query_end,

            slope      = blk.slope,
            block_id = blk.reference + "_" + str(i),
            original   = blk

        )
        blocks.append(b)
    
    return blocks

# ------------------------------------------------------------
# 3. Find block containing a reference coordinate
# ------------------------------------------------------------
## missing the case of multiple blocks containing the ref
def find_block_for_ref(blocks: List[SyntenyBlock], ref_pos: int) -> Optional[SyntenyBlock]: #optional is eq to union(x, None)
    for b in blocks:
        if b.contains_ref(ref_pos):
            return b
    return None


# ------------------------------------------------------------
# 4. Find syntenic orthologue
# ------------------------------------------------------------

@dataclass
class PAMResult:
    status: str                # SYNTENIC / LOST / NO_BLOCK / NON_SYNTENIC / MULTI
    best: Optional[int]        # best query coord or None
    expected: Optional[float]  # float projected position
    multi_hits: Optional[List[int]] = None


def find_syntenic_pam(block: SyntenyBlock,
                      ref_pos: int,
                      query_hits: List[int],
                      tol: int = 1000) -> PAMResult:
    """
    choose the syntenic hit in the block closest to projected position.
    """

    # 1. Expected query coordinate
    expected = block.project_ref_to_query(ref_pos)

    # 2. Filter hits within the block (+/- tol)
    candidates = [q for q in query_hits if block.contains_query(q, tol=tol)] # one to many safe

    if len(candidates) == 0:
        return PAMResult(status="LOST", best=None, expected=expected)

    # 3. Geometric closest-to-expected selection
    diffs = [(abs(q - expected), q) for q in candidates]
    diffs.sort(key=lambda x: x[0])
    best_q = diffs[0][1]

    # Optional: track multi-hit info
    multi = [int(round(q)) for q in candidates] if len(candidates) > 1 else None

    return PAMResult(
        status="SYNTENIC",
        best=int(round(best_q)),
        expected=expected,
        multi_hits=multi
    )


# ------------------------------------------------------------
# 5. Top-level convenience function
# ------------------------------------------------------------

def classify_pam(ref_pos: int,
                 query_hits: List[int],
                 blocks: List[SyntenyBlock],
                 tol: int = 1000) -> PAMResult:
    """
    High-level call:
      - finds the block
      - performs synteny-aware orthologue detection
    """

    block = find_block_for_ref(blocks, ref_pos)

    if block is None:
        return PAMResult(status="NO_BLOCK", best=None, expected=None)

    return find_syntenic_pam(block, ref_pos, query_hits, tol=tol)

def find_block_for_ref(blocks, chrom, pos):
    # return *all* matching block IDs
    return [b.block_id for b in blocks
            if b.ref_chr == chrom and b.contains_ref(pos)]

def find_block_for_query(blocks, chrom, pos):
    return [b.block_id for b in blocks
            if b.query_chr == chrom and b.contains_query(pos)]


