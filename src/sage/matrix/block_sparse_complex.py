from bisect import bisect_right
from sage.matrix.constructor import Matrix
from sage.rings.integer_ring import ZZ

class BlockSparseMatrixComplex:
    def __init__(self, shape, ring=ZZ):
        """
        Block-sparse matrix with arbitrary block sizes (non-overlapping).
        shape: (rows, cols) of full matrix
        ring: coefficient ring
        """
        self.shape = shape
        self.ring = ring
        # dict: row_start -> (col_starts, blocks)
        #   col_starts: sorted list of c0
        #   blocks: list of (c0, r1, c1, block) sorted by c0
        self.row_blocks = {}
        self.row_starts = []  # sorted list of row starts for binary search

    def set_block(self, i, j, block):
        """
        Place a block with top-left corner at (i, j).
        """
        nrows, ncols = block.nrows(), block.ncols()

        if i + nrows > self.shape[0] or j + ncols > self.shape[1]:
            raise ValueError("Block exceeds matrix bounds")

        if i not in self.row_blocks:
            self.row_blocks[i] = ([], [])
            self.row_starts.append(i)
            self.row_starts.sort()

        col_starts, blocks = self.row_blocks[i]

        # overlap check (linear within that row band, but still fast)
        for c0, r1, c1, _ in blocks:
            if not (j + ncols <= c0 or c1 <= j):
                raise ValueError("Block overlaps with existing block in same row band")

        # insert in sorted order by c0
        pos = bisect_right(col_starts, j)
        col_starts.insert(pos, j)
        blocks.insert(pos, (j, i + nrows, j + ncols, block))

    def _get_unsafe(self, i, j):
        """
        Return (i,j) entry (no range checking) using two-level binary search.
        """
        # 1. Find row band
        row_idx = bisect_right(self.row_starts, i) - 1
        if row_idx < 0:
            return self.ring(0)

        row_start = self.row_starts[row_idx]
        col_starts, blocks = self.row_blocks[row_start]

        # Quick check: if i not inside any block's row span, skip
        if not blocks or i >= blocks[-1][1]:  # compare to r1 of last block
            return self.ring(0)

        # 2. Find candidate block by col
        col_idx = bisect_right(col_starts, j) - 1
        if col_idx < 0:
            return self.ring(0)

        c0, r1, c1, block = blocks[col_idx]
        if row_start <= i < r1 and c0 <= j < c1:
            return block[i - row_start, j - c0]

        return self.ring(0)

    def __getitem__(self, idx):
        i, j = idx
        if not (0 <= i < self.shape[0] and 0 <= j < self.shape[1]):
            raise IndexError("matrix index out of range")
        return self._get_unsafe(i, j)

    def to_dense(self):
        M = Matrix(self.ring, self.shape[0], self.shape[1])
        for r0, (col_starts, blocks) in self.row_blocks.items():
            for c0, r1, c1, B in blocks:
                for ii in range(r1 - r0):
                    for jj in range(c1 - c0):
                        M[r0 + ii, c0 + jj] = B[ii, jj]
        return M

    def __repr__(self):
        result = []
        for r0, (col_starts, blocks) in self.row_blocks.items():
            for c0, r1, c1, B in blocks:
                result.append(f"Block at ({r0},{c0}), shape=({r1-r0},{c1-c0}):\n{B}\n")
        return "\n".join(result)
