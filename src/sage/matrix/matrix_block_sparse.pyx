# sage_setup: distribution = sagemath-modules

from bisect import bisect_right
from cysignals.memory cimport check_allocarray, sig_free
cimport sage.matrix.matrix_sparse as matrix_sparse
from sage.matrix.constructor import Matrix
from sage.matrix.args cimport MatrixArgs_init

cdef class Matrix_block_sparse(matrix_sparse.Matrix_sparse):
    """
    Level 1 block-sparse matrix implementation using binary search.
    """

    def __cinit__(self):
        # dict: row_start -> (col_starts, blocks)
        #   col_starts: sorted list of c0
        #   blocks: list of (c0, r1, c1, block) sorted by c0
        self.row_blocks = {}
        self.row_starts = []
        self._zero = None

    def __init__(self, parent, entries, copy, bint coerce):
        MatrixArgs_init(parent, entries)
        self._zero = self._base_ring.zero()
        self.row_blocks = {}
        self.row_starts = []

        if entries is None:
            return

        if isinstance(entries, dict):
            # keys are (r0, c0), values are blocks
            for (r0, c0), block in entries.items():
                self._insert_block(r0, c0, block)
            return

        if isinstance(entries, Matrix_block_sparse):
            self.row_blocks = dict(entries.row_blocks)
            self.row_starts = list(entries.row_starts)
            return

        raise TypeError("Unsupported entries type for Matrix_block_sparse: %s" % type(entries))

    def __dealloc__(self):
        try:
            if self.row_blocks is not None:
                self.row_blocks.clear()
            if self.row_starts is not None:
                self.row_starts[:] = []
            self._zero = None
        except Exception:
            pass

    cdef void _insert_block(self, Py_ssize_t r0, Py_ssize_t c0, object block):
        # This stuff kinda chucks efficiency out of the window for now...
        cdef Py_ssize_t nrows = block.nrows()
        cdef Py_ssize_t ncols = block.ncols()
        cdef Py_ssize_t r1 = r0 + nrows
        cdef Py_ssize_t c1 = c0 + ncols

        # Find row bands that might overlap
        start_idx = bisect_right(self.row_starts, r0) - 1
        if start_idx < 0:
            start_idx = 0
        end_idx = bisect_right(self.row_starts, r1 - 1)

        # Check overlaps for all row bands that intersect [r0, r1)
        for idx in range(start_idx, end_idx):
            r_band = self.row_starts[idx]
            col_starts, blocks = self.row_blocks[r_band]
            for c0_existing, r1_existing, c1_existing, _ in blocks:
                # check if rectangles overlap
                if not (r1 <= r_band or r1_existing <= r0 or c1 <= c0_existing or c1_existing <= c0):
                    raise ValueError(f"Block overlaps with existing block at row {r_band}, col {c0_existing}")
        
        # Keep dims immutable
        # reinit matrix space + return new object on resize
        # reinitialize with prev row_blocks and row_starts (or generate entries dict if reasonable)

        # Resize if necessary
        if r1 > self._nrows:
            self._nrows = r1
        if c1 > self._ncols:
            self._ncols = c1

        # Insert new row band if needed
        if r0 not in self.row_blocks:
            self.row_blocks[r0] = ([], [])
            posr = bisect_right(self.row_starts, r0)
            self.row_starts.insert(posr, r0)

        col_starts, blocks = self.row_blocks[r0]

        # Insert in sorted order by column
        posc = bisect_right(col_starts, c0)
        col_starts.insert(posc, c0)
        blocks.insert(posc, (c0, r1, c1, block))

    cpdef set_block(self, Py_ssize_t r0, Py_ssize_t c0, object block):
        self._insert_block(r0, c0, block)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Return (i,j) entry using two-level binary search.
        """

        # check over a range starting from row_idx_min to row_idx_max found using bisect
        row_idx_min = bisect_right(self.row_starts, i) - 1
        if row_idx_min < 0:
            return self._zero

        row_idx_max = bisect_right(self.row_starts, i) + 1

        for r0 in range(row_idx_min, row_idx_max):
            if r0 in self.row_blocks:
                col_starts, blocks = self.row_blocks[r0]

                col_idx = bisect_right(col_starts, j) - 1
                if col_idx < 0:
                    return self._zero

                c0, r1, c1, block = blocks[col_idx]
                if r0 <= i < r1 and c0 <= j < c1:
                    return block[i - r0, j - c0]

        return self._zero

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        """
        Set (i,j) entry using two-level binary search.
        """

        # Binary search for the row that the block starts on
        row_idx = bisect_right(self.row_starts, i) - 1
        if row_idx >= 0:
            r0 = self.row_starts[row_idx]
            col_starts, blocks = self.row_blocks[r0]

            # binary search within the row for the col the block starts on
            col_idx = bisect_right(col_starts, j) - 1
            if col_idx >= 0:
                c0, r1, c1, block = blocks[col_idx]
                if r0 <= i < r1 and c0 <= j < c1:
                    block[i - r0, j - c0] = x
                    return

        # Not inside existing block, create 1x1 block
        B = Matrix(self._base_ring, 1, 1)
        B[0, 0] = x
        self._insert_block(i, j, B)

    def __richcmp__(self, other, int op):
        if op != 2:
            return NotImplemented
        if not isinstance(other, Matrix_block_sparse):
            return False
        if self._nrows != other._nrows or self._ncols != other._ncols:
            return False
        if set(self.row_blocks.keys()) != set(other.row_blocks.keys()):
            return False
        for key in self.row_blocks:
            if self.row_blocks[key] != other.row_blocks[key]:
                return False
        return True


