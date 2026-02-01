# sage_setup: distribution = sagemath-modules

cimport sage.matrix.matrix_sparse as matrix_sparse

cdef class Matrix_block_sparse(matrix_sparse.Matrix_sparse):
    cdef dict row_blocks
    cdef list row_starts
    cdef object _zero


    cdef void _insert_block(self, Py_ssize_t r0, Py_ssize_t c0, object block)
    cpdef set_block(self, Py_ssize_t r0, Py_ssize_t c0, object block)