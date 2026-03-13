# sage_setup: distribution = sagemath-modules

cdef class Network_flow_column:
    cdef Py_ssize_t height
    cdef Py_ssize_t mod
    cdef list idx_map
    cdef list val_map

    cdef Py_hash_t _hash
    cdef bint _hash_valid

    cpdef set_col(self, object col)