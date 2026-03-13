# sage_setup: distribution = sagemath-categories
"""
Indexed container data structure

AUTHORS:

- Avery Li
"""

cdef class Arbitrary_Indexed_Container_Sqrt:
    cdef Py_ssize_t _length
    cdef dict _key_to_position

    cdef list _region_values
    cdef list _region_keys
    cdef list _region_size

    cdef Py_ssize_t _sqrtval
    cdef Py_ssize_t _num_regions
    cdef Py_ssize_t _reinitialize_ratio

    cdef bint _lazy_reinitialize_flag

    cdef void _initialize_container(self, list key_array, list value_array)

    cpdef object get_contiguous_array(self)
    cpdef void append(self, object key, object value)
    cpdef void extend(self, list key_array, list value_array)
    cpdef void insert(self, Py_ssize_t idx, object key, object value)
    cpdef void reinitialize_container(self)
    cpdef int index(self, object key)