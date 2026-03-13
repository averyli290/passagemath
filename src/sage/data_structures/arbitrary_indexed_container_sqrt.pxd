# sage_setup: distribution = sagemath-categories
"""
Indexed container data structure

AUTHORS:

- Avery Li
"""

cdef class Arbitrary_Indexed_Container_Sqrt:
    cdef Py_ssize_t _length
    cdef dict _key_to_repr_key
    cdef dict _repr_key_to_region_start_idx

    cdef list _regions
    cdef list _regions_key_values
    cdef list _region_size

    cdef Py_ssize_t _sqrtval
    cdef Py_ssize_t _num_regions
    cdef Py_ssize_t _reinitialize_ratio

    cdef void _initialize_container(self, object key_array, object value_array)

    cpdef object get_contiguous_array(self)
    cpdef void append(self, object key, object value)
    cpdef void insert(self, Py_ssize_t idx, object key, object value)
    cdef void _assert_no_duplicate_key(self, object key)
    cpdef void reinitialize_container(self)
    cpdef int index(self, object key)
'''

cdef class Arbitrary_Indexed_Container_Sqrt:
    cdef Py_ssize_t _length
    cdef dict _key_to_position

    cdef list _regions
    cdef list _regions_key_values
    cdef list _region_size

    cdef Py_ssize_t _sqrtval
    cdef Py_ssize_t _num_regions
    cdef Py_ssize_t _reinitialize_ratio

    cdef void _initialize_container(self, object key_array, object value_array)

    cpdef object get_contiguous_array(self)
    cpdef void append(self, object key, object value)
    cpdef void insert(self, Py_ssize_t idx, object key, object value)
    # cdef void _assert_no_duplicate_key(self, object key)
    cpdef void reinitialize_container(self)
    cpdef int index(self, object key)
'''