# sage_setup: distribution = sagemath-categories
"""
Indexed container data structure

AUTHORS:

- Avery Li
"""

cdef class Arbitrary_Indexed_Container_Base:
    cdef Py_ssize_t _length

    cdef void _initialize_container(self, object obj_array, object value_array)

    cpdef void append(self, object value)
    cpdef void insert(self, Py_ssize_t idx, object value)
    cpdef void extend(self, object obj_array, object value_array)
    cpdef object get_contiguous_array(self)