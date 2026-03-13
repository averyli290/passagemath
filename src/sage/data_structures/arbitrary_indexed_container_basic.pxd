# sage_setup: distribution = sagemath-categories
"""
Indexed container data structure

AUTHORS:

- Avery Li
"""

cdef class Arbitrary_Indexed_Container_Basic:
    cdef Py_ssize_t _length
    cdef list _value_array
    cdef dict _obj_to_idx

    cdef void _initialize_container(self, list obj_array, list value_array)

    cpdef void append(self, object obj, object value)
    cpdef void extend(self, object obj_array, object value_array)
    cpdef object get_contiguous_array(self)


    # __setitem__(self, object obj) # perhaps default to end if convenient/makes sense
    # obj, value -> key, value 