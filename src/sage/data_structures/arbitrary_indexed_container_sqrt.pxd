from cpython cimport Py_ssize_t
from sage.data_structures.arbitrary_indexed_container_base cimport Arbitrary_Indexed_Container_Base


cdef class Arbitrary_Indexed_Container_Sqrt(Arbitrary_Indexed_Container_Base):

    cdef dict _obj_to_repr_obj
    cdef dict _repr_obj_to_region_start_idx

    cdef list _regions
    cdef list _regions_obj_values
    cdef list _region_size

    cdef Py_ssize_t _sqrtval
    cdef Py_ssize_t _num_regions
    cdef Py_ssize_t _reinitialize_ratio

    cdef void _initialize_container(self, object obj_array, object _value_array)

    cpdef object __getitem__(self, object key)
    cpdef object get_contiguous_array(self)
    cpdef void append(self, object obj, object value)
    cpdef void insert(self, Py_ssize_t idx, object obj, object value)
    cdef void _check_duplicate_obj(self, object obj)
    cdef void reinitialize_container(self)
    cdef int index(self, object obj)