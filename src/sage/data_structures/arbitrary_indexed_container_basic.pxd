from cpython cimport Py_ssize_t
from sage.data_structures.arbitrary_indexed_container_base cimport (
    Arbitrary_Indexed_Container_Base
)

cdef class Arbitrary_Indexed_Container_Basic(Arbitrary_Indexed_Container_Base):

    cdef dict _obj_to_idx

    cdef void _initialize_container(self, object obj_array, object _value_array)
    cpdef object __getitem__(self, object key)
    cpdef void append(self, object obj, object value)
    cpdef void insert(self, Py_ssize_t idx, object obj, object value)