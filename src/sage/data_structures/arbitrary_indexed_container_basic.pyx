# sage_setup: distribution = sagemath-categories

from cpython cimport Py_ssize_t
from sage.data_structures.arbitrary_indexed_container_base cimport (
    Arbitrary_Indexed_Container_Base
)

cdef class Arbitrary_Indexed_Container_Basic(Arbitrary_Indexed_Container_Base):

    def __init__(self, obj_array=None, _value_array=None):
        self._obj_to_idx = {}
        super().__init__(obj_array=obj_array,
                         _value_array=_value_array)

    cdef void _initialize_container(
            self,
            object obj_array,
            object _value_array):

        if (obj_array is None) != (_value_array is None):
            raise AssertionError("obj_array and _value_array must both be None or both not None")

        if obj_array is None or len(obj_array) == 0:
            self._value_array = []
            self._length = 0
            return

        if len(obj_array) != len(_value_array):
            raise AssertionError("Lengths must match")

        cdef Py_ssize_t i
        self._length = len(obj_array)
        self._value_array = list(_value_array)

        for i in range(self._length):
            self._obj_to_idx[obj_array[i]] = i

    cpdef object __getitem__(self, object key):

        cdef Py_ssize_t idx

        if isinstance(key, slice):
            raise NotImplementedError("Slicing not yet implemented for Basic container.")

        if key not in self._obj_to_idx:
            raise AssertionError(f"Object {key} not found.")

        idx = self._obj_to_idx[key]

        if idx < 0 or idx >= self._length:
            raise AssertionError(f"Index {idx} out of bounds.")

        return self._value_array[idx]

    cpdef void append(self, object obj, object value):

        if obj in self._obj_to_idx:
            raise AssertionError(f"Object {obj} already assigned.")

        self._obj_to_idx[obj] = self._length
        self._value_array.append(value)
        self._length += 1

    cpdef void insert(self,
                      Py_ssize_t idx,
                      object obj,
                      object value):
        raise NotImplementedError("Insert not implemented for Basic container.")