# sage_setup: distribution = sagemath-categories

cdef class Arbitrary_Indexed_Container_Base:
    def __init__(self, obj_array=None, _value_array=None):
        self._length = 0

    cdef void _initialize_container(self, object obj_array, object _value_array):
        if (obj_array is None) != (_value_array is None):
            raise AssertionError(
                f"Either obj_array and _value_array must both not be None, "
                f"or they must both be None, "
                f"but obj_array has type {type(obj_array)} "
                f"and _value_array has type {type(_value_array)}"
            )

        if obj_array is None or len(obj_array) == 0:
            self._value_array = []
            return

        if len(obj_array) != len(_value_array):
            raise AssertionError(
                f"_length of obj_array must equal _length of _value_array, "
                f"but len(obj_array)={len(obj_array)} "
                f"and len(_value_array)={len(_value_array)}."
            )

        raise NotImplementedError(
            "_initialize_container not implemented for non-empty arrays."
        )

    def __getitem__(self, object idx):
        raise NotImplementedError("__getitem__ not implemented.")

    cpdef void append(self, object value):
        raise NotImplementedError("append not implemented.")

    cpdef void insert(self, Py_ssize_t idx, object value):
        raise NotImplementedError("insert not implemented.")

    cpdef void extend(self, object obj_array, object _value_array):
        raise NotImplementedError("extend not implemented.")

    cpdef object get_contiguous_array(self):
        raise NotImplementedError("get_contiguous_array not implemented.")

    def __len__(self):
        return self._length