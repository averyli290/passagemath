# sage_setup: distribution = sagemath-categories

from cpython cimport Py_ssize_t
from sage.data_structures.arbitrary_indexed_container_base cimport Arbitrary_Indexed_Container_Base

cdef class Arbitrary_Indexed_Container_Sqrt(Arbitrary_Indexed_Container_Base):

    def __init__(self, obj_array=None, _value_array=None):

        self._sqrtval = 0
        self._num_regions = 0
        self._reinitialize_ratio = 10

        self._obj_to_repr_obj = {}
        self._repr_obj_to_region_start_idx = {}

        self._regions = []
        self._regions_obj_values = []
        self._region_size = []

        super().__init__(obj_array=obj_array,
                         _value_array=_value_array)

    cdef void _initialize_container(self, object obj_array, object _value_array):
        """
        Initialize regions and mapping for sqrt-decomposed container
        """
        if (obj_array is None) != (_value_array is None):
            raise AssertionError("Arrays must both be None or both not None")

        if obj_array is None or len(obj_array) == 0:
            self._length = 0
            return

        if len(obj_array) != len(_value_array):
            raise AssertionError("Lengths of obj_array and _value_array must match")

        cdef Py_ssize_t i
        cdef Py_ssize_t region_length
        self._length = len(obj_array)
        region_length = max(10, <Py_ssize_t>(self._length ** 0.5))
        self._sqrtval = region_length
        self._reinitialize_ratio = region_length

        for i in range(self._length):
            obj = obj_array[i]

            if obj in self._obj_to_repr_obj:
                # Clear and reset state
                self._length = 0
                self._sqrtval = 0
                self._num_regions = 0
                self._regions.clear()
                self._region_size.clear()
                self._regions_obj_values.clear()
                self._obj_to_repr_obj.clear()
                self._repr_obj_to_region_start_idx.clear()
                raise AssertionError(f"Duplicate object {obj} during initialization")

            if i % region_length == 0:
                self._regions.append([])
                self._regions_obj_values.append([])
                self._region_size.append(0)
                self._num_regions += 1
                self._obj_to_repr_obj[obj] = obj
                self._repr_obj_to_region_start_idx[obj] = i
            else:
                self._obj_to_repr_obj[obj] = self._regions_obj_values[-1][0]

            self._regions[-1].append(_value_array[i])
            self._regions_obj_values[-1].append(obj)
            self._region_size[-1] += 1

    cdef void _check_duplicate_obj(self, object obj):
        if obj in self._obj_to_repr_obj:
            raise AssertionError(f"The object {obj} already exists.")

    cdef void reinitialize_container(self):
        """
        Rebuild container after regions grow too large
        """
        cdef Py_ssize_t i, j, idx = 0
        obj_array = [None] * self._length
        value_array = [None] * self._length

        for i in range(self._num_regions):
            for j in range(self._region_size[i]):
                obj_array[idx] = self._regions_obj_values[i][j]
                value_array[idx] = self._regions[i][j]
                idx += 1

        # Reset container
        self._length = 0
        self._sqrtval = 0
        self._num_regions = 0
        self._regions.clear()
        self._region_size.clear()
        self._regions_obj_values.clear()
        self._obj_to_repr_obj.clear()
        self._repr_obj_to_region_start_idx.clear()

        self._initialize_container(obj_array, value_array)

    cdef Py_ssize_t index(self, object obj):
        """
        Return the linear index for object obj using sqrt decomposition.
        """
        if obj not in self._obj_to_repr_obj:
            raise LookupError(f"Representative object of {obj} not found")

        repr_obj = self._obj_to_repr_obj[obj]
        region_start = self._repr_obj_to_region_start_idx[repr_obj]

        # locate region
        for i in range(self._num_regions):
            if self._repr_obj_to_region_start_idx[self._obj_to_repr_obj[self._regions_obj_values[i][0]]] == region_start:
                return self._regions_obj_values[i].index(obj)
        raise LookupError(f"Object {obj} not found in any region")

    cpdef object __getitem__(self, object key):
        """
        Get value by object or slice
        """
        cdef Py_ssize_t start_idx, stop_idx, idx, cur_region_idx, num_remaining
        result = []

        if isinstance(key, slice):
            # No stepped slicing
            if key.step is not None:
                raise TypeError("Stepped slicing not supported")

            start_idx = 0 if key.start is None else self.index(key.start)
            stop_idx = self._length if key.stop is None else self.index(key.stop)

            idx = start_idx
            cur_region_idx = 0
            while idx + self._region_size[cur_region_idx] < stop_idx:
                result.extend(self._regions[cur_region_idx])
                idx += self._region_size[cur_region_idx]
                cur_region_idx += 1

            if idx < stop_idx:
                num_remaining = stop_idx - idx
                result.extend(self._regions[cur_region_idx][:num_remaining])
            return result

        else:
            obj = key
            if obj not in self._obj_to_repr_obj:
                raise LookupError(f"Representative object of {obj} not found")

            repr_obj = self._obj_to_repr_obj[obj]
            region_start = self._repr_obj_to_region_start_idx[repr_obj]

            # locate region
            for i in range(self._num_regions):
                if self._repr_obj_to_region_start_idx[self._obj_to_repr_obj[self._regions_obj_values[i][0]]] == region_start:
                    # search within region
                    for j in range(self._region_size[i]):
                        if self._regions_obj_values[i][j] == obj:
                            return self._regions[i][j]
            raise LookupError(f"Object {obj} not found in any region")

    cpdef object get_contiguous_array(self):
        """
        Return flattened array and mapping from obj to linear index
        """
        cdef Py_ssize_t idx, i, j
        idx = 0
        array = [None] * self._length
        mapping = {}

        for i in range(self._num_regions):
            for j in range(self._region_size[i]):
                array[idx] = self._regions[i][j]
                mapping[self._regions_obj_values[i][j]] = idx
                idx += 1
        return array, mapping

    cpdef void append(self, object obj, object value):
        self._check_duplicate_obj(obj)
        self._obj_to_repr_obj[obj] = self._obj_to_repr_obj[self._regions_obj_values[-1][0]]

        self._regions[-1].append(value)
        self._regions_obj_values[-1].append(obj)
        self._region_size[-1] += 1
        self._length += 1

        if self._region_size[-1] > self._reinitialize_ratio * self._sqrtval:
            self.reinitialize_container()

    cpdef void insert(self, Py_ssize_t idx, object obj, object value):
        """
        Insert obj,value at linear index idx
        """
        cdef Py_ssize_t region_idx, insert_idx, region_start_idx
        if idx < 0 or idx > self._length:
            raise AssertionError(f"Cannot insert at index {idx}")

        if idx == self._length:
            self.append(obj, value)
            return

        self._check_duplicate_obj(obj)

        region_idx = 0
        while region_idx < self._num_regions:
            region_start_idx = self._repr_obj_to_region_start_idx[self._obj_to_repr_obj[self._regions_obj_values[region_idx][0]]]
            if region_start_idx <= idx < region_start_idx + self._region_size[region_idx]:
                insert_idx = idx - region_start_idx
                self._obj_to_repr_obj[obj] = self._obj_to_repr_obj[self._regions_obj_values[region_idx][0]]
                self._regions[region_idx].insert(insert_idx, value)
                self._regions_obj_values[region_idx].insert(insert_idx, obj)
                self._region_size[region_idx] += 1
                self._length += 1
                break
            region_idx += 1

        target_region_idx = region_idx
        if target_region_idx >= self._num_regions:
            raise IndexError(f"Could not insert value at index {idx}")

        # shift region start indices for subsequent regions
        region_idx += 1
        while region_idx < self._num_regions:
            repr_obj = self._obj_to_repr_obj[self._regions_obj_values[region_idx][0]]
            self._repr_obj_to_region_start_idx[repr_obj] += 1
            region_idx += 1

        if self._region_size[target_region_idx] > self._reinitialize_ratio * self._sqrtval:
            self.reinitialize_container()