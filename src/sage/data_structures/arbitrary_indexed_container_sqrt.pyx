# sage_setup: distribution = sagemath-categories

cdef class Arbitrary_Indexed_Container_Sqrt:
    def __init__(self, key_array=None, value_array=None):
        r"""
        Create a sparse matrix with 8-bit integer entries implemented in CMR.

        INPUT:

        - ``parent`` -- a matrix space

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``False``, assume without checking that the
          entries lie in the base ring

        - ``immutable`` -- ignored (for backwards compatibility)?

        (TO REPLACE)
        TESTS::

            sage: from sage.matrix.matrix_cmr_sparse import Matrix_cmr_chr_sparse
            sage: M = Matrix_cmr_chr_sparse(MatrixSpace(ZZ, 2, 3, sparse=True),
            ....:                           [[1, 2, 3], [4, 0, 6]]); M
            [1 2 3]
            [4 0 6]
            sage: TestSuite(M).run()
        """

        self._sqrtval = 0
        self._num_regions = 0
        self._reinitialize_ratio = 10

        self._key_to_repr_key = {}
        self._repr_key_to_region_start_idx = {}

        self._regions = []
        self._regions_key_values = []
        self._region_size = []

        self._initialize_container(key_array=key_array, value_array=value_array)

    cdef void _initialize_container(self, object key_array, object value_array):
        """
        Initialize regions and mapping for sqrt-decomposed container
        """
        if (key_array is None) != (value_array is None):
            raise AssertionError("Arrays must both be None or both not None")

        if key_array is None or len(key_array) == 0:
            self._length = 0
            return

        if len(key_array) != len(value_array):
            raise AssertionError("Lengths of key_array and value_array must match")

        cdef Py_ssize_t i
        cdef Py_ssize_t region_length
        self._length = len(key_array)
        region_length = max(10, <Py_ssize_t>(self._length ** 0.5))
        self._sqrtval = region_length
        self._reinitialize_ratio = region_length

        for i in range(self._length):
            key = key_array[i]

            if key in self._key_to_repr_key:
                # Clear and reset state
                self._length = 0
                self._sqrtval = 0
                self._num_regions = 0
                self._regions.clear()
                self._region_size.clear()
                self._regions_key_values.clear()
                self._key_to_repr_key.clear()
                self._repr_key_to_region_start_idx.clear()
                raise AssertionError(f"Duplicate object {key} during initialization")

            if i % region_length == 0:
                self._regions.append([])
                self._regions_key_values.append([])
                self._region_size.append(0)
                self._num_regions += 1
                self._key_to_repr_key[key] = key
                self._repr_key_to_region_start_idx[key] = i
            else:
                self._key_to_repr_key[key] = self._regions_key_values[-1][0]

            self._regions[-1].append(value_array[i])
            self._regions_key_values[-1].append(key)
            self._region_size[-1] += 1

    def __len__(self):
        return self._length

    cdef void _assert_no_duplicate_key(self, object key):
        if key in self._key_to_repr_key:
            raise AssertionError(f"The object {key} already exists.")

    cpdef void reinitialize_container(self):
        """
        Rebuild container after regions grow too large
        """
        cdef Py_ssize_t i, j, idx = 0
        key_array = [None] * self._length
        value_array = [None] * self._length

        for i in range(self._num_regions):
            for j in range(self._region_size[i]):
                key_array[idx] = self._regions_key_values[i][j]
                value_array[idx] = self._regions[i][j]
                idx += 1

        # Reset container
        self._length = 0
        self._sqrtval = 0
        self._num_regions = 0
        self._regions.clear()
        self._region_size.clear()
        self._regions_key_values.clear()
        self._key_to_repr_key.clear()
        self._repr_key_to_region_start_idx.clear()

        self._initialize_container(key_array, value_array)

    cpdef int index(self, object key):
        """
        Return the linear index for object key using sqrt decomposition.
        """
        if key not in self._key_to_repr_key:
            raise LookupError(f"Representative object of {key} not found")

        repr_key = self._key_to_repr_key[key]
        region_start = self._repr_key_to_region_start_idx[repr_key]

        # locate region
        for i in range(self._num_regions):
            if self._repr_key_to_region_start_idx[self._key_to_repr_key[self._regions_key_values[i][0]]] == region_start:
                return self._regions_key_values[i].index(key)
        raise LookupError(f"object {key} not found in any region")

    def __setitem__(self, object key, object value):
        """
        If key already exists, set the item at the index represented by key to value.
        Otherwise, append the item.
        """
        if key in self._key_to_repr_key:
            # Find key
            repr_key = self._key_to_repr_key[key]
            region_start = self._repr_key_to_region_start_idx[repr_key]

            # locate region and set value when key found
            for i in range(self._num_regions):
                if self._repr_key_to_region_start_idx[self._key_to_repr_key[self._regions_key_values[i][0]]] == region_start:
                    for j in range(self._region_size[i]):
                        if self._regions_key_values[i][j] == key:
                            self._regions[i][j] = value
                            return
            raise LookupError(f"key {key} not found in any region, but it has a representative key")
        else:
            self.append(key, value)

    def __getitem__(self, object key):
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
            
            if key not in self._key_to_repr_key:
                raise LookupError(f"Representative object of {key} not found")

            repr_key = self._key_to_repr_key[key]
            region_start = self._repr_key_to_region_start_idx[repr_key]

            # locate region
            for i in range(self._num_regions):
                if self._repr_key_to_region_start_idx[self._key_to_repr_key[self._regions_key_values[i][0]]] == region_start:
                    # search within region
                    for j in range(self._region_size[i]):
                        if self._regions_key_values[i][j] == key:
                            return self._regions[i][j]
            raise LookupError(f"object {key} not found in any region")

    cpdef object get_contiguous_array(self):
        """
        Return flattened array and mapping from key to linear index
        """
        cdef Py_ssize_t idx, i, j
        idx = 0
        array = [None] * self._length
        mapping = {}

        for i in range(self._num_regions):
            for j in range(self._region_size[i]):
                array[idx] = self._regions[i][j]
                mapping[self._regions_key_values[i][j]] = idx
                idx += 1
        return array, mapping

    cpdef void append(self, object key, object value):
        self._assert_no_duplicate_key(key)

        if len(self._regions_key_values) == 0:
            self._regions.append([])
            self._regions_key_values.append([])
            self._region_size.append(0)
            self._num_regions += 1

        if len(self._regions_key_values[-1]) == 0:
            self._key_to_repr_key[key] = key
        else:
            self._key_to_repr_key[key] = self._key_to_repr_key[self._regions_key_values[-1][0]]

        self._regions[-1].append(value)
        self._regions_key_values[-1].append(key)
        self._region_size[-1] += 1
        self._length += 1

        if self._region_size[-1] > self._reinitialize_ratio * self._sqrtval:
            self.reinitialize_container() 

    cpdef void insert(self, Py_ssize_t idx, object key, object value):
        """
        Insert key,value at linear index idx
        """
        cdef Py_ssize_t region_idx, insert_idx, region_start_idx
        if idx < 0 or idx > self._length:
            raise AssertionError(f"Cannot insert at index {idx}")

        if idx == self._length:
            self.append(key, value)
            return

        self._assert_no_duplicate_key(key)

        region_idx = 0
        while region_idx < self._num_regions:
            region_start_idx = self._repr_key_to_region_start_idx[self._key_to_repr_key[self._regions_key_values[region_idx][0]]]
            if region_start_idx <= idx < region_start_idx + self._region_size[region_idx]:
                insert_idx = idx - region_start_idx
                self._key_to_repr_key[key] = self._key_to_repr_key[self._regions_key_values[region_idx][0]]
                self._regions[region_idx].insert(insert_idx, value)
                self._regions_key_values[region_idx].insert(insert_idx, key)
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
            repr_key = self._key_to_repr_key[self._regions_key_values[region_idx][0]]
            self._repr_key_to_region_start_idx[repr_key] += 1
            region_idx += 1

        if self._region_size[target_region_idx] > self._reinitialize_ratio * self._sqrtval:
            self.reinitialize_container()
'''
# sage_setup: distribution = sagemath-categories

cdef class Arbitrary_Indexed_Container_Sqrt:

    def __init__(self, key_array=None, value_array=None):

        self._sqrtval = 0
        self._num_regions = 0
        self._reinitialize_ratio = 10

        self._key_to_position = {}

        self._regions = []
        self._regions_key_values = []
        self._region_size = []

        self._initialize_container(key_array=key_array, value_array=value_array)

    cdef void _initialize_container(self, object key_array, object value_array):

        if (key_array is None) != (value_array is None):
            raise AssertionError("Arrays must both be None or both not None")

        if key_array is None or len(key_array) == 0:
            self._length = 0
            return

        if len(key_array) != len(value_array):
            raise AssertionError("Lengths must match")

        cdef Py_ssize_t i, region_length

        self._length = len(key_array)

        region_length = max(10, <Py_ssize_t>(self._length ** 0.5))
        self._sqrtval = region_length
        self._reinitialize_ratio = region_length

        for i in range(self._length):

            if i % region_length == 0:
                self._regions.append([])
                self._regions_key_values.append([])
                self._region_size.append(0)
                self._num_regions += 1

            key = key_array[i]

            if key in self._key_to_position:
                raise AssertionError(f"Duplicate key {key}")

            region_idx = self._num_regions - 1
            offset = self._region_size[region_idx]

            self._regions[region_idx].append(value_array[i])
            self._regions_key_values[region_idx].append(key)
            self._region_size[region_idx] += 1

            self._key_to_position[key] = (region_idx, offset)

    def __len__(self):
        return self._length

    cpdef void reinitialize_container(self):

        cdef Py_ssize_t i, j, idx = 0

        key_array = [None] * self._length
        value_array = [None] * self._length

        for i in range(self._num_regions):
            for j in range(self._region_size[i]):
                key_array[idx] = self._regions_key_values[i][j]
                value_array[idx] = self._regions[i][j]
                idx += 1

        self._length = 0
        self._sqrtval = 0
        self._num_regions = 0

        self._regions.clear()
        self._region_size.clear()
        self._regions_key_values.clear()
        self._key_to_position.clear()

        self._initialize_container(key_array, value_array)

    cpdef int index(self, object key):

        if key not in self._key_to_position:
            raise LookupError(f"{key} not found")

        region_idx, offset = self._key_to_position[key]

        cdef Py_ssize_t i, idx = 0

        for i in range(region_idx):
            idx += self._region_size[i]

        return idx + offset

    def __getitem__(self, object key):

        if key not in self._key_to_position:
            raise LookupError(f"{key} not found")

        region_idx, offset = self._key_to_position[key]
        return self._regions[region_idx][offset]

    def __setitem__(self, object key, object value):

        if key in self._key_to_position:

            region_idx, offset = self._key_to_position[key]
            self._regions[region_idx][offset] = value
            return

        self.append(key, value)

    cpdef object get_contiguous_array(self):

        cdef Py_ssize_t idx = 0
        cdef Py_ssize_t i, j

        array = [None] * self._length
        mapping = {}

        for i in range(self._num_regions):
            for j in range(self._region_size[i]):
                key = self._regions_key_values[i][j]
                array[idx] = self._regions[i][j]
                mapping[key] = idx
                idx += 1

        return array, mapping

    cpdef void append(self, object key, object value):

        if key in self._key_to_position:
            raise AssertionError(f"{key} already exists")

        if len(self._regions) == 0:

            self._regions.append([])
            self._regions_key_values.append([])
            self._region_size.append(0)
            self._num_regions += 1

        region_idx = self._num_regions - 1
        offset = self._region_size[region_idx]

        self._regions[region_idx].append(value)
        self._regions_key_values[region_idx].append(key)

        self._key_to_position[key] = (region_idx, offset)

        self._region_size[region_idx] += 1
        self._length += 1

        if self._region_size[region_idx] > self._reinitialize_ratio * self._sqrtval:
            self.reinitialize_container()

    cpdef void insert(self, Py_ssize_t idx, object key, object value):

        cdef Py_ssize_t region_idx = 0
        cdef Py_ssize_t region_start = 0
        cdef Py_ssize_t insert_offset

        if idx < 0 or idx > self._length:
            raise AssertionError("Invalid index")

        if key in self._key_to_position:
            raise AssertionError(f"{key} already exists")

        if idx == self._length:
            self.append(key, value)
            return

        while region_idx < self._num_regions:

            if region_start + self._region_size[region_idx] > idx:
                break

            region_start += self._region_size[region_idx]
            region_idx += 1

        insert_offset = idx - region_start

        self._regions[region_idx].insert(insert_offset, value)
        self._regions_key_values[region_idx].insert(insert_offset, key)

        self._region_size[region_idx] += 1
        self._length += 1

        # update offsets inside region
        cdef Py_ssize_t j
        for j in range(insert_offset, self._region_size[region_idx]):

            k = self._regions_key_values[region_idx][j]
            self._key_to_position[k] = (region_idx, j)      # Update offset to new index

        if self._region_size[region_idx] > self._reinitialize_ratio * self._sqrtval:
            self.reinitialize_container()
'''