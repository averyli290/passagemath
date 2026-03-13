# sage_setup: distribution = sagemath-categories


"""
Possible changes
- prefix sums (depending on access/update frequencies)
"""


cdef class Arbitrary_Indexed_Container_Sqrt:

    def __init__(self, key_array=None, value_array=None):
        """
        Initialize a sqrt-decomposed indexed container.

        The container stores key/value pairs while maintaining the order
        of insertion. Internally the data is partitioned into regions of
        size approximately ``sqrt(n)`` to allow efficient range queries
        and insertions.

        INPUT:

        - ``key_array`` -- list of keys (optional)
        - ``value_array`` -- list of values corresponding to ``key_array``

        If both arrays are provided, the container is initialized with
        those entries.

        EXAMPLES::

            sage: from sage.data_structures.arbitrary_indexed_container_sqrt import Arbitrary_Indexed_Container_Sqrt
            sage: keys = [(0,0),(1,1),(2,2)]
            sage: vals = [10,20,30]
            sage: c = Arbitrary_Indexed_Container_Sqrt(keys, vals)
            sage: len(c)
            3

        TESTS::

            sage: c = Arbitrary_Indexed_Container_Sqrt()
            sage: len(c)
            0
        """

        self._sqrtval = 0
        self._num_regions = 0
        self._reinitialize_ratio = 10

        self._key_to_position = {}

        self._region_values = []                  # 2d array storing values
        self._region_keys = []              # 2d array storing keys
        self._region_size = []

        # Whether or not the container needs to be reinitialized
        # Reinitialization required when container imbalanced and
        # and on __getitem__(slice) or insert
        # append, _initialize_container, reinitialize_container, __len__, index, __getitem__(index)
        self._lazy_reinitialize_flag = False 

        self._initialize_container(key_array=key_array, value_array=value_array)

    cdef void _initialize_container(self, list key_array, list value_array):
        """
        Internal method to construct the container from arrays.

        The data is partitioned into regions of size approximately
        ``sqrt(n)`` and mappings from keys to their region position
        are created.

        INPUT:

        - ``key_array`` -- list of keys
        - ``value_array`` -- list of values

        EXAMPLES::

            sage: keys = [(0,0),(1,1),(2,2)]
            sage: vals = [1,2,3]
            sage: c = Arbitrary_Indexed_Container_Sqrt(keys, vals)
            sage: c[(1,1)]
            2

        TESTS::

            sage: keys = [(0,0),(0,0)]
            sage: vals = [1,2]
            sage: Arbitrary_Indexed_Container_Sqrt(keys, vals)
            Traceback (most recent call last):
            ...
            AssertionError: Key array and value array lengths must match
        """

        if (key_array is None) != (value_array is None):
            raise AssertionError("Arrays must both be None or both not None")

        if key_array is None or len(key_array) == 0:
            self._length = 0
            return

        if len(key_array) != len(value_array):
            raise AssertionError("Key array and value array lengths must match")

        cdef Py_ssize_t i, region_length

        self._length = len(key_array)

        region_length = max(10, <Py_ssize_t>(self._length ** 0.5))
        self._sqrtval = region_length
        self._reinitialize_ratio = region_length

        for i in range(self._length):

            if i % region_length == 0:
                self._region_values.append([])
                self._region_keys.append([])
                self._region_size.append(0)
                self._num_regions += 1

            key = key_array[i]

            if key in self._key_to_position:
                raise AssertionError(f"Duplicate key {key}")

            region_idx = self._num_regions - 1
            offset = self._region_size[region_idx]

            self._region_values[region_idx].append(value_array[i])
            self._region_keys[region_idx].append(key)
            self._region_size[region_idx] += 1

            self._key_to_position[key] = (region_idx, offset)
        
        self._lazy_reinitialize_flag = False

    def __len__(self):
        return self._length

    cpdef void reinitialize_container(self):
        """
        Rebuild the internal region structure.

        This method reconstructs the container so that regions are
        balanced again. It is used when regions become too large
        after many insertions or when accessing and the lazy
        rebuild flag has been set.

        EXAMPLES::

            sage: c = Arbitrary_Indexed_Container_Sqrt([(0,0)], [1])
            sage: for i in range(20):
            ....:     c.append((i+1,i+1), i)
            sage: c.reinitialize_container()
            sage: len(c) > 0
            True

        TESTS::

            sage: c = Arbitrary_Indexed_Container_Sqrt([(0,0)], [1])
            sage: c.reinitialize_container()
            sage: len(c)
            1
        """

        cdef Py_ssize_t i, j, idx = 0

        key_array = [None] * self._length
        value_array = [None] * self._length

        for i in range(self._num_regions):
            for j in range(self._region_size[i]):
                key_array[idx] = self._region_keys[i][j]
                value_array[idx] = self._region_values[i][j]
                idx += 1

        self._length = 0
        self._sqrtval = 0
        self._num_regions = 0

        self._region_values.clear()
        self._region_size.clear()
        self._region_keys.clear()
        self._key_to_position.clear()

        self._initialize_container(key_array, value_array)

    cpdef int index(self, object key):
        """
        Return the integer position of a key in the container.

        INPUT:

        - ``key`` -- key stored in the container

        OUTPUT:

        - integer index of the key in the ordered container

        EXAMPLES::

            sage: keys = [(0,0),(1,1),(2,2)]
            sage: vals = [10,20,30]
            sage: c = Arbitrary_Indexed_Container_Sqrt(keys, vals)
            sage: c.index((1,1))
            1

        TESTS::

            sage: c = Arbitrary_Indexed_Container_Sqrt([(0,0)], [5])
            sage: c.index((1,1))
            Traceback (most recent call last):
            ...
            LookupError
        """

        if key not in self._key_to_position:
            raise LookupError(f"{key} not found")

        region_idx, offset = self._key_to_position[key]

        cdef Py_ssize_t i, idx = 0

        for i in range(region_idx):
            idx += self._region_size[i]

        return idx + offset

    def __getitem__(self, object key):
        """
        Retrieve a value or a slice of values.

        If ``key`` is a stored key, the corresponding value is returned.
        If ``key`` is a slice with keys as ``start`` and ``stop``,
        a list of values in that range is returned.

        EXAMPLES::

            sage: keys = [(0,0),(1,1),(2,2)]
            sage: vals = [5,6,7]
            sage: c = Arbitrary_Indexed_Container_Sqrt(keys, vals)
            sage: c[(1,1)]
            6
            sage: c[(0,0):(2,2)]
            [5, 6]

        TESTS::

            sage: c = Arbitrary_Indexed_Container_Sqrt([(0,0)], [1])
            sage: c[(1,1)]
            Traceback (most recent call last):
            ...
            LookupError
        """
        cdef Py_ssize_t start_region_idx, start_region_internal_idx, stop_region_idx, stop_region_internal_idx

        if not isinstance(key, slice):
            if key not in self._key_to_position:
                raise LookupError(f"{key} not found")

            region_idx, offset = self._key_to_position[key]
            return self._region_values[region_idx][offset]
        else:
            # Check if container should be rebuilt first
            if self._lazy_reinitialize_flag:
                self.reinitialize_container()

            result = []

            start_region_idx, start_region_internal_idx = self._key_to_position[key.start]
            stop_region_idx, stop_region_internal_idx = self._key_to_position[key.stop]

            # If start is exactly equal to stop (faster to compare indices than keys for some large datatypes)
            # then return the empty list
            if start_region_idx == stop_region_idx and start_region_internal_idx == stop_region_internal_idx:
                return result

            if key.step is not None and key.step == 0:
                raise ValueError("Slice step cannot be zero")

            # Ensure stop position is after start position
            if ((start_region_idx > stop_region_idx) or
                (start_region_idx == stop_region_idx and start_region_internal_idx > stop_region_internal_idx)):
                raise LookupError(
                                    f"Start position {(start_region_idx, start_region_internal_idx)}"
                                    f"represented by {key.start} is equal or later than end position"
                                    f"{(start_region_idx, start_region_internal_idx)} reprsented"
                                    f"by {key.stop}")
            
            # Add regions together to get resulting range query
            # Two cases:
            #   1. start region < stop region
            #   2. start region == stop region
            if start_region_idx < stop_region_idx:
                # Add values in start region
                result.extend(self._region_values[start_region_idx][start_region_internal_idx:])
                
                # Add values in regions between start and stop regions
                for i in range(start_region_idx + 1, stop_region_idx):
                    result.extend(self._region_values[i])
                
                # Add values in stop region
                result.extend(self._region_values[stop_region_idx][:stop_region_internal_idx])

            else:
                # start and stop region are same, get subarray of that region
                result = self._region_values[start_region_idx][start_region_internal_idx:stop_region_internal_idx]
            
            # Return stepped result
            if key.step is None or key.step == 1:
                return result
            else:
                return result[::key.step]

    def __setitem__(self, object key, object value):
        """
        Set the value associated with a key.

        If the key already exists its value is updated.
        Otherwise the key/value pair is appended to the container.

        EXAMPLES::

            sage: c = Arbitrary_Indexed_Container_Sqrt([(0,0)], [1])
            sage: c[(0,0)] = 5
            sage: c[(0,0)]
            5

        TESTS::

            sage: c = Arbitrary_Indexed_Container_Sqrt()
            sage: c[(1,1)] = 7
            sage: c[(1,1)]
            7
        """

        # If key already exists, set the item
        if key in self._key_to_position:

            region_idx, offset = self._key_to_position[key]
            self._region_values[region_idx][offset] = value
            return

        # Otherwise, append it to the end
        self.append(key, value)

    cpdef object get_contiguous_array(self):
        """
        Return a flattened representation of the container.

        OUTPUT:

        - tuple ``(array, mapping)``
        - ``array`` -- list of values in order
        - ``mapping`` -- dictionary mapping keys to indices

        EXAMPLES::

            sage: c = Arbitrary_Indexed_Container_Sqrt([(0,0),(1,1)], [10,20])
            sage: arr, mp = c.get_contiguous_array()
            sage: arr
            [10, 20]
            sage: mp[(1,1)]
            1

        TESTS::

            sage: c = Arbitrary_Indexed_Container_Sqrt()
            sage: arr, mp = c.get_contiguous_array()
            sage: arr
            []
        """

        cdef Py_ssize_t idx = 0
        cdef Py_ssize_t i, j

        array = [None] * self._length
        mapping = {}

        for i in range(self._num_regions):
            for j in range(self._region_size[i]):
                key = self._region_keys[i][j]
                array[idx] = self._region_values[i][j]
                mapping[key] = idx
                idx += 1

        return array, mapping

    cpdef void append(self, object key, object value):
        """
        Append a key/value pair to the end of the container.

        INPUT:

        - ``key`` -- unique key
        - ``value`` -- associated value

        EXAMPLES::

            sage: c = Arbitrary_Indexed_Container_Sqrt()
            sage: c.append((0,0), 10)
            sage: c[(0,0)]
            10

        TESTS::

            sage: c = Arbitrary_Indexed_Container_Sqrt([(0,0)], [1])
            sage: c.append((0,0), 2)
            Traceback (most recent call last):
            ...
            AssertionError
        """

        if key in self._key_to_position:
            raise AssertionError(f"{key} already exists")

        if len(self._region_values) == 0:

            self._region_values.append([])
            self._region_keys.append([])
            self._region_size.append(0)
            self._num_regions += 1

        region_idx = self._num_regions - 1
        offset = self._region_size[region_idx]

        self._region_values[region_idx].append(value)
        self._region_keys[region_idx].append(key)

        self._key_to_position[key] = (region_idx, offset)

        self._region_size[region_idx] += 1
        self._length += 1

        # Check if data structure is inefficient
        if self._region_size[region_idx] > self._reinitialize_ratio * self._sqrtval:
            # Indicate that reinialization is required now
            # that the data structure has become inefficient
            self._lazy_reinitialize_flag = True
    
    cpdef void extend(self, list key_array, list value_array):
        """
        Append multiple key/value pairs to the container.

        The elements are appended directly to the last region.
        If the region becomes too large the container will be
        lazily reinitialized.

        INPUT:

        - ``key_array`` -- list of keys
        - ``value_array`` -- list of values

        EXAMPLES::

            sage: c = Arbitrary_Indexed_Container_Sqrt()
            sage: c.extend([(0,0),(1,1)], [10,20])
            sage: c[(1,1)]
            20

        TESTS::

            sage: c = Arbitrary_Indexed_Container_Sqrt([(0,0)], [1])
            sage: c.extend([(0,0)], [5])
            Traceback (most recent call last):
            ...
            AssertionError
        """

        cdef Py_ssize_t n, i, region_idx, offset
        n = len(key_array)
        if n != len(value_array):
            raise AssertionError("Lengths of key_array and value_array must match")

        # Check for duplicates before appending
        for i in range(n):
            if key_array[i] in self._key_to_position:
                raise AssertionError(f"{key_array[i]} already exists")

        # Make sure at least one region exists
        if len(self._region_values) == 0:
            self._region_values.append([])
            self._region_keys.append([])
            self._region_size.append(0)
            self._num_regions += 1

        # Get last region index and the offset
        region_idx = max(0, self._num_regions - 1)
        offset = self._region_size[region_idx]

        # Extend the last region directly
        self._region_values[region_idx].extend(value_array)
        self._region_keys[region_idx].extend(key_array)
        self._region_size[region_idx] += n
        self._length += n

        # Add _key_to_position mapping for new keys
        for i in range(n):
            self._key_to_position[key_array[i]] = (region_idx, offset)
            offset += 1

        # Mark lazy reinitialization if region grew too large
        if self._region_size[region_idx] > self._reinitialize_ratio * self._sqrtval:
            self._lazy_reinitialize_flag = True

    cpdef void insert(self, Py_ssize_t idx, object key, object value):
        """
        Insert a key/value pair at a given index.

        The container locates the region containing the index
        and inserts the element there.

        INPUT:

        - ``idx`` -- integer position
        - ``key`` -- unique key
        - ``value`` -- associated value

        EXAMPLES::

            sage: c = Arbitrary_Indexed_Container_Sqrt([(0,0),(2,2)], [1,3])
            sage: c.insert(1,(1,1),2)
            sage: c[(1,1)]
            2
            sage: c.index((1,1))
            1

        TESTS::

            sage: c = Arbitrary_Indexed_Container_Sqrt([(0,0)], [1])
            sage: c.insert(5,(1,1),2)
            Traceback (most recent call last):
            ...
            AssertionError
        """

        cdef Py_ssize_t region_idx = 0
        cdef Py_ssize_t region_start = 0
        cdef Py_ssize_t insert_offset

        # Check if container should be rebuilt first
        if self._lazy_reinitialize_flag:
            self.reinitialize_container()

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

        self._region_values[region_idx].insert(insert_offset, value)
        self._region_keys[region_idx].insert(insert_offset, key)

        self._region_size[region_idx] += 1
        self._length += 1

        # update offsets inside region
        cdef Py_ssize_t j
        for j in range(insert_offset, self._region_size[region_idx]):

            k = self._region_keys[region_idx][j]
            self._key_to_position[k] = (region_idx, j)      # Update offset to new index

        # Check if data structure is inefficient
        if self._region_size[region_idx] > self._reinitialize_ratio * self._sqrtval:
            # Indicate that reinialization is required now
            # that the data structure has become inefficient
            self._lazy_reinitialize_flag = True