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

        # Whether or not the container needs to be reinitialized
        # Reinitialization required when container imbalanced and
        # and on __getitem__(slice) or insert
        # append, _initialize_container, reinitialize_container, __len__, index, __getitem__(index)
        self._lazy_reinitialize_required = False 

        self._initialize_container(key_array=key_array, value_array=value_array)

    cdef void _initialize_container(self, list key_array, list value_array):

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
        
        self._lazy_reinitialize_required = False

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
        cdef Py_ssize_t start_region_idx, start_region_internal_idx, stop_region_idx, stop_region_internal_idx

        if not isinstance(key, slice):
            if key not in self._key_to_position:
                raise LookupError(f"{key} not found")

            region_idx, offset = self._key_to_position[key]
            return self._regions[region_idx][offset]
        else:
            # Check if container should be rebuilt first
            if self._lazy_reinitialize_required:
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
                result.extend(self._regions[start_region_idx][start_region_internal_idx:])
                
                # Add values in regions between start and stop regions
                for i in range(start_region_idx + 1, stop_region_idx):
                    result.extend(self._regions[i])
                
                # Add values in stop region
                result.extend(self._regions[stop_region_idx][:stop_region_internal_idx])

            else:
                # start and stop region are same, get subarray of that region
                result = self._regions[start_region_idx][start_region_internal_idx:stop_region_internal_idx]
            
            # Return stepped result
            if key.step is None or key.step == 1:
                return result
            else:
                return result[::key.step]

    def __setitem__(self, object key, object value):

        # If key already exists, set the item
        if key in self._key_to_position:

            region_idx, offset = self._key_to_position[key]
            self._regions[region_idx][offset] = value
            return

        # Otherwise, append it to the end
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

        # Check if data structure is inefficient
        if self._region_size[region_idx] > self._reinitialize_ratio * self._sqrtval:
            # Indicate that reinialization is required now
            # that the data structure has become inefficient
            self._lazy_reinitialize_required = True
    
    cpdef void extend(self, list key_array, list value_array):
        """
        Append multiple key/value pairs directly into the last region.
        Marks the container for lazy reinitialization if the region grows too large.

        INPUT:
            key_array   -- list of keys
            value_array -- list of corresponding values

        Raises AssertionError if any key already exists.
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
        if len(self._regions) == 0:
            self._regions.append([])
            self._regions_key_values.append([])
            self._region_size.append(0)
            self._num_regions += 1

        # Get last region index and the offset
        region_idx = max(0, self._num_regions - 1)
        offset = self._region_size[region_idx]

        # Extend the last region directly
        self._regions[region_idx].extend(value_array)
        self._regions_key_values[region_idx].extend(key_array)
        self._region_size[region_idx] += n
        self._length += n

        # Add _key_to_position mapping for new keys
        for i in range(n):
            self._key_to_position[key_array[i]] = (region_idx, offset)
            offset += 1

        # Mark lazy reinitialization if region grew too large
        if self._region_size[region_idx] > self._reinitialize_ratio * self._sqrtval:
            self._lazy_reinitialize_required = True

    cpdef void insert(self, Py_ssize_t idx, object key, object value):

        cdef Py_ssize_t region_idx = 0
        cdef Py_ssize_t region_start = 0
        cdef Py_ssize_t insert_offset

        # Check if container should be rebuilt first
        if self._lazy_reinitialize_required:
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

        self._regions[region_idx].insert(insert_offset, value)
        self._regions_key_values[region_idx].insert(insert_offset, key)

        self._region_size[region_idx] += 1
        self._length += 1

        # update offsets inside region
        cdef Py_ssize_t j
        for j in range(insert_offset, self._region_size[region_idx]):

            k = self._regions_key_values[region_idx][j]
            self._key_to_position[k] = (region_idx, j)      # Update offset to new index

        # Check if data structure is inefficient
        if self._region_size[region_idx] > self._reinitialize_ratio * self._sqrtval:
            # Indicate that reinialization is required now
            # that the data structure has become inefficient
            self._lazy_reinitialize_required = True