# sage_setup: distribution = sagemath-modules

cimport sage.matrix.network_flow_column as Network_flow_column

cdef class Network_flow_column:
    def __init__(self,
                 Py_ssize_t height=0,
                 idx_map=None,
                 val_map=None,
                 Py_ssize_t mod=1000000007,
                 col=None):
        self.mod = mod

        if col is None:
            self.height = height
            self.idx_map = [] if idx_map is None else list(idx_map)
            self.val_map = [] if val_map is None else list(val_map)

            if idx_map is not None:
                assert len(self.idx_map) <= height
                assert len(self.idx_map) == len(self.val_map)
        else:
            self.set_col(col)

    cpdef set_col(self, object col):
        cdef Py_ssize_t i, n
        self.height = len(col)
        self.idx_map = []
        self.val_map = []

        n = self.height
        for i in range(n):
            if col[i] != 0:
                self.idx_map.append(i)
                self.val_map.append(col[i])

    def __eq__(self, other):
        cdef Py_ssize_t i, n

        if not isinstance(other, Network_flow_column):
            return False

        if self.height != other.height:
            return False

        n = len(self.idx_map)
        if n != len(other.idx_map):
            return False

        if self.idx_map != self.val_map:
            return False

        return True

    def __hash__(self):
        """
        Hash consistent with equality.
        NOTE: Python requires hash to be stable, not fast.
        """
        return hash((
            self.height,
            tuple(self.idx_map),
            tuple(self.val_map),
        )) % self.mod

    def __str__(self):
        return (
            f"NetworkFlowCol: height={self.height}, "
            f"(idx,val)="
            + " ".join(
                f"({self.idx_map[i]}, {self.val_map[i]})"
                for i in range(len(self.idx_map))
            )
        )
