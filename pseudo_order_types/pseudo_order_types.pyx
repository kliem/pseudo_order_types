from itertools import combinations, permutations, product, chain
from memory_allocator.memory_allocator cimport MemoryAllocator
from libc.string cimport memset, memcpy
from libc.stdio cimport FILE, fopen, fclose, fwrite, fread
cimport python_exc


cdef struct impli:
    # Read like this: if all pairs in ``same`` are the same
    # and all pairs in ``opposite`` are opposite,
    # then ``max_index`` must be likewise/opposite to ``other`` according to ``typ``.
    int typ  # 1 for likewise, -1 for opposite.
    int max_index
    int other
    int n_same
    int n_opposite
    int* same
    int* opposite


cdef struct order_type_iter:
    int end
    impli** implications
    int* n_implications
    bint* skips
    int* new_order
    char** choices
    int** lambda_matrix
    int* lambda_matrix_entries
    int pos
    int minimal_pos
    int n


cdef order_type_iter* prepare_order_type_iter(int n, MemoryAllocator mem) except NULL:
    if n < 3:
        raise ValueError("at least 3 points needed")
    if n == 3:
        raise ValueError("the special case of 3 points needs to be treated seperately")

    cdef int i, j, k
    combs = tuple(combinations(range(n), 3))
    combs_inv = {comb: i for i, comb in enumerate(combs)}
    cdef order_type_iter* iterator = <order_type_iter*> mem.allocarray(1, sizeof(order_type_iter))
    cdef int end = len(combs)

    # Only needed locally.
    cdef MemoryAllocator mem2 = MemoryAllocator()
    cdef int *size_implications = <int*> mem2.calloc(end, sizeof(int))

    cdef impli** implications = <impli**> mem.calloc(end, sizeof(impli*))
    cdef int *n_implications = <int*> mem.calloc(end, sizeof(int))
    cdef int new_size
    cdef impli* current_implication

    for same, opposite in all_implications(n):
        maxi = ()
        prev_max = ()
        for x, y in chain(same, opposite):
            if y > maxi:
                if maxi > prev_max:
                    prev_max = maxi
                maxi = y
            if maxi > y > prev_max:
                prev_max = y
            if maxi > x > prev_max:
                prev_max = x
        max_index = combs_inv[maxi]
        prevmax_index = combs_inv[prev_max]

        if n_implications[prevmax_index] == size_implications[prevmax_index]:
            new_size = 2 * size_implications[prevmax_index]
            if new_size == 0:
               new_size = 1

            implications[prevmax_index] = <impli*> mem.reallocarray(implications[prevmax_index], new_size, sizeof(impli))
            size_implications[prevmax_index] = new_size

        current_implication = &implications[prevmax_index][n_implications[prevmax_index]]
        n_implications[prevmax_index] += 1

        current_implication.same = <int*> mem.allocarray(2*len(same), sizeof(int))
        current_implication.opposite = <int*> mem.allocarray(2*len(opposite), sizeof(int))
        current_implication.n_same = 0
        current_implication.n_opposite = 0
        current_implication.max_index = max_index

        for i, (x, y) in enumerate(same):
            if y == maxi:
                current_implication.other = combs_inv[x]
                current_implication.typ = -1
            else:
                current_implication.same[2*current_implication.n_same] = combs_inv[x]
                current_implication.same[2*current_implication.n_same + 1] = combs_inv[y]
                current_implication.n_same += 1

        for i, (x, y) in enumerate(opposite):
            if y == maxi:
                current_implication.other = combs_inv[x]
                current_implication.typ = 1
            else:
                current_implication.opposite[2*current_implication.n_opposite] = combs_inv[x]
                current_implication.opposite[2*current_implication.n_opposite + 1] = combs_inv[y]
                current_implication.n_opposite += 1

    # We may assume that ``0`` is a vertex of the convex hull and
    # ``1, ..., n-1`` are sorted counter-clockwise.
    # This is exactly where we assume that the chirotope is acyclic.
    cdef int minimal_pos = combs_inv[1, 2, 3]
    cdef int pos = minimal_pos
    cdef char** choices = <char**> mem.allocarray(end, sizeof(char*))
    for i in range(end):
        choices[i] = <char*> mem.allocarray(end, sizeof(char))
    for i in range(minimal_pos):
        choices[pos][i] = -1
        choices[pos-1][i] = -1
    for i in range(minimal_pos, end):
        choices[pos][i] = 0
        choices[pos-1][i] = 0

    # Mark the positions, for which the sign was set by implication.
    cdef bint *skips = <bint*> mem.calloc(end, sizeof(bint))

    for i in range(end):
        skips[i] = False

    cdef int **lambda_matrix = <int**> mem.allocarray(n, sizeof(int*))
    for i in range(n):
        lambda_matrix[i] = <int*> mem.allocarray(n, sizeof(int))

    # For each combinations we save which entries of the lambda matrix should be set, if the entry is positive.
    cdef int *lambda_matrix_entries = <int*> mem.allocarray(3*end, sizeof(int))
    for i, (a, b, c) in enumerate(combs):
        lambda_matrix_entries[i*3] = a
        lambda_matrix_entries[i*3+1] = b
        lambda_matrix_entries[i*3+2] = c

    # We will only give those chirotopes with minimal lambda matrix with respect translation
    # and permutation rows and columns simultanously.
    # Again, for any other choice of a first vertex, we may assume the other points to be sorted
    # counter-clockwise.
    cdef int* new_order = <int*> mem.allocarray(n, sizeof(int))

    iterator.end = end
    iterator.implications = implications
    iterator.n_implications = n_implications
    iterator.skips = skips
    iterator.new_order = new_order
    iterator.choices = choices
    iterator.lambda_matrix = lambda_matrix
    iterator.lambda_matrix_entries = lambda_matrix_entries
    iterator.pos = pos
    iterator.minimal_pos = minimal_pos
    iterator.n = n

    return iterator


cdef bint lambda_matrix_is_minimal(order_type_iter* it):
    cdef int i, j, k
    cdef int n = it.n
    for i in range(n):
        memset(it.lambda_matrix[i], 0, n*sizeof(int))

    for i in range(it.end):
        if it.choices[it.pos][i] == 1:
            it.lambda_matrix[it.lambda_matrix_entries[i*3]][it.lambda_matrix_entries[i*3+1]] += 1
            it.lambda_matrix[it.lambda_matrix_entries[i*3+1]][it.lambda_matrix_entries[i*3+2]] += 1
            it.lambda_matrix[it.lambda_matrix_entries[i*3+2]][it.lambda_matrix_entries[i*3]] += 1
        else:
            it.lambda_matrix[it.lambda_matrix_entries[i*3+2]][it.lambda_matrix_entries[i*3+1]] += 1
            it.lambda_matrix[it.lambda_matrix_entries[i*3]][it.lambda_matrix_entries[i*3+2]] += 1
            it.lambda_matrix[it.lambda_matrix_entries[i*3+1]][it.lambda_matrix_entries[i*3]] += 1

    # See, if the lambda matrix is minimal.
    for i in range(1, n):
        it.new_order[n-1] = -1
        for j in range(n):
            if j == i:
                it.new_order[0] = i
            else:
                it.new_order[it.lambda_matrix[i][j] + 1] = j

        if it.new_order[n-1] == -1:
            # the corresponding point does not lie on the boundary of the convex hull
            continue

        for k in range(n):
            for j in range(n):
                foo = it.lambda_matrix[it.new_order[j]][it.new_order[k]] - it.lambda_matrix[j][k]
                if foo:
                    break
            else:
                continue
            break
        if foo < 0:
            # Turns out our matrix is not the smallest one.
            return False
    else:
        # Now the transposed.
        for i in range(n):
            it.new_order[n-1] = -1
            for j in range(n):
                if j == i:
                    it.new_order[0] = i
                else:
                    it.new_order[it.lambda_matrix[j][i] + 1] = j

            if it.new_order[n-1] == -1:
                # the corresponding point does not lie on the boundary of the convex hull
                continue

            for k in range(n):
                for j in range(n):
                    foo = it.lambda_matrix[it.new_order[k]][it.new_order[j]] - it.lambda_matrix[j][k]
                    if foo:
                        break
                else:
                    continue
                break
            if foo < 0:
                # Turns out our matrix is not the smallest one.
                return False

    return True


cdef bint next_order_type(order_type_iter* it):
    cdef int foo
    cdef int designated
    cdef impli* implication
    cdef int i, j
    cdef bint minimal

    while it.pos >= it.minimal_pos:
        if not it.skips[it.pos]:
            if it.choices[it.pos][it.pos] == 0:
                # memcpy was done below.
                it.choices[it.pos][it.pos] = 1
            elif it.choices[it.pos][it.pos] == 1:
                memcpy(it.choices[it.pos], it.choices[it.pos-1], it.end*sizeof(char))
                it.choices[it.pos][it.pos] = -1
            else:
                it.pos -= 1
                while it.skips[it.pos] and it.pos >= it.minimal_pos:
                    it.pos -= 1
                continue
        for i in range(it.n_implications[it.pos]):
            implication = &(it.implications[it.pos][i])

            if (all(it.choices[it.pos][implication.same[2*j]] == it.choices[it.pos][implication.same[2*j + 1]]
                    for j in range(implication.n_same)) and
                all(it.choices[it.pos][implication.opposite[2*j]] != it.choices[it.pos][implication.opposite[2*j + 1]]
                    for j in range(implication.n_opposite))):
                designated = it.choices[it.pos][implication.other]*implication.typ
                if it.choices[it.pos][implication.max_index] == 0:
                    it.choices[it.pos][implication.max_index] = designated
                elif it.choices[it.pos][implication.max_index] != designated:
                    # Not a valid choice.
                    # Only track back, if the last choice was forced.
                    # Otherwise just undo the last choice.
                    while it.skips[it.pos] and it.pos >= it.minimal_pos:
                        it.pos -= 1
                    break
        else:
            if it.pos == it.end - 1:
                # Valid choice so far.

                minimal = lambda_matrix_is_minimal(it)

                # Only track back, if the last choice was forced.
                while it.skips[it.pos] and it.pos >= it.minimal_pos:
                    it.pos -= 1

                if minimal:
                    return True
            else:
                it.pos += 1
                memcpy(it.choices[it.pos], it.choices[it.pos-1], it.end*sizeof(char))
                it.skips[it.pos] = it.choices[it.pos][it.pos] != 0

    return False


def pseudo_order_type_writer(int n, path):
    """
    Write all non-degenerate acyclic chirotopes with base set 0,..,n-1 of rank 3 to a file.

    The following symmetries are quotiened out:
    - relabeling of the points,
    - replacing ``chi`` by ``-chi``; switching all signs.

    INPUT:

    - ``n`` -- integer greater or equal to ``3``
    - ``path`` -- string

    .. SEEALSO::

        :func:`pseudo_order_type_iterator`

    EXAMPLES::

        >>> from pseudo_order_types.pseudo_order_types import pseudo_order_type_writer
        >>> from pseudo_order_types.pseudo_order_types import pseudo_order_type_iterator
        >>> it = pseudo_order_type_iterator(7)
        >>> pseudo_order_type_writer(7, 'pseudo_order_types_tmp_file')
        >>> it2 = pseudo_order_type_iterator(7, 'pseudo_order_types_tmp_file')
        >>> all(next(it) == a for a in it2)
        True
        >>> import os
        >>> os.remove("pseudo_order_types_tmp_file")
    """
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef order_type_iter* it
    cdef char trivial = -1
    if n != 3:
        it = prepare_order_type_iter(n, mem)

    cdef FILE* fp
    fp = fopen(path.encode('utf-8'), "w")
    if (fp==NULL):
        raise IOError("cannot open file {}".format(path))

    if n == 3:
        fwrite(&trivial, sizeof(char), 1, fp)
        fclose(fp)
        return

    try:
        while next_order_type(it):
            if it.end != fwrite(it.choices[it.end - 1], sizeof(char), it.end, fp):
                raise IOError("could not write to file")
            if python_exc.PyErr_CheckSignals():
                return
    finally:
        fclose(fp)


def _pseudo_order_type_iterator(int n):
    """
    Helper function to iterate of non-degenerate acyclic chirotopes of rank 3.

    .. SEEALSO::

        :func:`pseudo_order_type_iterator`
    """
    if n == 3:
        yield (-1,)
        return

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef order_type_iter* it = prepare_order_type_iter(n, mem)
    cdef int i

    while next_order_type(it):
        yield tuple(it.choices[it.end - 1][i] for i in range(it.end))


def pseudo_order_type_iterator(int n, path=None):
    """
    Iterate over all non-degenerate acyclic chirotopes with base set 0,..,n-1 of rank 3.

    The following symmetries are quotiened out:
    - relabeling of the points,
    - replacing ``chi`` by ``-chi``; switching all signs.

    INPUT:

    - ``n`` -- integer greater or equal to ``3``
    - ``path`` -- string (default: ``None``); specify a file to which
      ``pseudo_order_type_writer`` as written the output

    OUTPUT:

    Each chirotope as tuple of signs ``-1`` or ``1``,
    one sign for each element in ``itertools.combinations(range(n), 3)``.

    EXAMPLES::

        >>> from pseudo_order_types.pseudo_order_types import pseudo_order_type_iterator
        >>> it = pseudo_order_type_iterator(5)
        >>> for i in it:
        ...     print(i)
        (-1, -1, -1, -1, -1, -1, 1, 1, 1, 1)
        (-1, -1, -1, -1, -1, -1, 1, 1, -1, -1)
        (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1)
        >>> it = pseudo_order_type_iterator(6)
        >>> next(it)
        (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
        >>> next(it)
        (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1)
        >>> next(it)
        (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1)
        >>> sum(1 for _ in pseudo_order_type_iterator(6))
        16
        >>> sum(1 for _ in pseudo_order_type_iterator(7))
        135
        >>> sum(1 for _ in pseudo_order_type_iterator(8))
        3315

    Store to a file with ``pseudo_order_type_writer``::

        >>> from pseudo_order_types.pseudo_order_types import pseudo_order_type_writer
        >>> from pseudo_order_types.pseudo_order_types import pseudo_order_type_iterator
        >>> pseudo_order_type_writer(6, 'pseudo_order_types_tmp_file')
        >>> it = pseudo_order_type_iterator(6, 'pseudo_order_types_tmp_file')
        >>> next(it)
        (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
        >>> import os
        >>> os.remove("pseudo_order_types_tmp_file")
    """
    if not path:
        yield from _pseudo_order_type_iterator(n)
        return
    cdef FILE* fp
    fp = fopen(path.encode('utf-8'), "r")
    if (fp == NULL):
        raise IOError("cannot open file {}".format(path))
    combs = tuple(combinations(range(n), 3))
    cdef size_t end = len(combs)
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef char* choices = <char*> mem.allocarray(len(combs), sizeof(char))
    cdef size_t count
    try:
        while True:
            count = fread(choices, sizeof(char), end, fp)
            if count == end:
                yield tuple(choices[i] for i in range(len(combs)))
            elif count == 0:
                return
            else:
                raise IOError("size of file is incorrect")
    finally:
        fclose(fp)


def all_implications(n, r=3):
    """
    Return all implications by the chirotope exchange axiom for rank ``r``.

    INPUT:

    - ``n`` -- integer greater or equal to ``3``; the size of the chirotope
    - ``r`` -- integer greater or equal to ``2``; the rank of the matroid

    OUTPUT:

    Each implication is given by as lists ``same, opposite``.

    Each element in ``same`` or ``opposite`` is given as tuples ``(A, B)``.
    The implication is of the form that it cannot hold that
    ``chi(A)*chi(B) >= 0`` for all ``(A, B)`` in ``same`` and
    ``chi(A)*chi(B) <= 0`` for all ``(A, B)`` in ``opposite``.

    So at least one of the inequalities cannot hold.

    EXAMPLES::

        >>> from pseudo_order_types.pseudo_order_types import all_implications, implication_printer
        >>> s = sorted(all_implications(5))
        >>> len(s)
        10
        >>> for x in s:
        ...     print(x)
        ((((0, 1, 2), (0, 3, 4)), ((0, 1, 4), (0, 2, 3))), (((0, 1, 3), (0, 2, 4)),))
        ((((0, 1, 2), (1, 3, 4)), ((0, 1, 4), (1, 2, 3))), (((0, 1, 3), (1, 2, 4)),))
        ((((0, 1, 2), (2, 3, 4)), ((0, 2, 4), (1, 2, 3))), (((0, 2, 3), (1, 2, 4)),))
        ((((0, 1, 3), (0, 2, 4)),), (((0, 1, 2), (0, 3, 4)), ((0, 1, 4), (0, 2, 3))))
        ((((0, 1, 3), (1, 2, 4)),), (((0, 1, 2), (1, 3, 4)), ((0, 1, 4), (1, 2, 3))))
        ((((0, 1, 3), (2, 3, 4)), ((0, 3, 4), (1, 2, 3))), (((0, 2, 3), (1, 3, 4)),))
        ((((0, 1, 4), (2, 3, 4)), ((0, 3, 4), (1, 2, 4))), (((0, 2, 4), (1, 3, 4)),))
        ((((0, 2, 3), (1, 2, 4)),), (((0, 1, 2), (2, 3, 4)), ((0, 2, 4), (1, 2, 3))))
        ((((0, 2, 3), (1, 3, 4)),), (((0, 1, 3), (2, 3, 4)), ((0, 3, 4), (1, 2, 3))))
        ((((0, 2, 4), (1, 3, 4)),), (((0, 1, 4), (2, 3, 4)), ((0, 3, 4), (1, 2, 4))))
        >>> for same, opposite in s:
        ...     implication_printer(same, opposite)
        One of the following does not hold:
        chi(0, 1, 2)*chi(0, 3, 4) >= 0
        chi(0, 1, 4)*chi(0, 2, 3) >= 0
        chi(0, 1, 3)*chi(0, 2, 4) <= 0
        One of the following does not hold:
        chi(0, 1, 2)*chi(1, 3, 4) >= 0
        chi(0, 1, 4)*chi(1, 2, 3) >= 0
        chi(0, 1, 3)*chi(1, 2, 4) <= 0
        One of the following does not hold:
        chi(0, 1, 2)*chi(2, 3, 4) >= 0
        chi(0, 2, 4)*chi(1, 2, 3) >= 0
        chi(0, 2, 3)*chi(1, 2, 4) <= 0
        One of the following does not hold:
        chi(0, 1, 3)*chi(0, 2, 4) >= 0
        chi(0, 1, 2)*chi(0, 3, 4) <= 0
        chi(0, 1, 4)*chi(0, 2, 3) <= 0
        One of the following does not hold:
        chi(0, 1, 3)*chi(1, 2, 4) >= 0
        chi(0, 1, 2)*chi(1, 3, 4) <= 0
        chi(0, 1, 4)*chi(1, 2, 3) <= 0
        One of the following does not hold:
        chi(0, 1, 3)*chi(2, 3, 4) >= 0
        chi(0, 3, 4)*chi(1, 2, 3) >= 0
        chi(0, 2, 3)*chi(1, 3, 4) <= 0
        One of the following does not hold:
        chi(0, 1, 4)*chi(2, 3, 4) >= 0
        chi(0, 3, 4)*chi(1, 2, 4) >= 0
        chi(0, 2, 4)*chi(1, 3, 4) <= 0
        One of the following does not hold:
        chi(0, 2, 3)*chi(1, 2, 4) >= 0
        chi(0, 1, 2)*chi(2, 3, 4) <= 0
        chi(0, 2, 4)*chi(1, 2, 3) <= 0
        One of the following does not hold:
        chi(0, 2, 3)*chi(1, 3, 4) >= 0
        chi(0, 1, 3)*chi(2, 3, 4) <= 0
        chi(0, 3, 4)*chi(1, 2, 3) <= 0
        One of the following does not hold:
        chi(0, 2, 4)*chi(1, 3, 4) >= 0
        chi(0, 1, 4)*chi(2, 3, 4) <= 0
        chi(0, 3, 4)*chi(1, 2, 4) <= 0
        >>> s = sorted(all_implications(5, 2))
        >>> implication_printer(*s[0])
        One of the following does not hold:
        chi(0, 1)*chi(2, 3) >= 0
        chi(0, 3)*chi(1, 2) >= 0
        chi(0, 2)*chi(1, 3) <= 0
        >>> s = sorted(all_implications(6, 4))
        >>> implication_printer(*s[0])
        One of the following does not hold:
        chi(0, 1, 2, 3)*chi(0, 1, 4, 5) >= 0
        chi(0, 1, 2, 5)*chi(0, 1, 3, 4) >= 0
        chi(0, 1, 2, 4)*chi(0, 1, 3, 5) <= 0
    """
    return tuple(set(_all_implications_iter(n, r)))


def _all_implications_iter(n, r=3):
    """
    Yield all implications by the chirotope exchange axiom for rank ``r`` with repeats.

    Each implication is given by as lists ``same, opposite``.

    Each element in ``same`` or ``opposite`` is given as tuples ``(A, B)``.
    The implication is of the form that it cannot hold that
    ``chi(A)*chi(B) >= 0`` for all ``(A, B)`` in ``same`` and
    ``chi(A)*chi(B) <= 0`` for all ``(A, B)`` in ``opposite``.

    So at least one of the inequalities cannot hold.
    """
    def sign(*args):
        for a, b in combinations(args, 2):
            if a == b:
                return 0
        counter = 0
        for i, j in combinations(range(r), 2):
            if args[i] > args[j]:
                counter += 1
        if counter % 2:
            return -1
        else:
            return 1

    # Note that we only store chi(a, b, c) for a < b < c.
    # So chi(a, b, c) = sign(a, b, c) * chi(sort(a, b, c)).

    def sort(*args):
        return tuple(sorted(args))

    # For each A and B, we get the impliciations by the exchange axiom.
    for A1 in combinations(range(n), r):
        for A in permutations(A1):
            for B1 in combinations(range(n), r):
                for B in permutations(B1):

                    # We rule out tautologies.
                    if A[0] in B:
                        continue
                    if all(i in B for i in A[1:]):
                        continue

                    same = []
                    opposite = []

                    # The chirotope exchange axiom for A, B is that from all the following,
                    # at least one has to be false:
                    # chi(a_1,a_2,...,a_n) * chi(b_1,b_2,b_3,...,b_n) <= 0
                    # chi(b_1,a_2,...,a_n) * chi(a_1,b_2,b_3,...,b_n) >= 0
                    # chi(b_2,a_2,...,a_n) * chi(b_1,a_1,b_3,...,b_n) >= 0
                    # ...
                    # chi(b_n,a_2,...,a_n) * chi(b_1,b_2,...,b_(n-1),a_1) >= 0

                    if sign(*A) * sign(*B) == 1:
                        opposite.append(sort(sort(*A), sort(*B)))
                    else:
                        same.append(sort(sort(*A), sort(*B)))

                    for i in range(r):
                        A2 = B[i:i+1] + A[1:]
                        B2 = B[:i] + A[0:1] + B[i+1:]
                        sig = sign(*A2) * sign(*B2)
                        if sig == 0:
                            continue
                        if sig == 1:
                            same.append(sort(sort(*A2), sort(*B2)))
                        else:
                            opposite.append(sort(sort(*A2), sort(*B2)))

                    yield (sort(*same), sort(*opposite))
                    if python_exc.PyErr_CheckSignals():
                        return


def implication_printer(same, opposite):
    r"""
    Print ``all_implications`` in a meaningful way.

    EXAMPLES::

        >>> from pseudo_order_types.pseudo_order_types import all_implications, implication_printer
        >>> s = sorted(all_implications(6))
        >>> implication_printer(*s[5])
        One of the following does not hold:
        chi(0, 1, 2)*chi(1, 4, 5) >= 0
        chi(0, 1, 5)*chi(1, 2, 4) >= 0
        chi(0, 1, 4)*chi(1, 2, 5) <= 0
    """
    output = "One of the following does not hold:"
    for a, b in same:
        output += "\nchi{}*chi{} >= 0".format(a, b)
    for a, b in opposite:
        output += "\nchi{}*chi{} <= 0".format(a, b)
    print(output)
