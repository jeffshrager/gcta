def _weightedrankedtau(ordered[:] x, ordered[:] y, intp_t[:] rank, weigher, bool additive):
    cdef intp_t i, first
    cdef float64_t t, u, v, w, s, sq
    cdef int64_t n = np.int64(len(x))
    cdef float64_t[::1] exchanges_weight = np.zeros(1, dtype=np.float64)
    # initial sort on values of x and, if tied, on values of y
    cdef intp_t[::1] perm = np.lexsort((y, x))
    cdef intp_t[::1] temp = np.empty(n, dtype=np.intp) # support structure

    if weigher is None:
        weigher = lambda x: 1./(1 + x)

    if rank is None:
        # To generate a rank array, we must first reverse the permutation
        # (to get higher ranks first) and then invert it.
        rank = np.empty(n, dtype=np.intp)
        rank[...] = perm[::-1]
        _invert_in_place(rank)

    # weigh joint ties
    first = 0
    t = 0
    w = weigher(rank[perm[first]])
    s = w
    sq = w * w

    for i in range(1, n):
        if x[perm[first]] != x[perm[i]] or y[perm[first]] != y[perm[i]]:
            t += s * (i - first - 1) if additive else (s * s - sq) / 2
            first = i
            s = sq = 0

        w = weigher(rank[perm[i]])
        s += w
        sq += w * w

    t += s * (n - first - 1) if additive else (s * s - sq) / 2

    # weigh ties in x
    first = 0
    u = 0
    w = weigher(rank[perm[first]])
    s = w
    sq = w * w

    for i in range(1, n):
        if x[perm[first]] != x[perm[i]]:
            u += s * (i - first - 1) if additive else (s * s - sq) / 2
            first = i
            s = sq = 0

        w = weigher(rank[perm[i]])
        s += w
        sq += w * w

    u += s * (n - first - 1) if additive else (s * s - sq) / 2
    if first == 0: # x is constant (all ties)
        return np.nan

    # this closure recursively sorts sections of perm[] by comparing
    # elements of y[perm[]] using temp[] as support

    def weigh(intp_t offset, intp_t length):
        cdef intp_t length0, length1, middle, i, j, k
        cdef float64_t weight, residual

        if length == 1:
            return weigher(rank[perm[offset]])
        length0 = length // 2
        length1 = length - length0
        middle = offset + length0
        residual = weigh(offset, length0)
        weight = weigh(middle, length1) + residual
        if y[perm[middle - 1]] < y[perm[middle]]:
            return weight

        # merging
        i = j = k = 0

        while j < length0 and k < length1:
            if y[perm[offset + j]] <= y[perm[middle + k]]:
                temp[i] = perm[offset + j]
                residual -= weigher(rank[temp[i]])
                j += 1
            else:
                temp[i] = perm[middle + k]
                exchanges_weight[0] += weigher(rank[temp[i]]) * (
                    length0 - j) + residual if additive else weigher(
                    rank[temp[i]]) * residual
                k += 1
            i += 1

        perm[offset+i:offset+i+length0-j] = perm[offset+j:offset+length0]
        perm[offset:offset+i] = temp[0:i]
        return weight

    # weigh discordances
    weigh(0, n)

    # weigh ties in y
    first = 0
    v = 0
    w = weigher(rank[perm[first]])
    s = w
    sq = w * w

    for i in range(1, n):
        if y[perm[first]] != y[perm[i]]:
            v += s * (i - first - 1) if additive else (s * s - sq) / 2
            first = i
            s = sq = 0

        w = weigher(rank[perm[i]])
        s += w
        sq += w * w

    v += s * (n - first - 1) if additive else (s * s - sq) / 2
    if first == 0: # y is constant (all ties)
        return np.nan

    # weigh all pairs
    s = sq = 0
    for i in range(n):
        w = weigher(rank[perm[i]])
        s += w
        sq += w * w

    tot = s * (n - 1) if additive else (s * s - sq) / 2

    tau = ((tot - (v + u - t)) - 2. * exchanges_weight[0]
           ) / np.sqrt(tot - u) / np.sqrt(tot - v)
    return min(1., max(-1., tau))
