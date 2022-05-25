from itertools import combinations, permutations

def goodman_kruskal_gamma(m, n):
    """ 
    compute the Goodman and Kruskal gamma rank correlation coefficient; 
    this statistic ignores ties is unsuitable when the number of ties in the
    data is high. it's also slow. 
    >>> x = [2, 8, 5, 4, 2, 6, 1, 4, 5, 7, 4]
    >>> y = [3, 9, 4, 3, 1, 7, 2, 5, 6, 8, 3]
    >>> goodman_kruskal_gamma(x, y)
    0.9166666666666666
    """
    print(m,n)
    num = 0
    den = 0
    for (i, j) in permutations(range(len(m)), 2):
        m_dir = m[i] - m[j]
        n_dir = n[i] - n[j]
        sign = m_dir * n_dir
        print(i,j,m_dir,n_dir,sign)
        if sign > 0:
            num += 1
            den += 1
        elif sign < 0:
            num -= 1
            den += 1
    return num / float(den)

x = [2, 8, 5, 4, 2, 6, 1, 4, 5, 7, 4]
y = [3, 9, 4, 3, 1, 7, 2, 5, 6, 8, 3]
print(goodman_kruskal_gamma(x, y))
t1 = [1,2,3,4]
t2 = [1,2,5,6]
t3 = [5,6,7,8]
print(goodman_kruskal_gamma(t1, t2))
print(goodman_kruskal_gamma(t1, t3))
print(goodman_kruskal_gamma(t2, t3))


