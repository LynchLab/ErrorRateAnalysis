def list_ep(als,step=1):
    """
    Return a nested list of ends of consecutive runs given an ordered list of numbers.

    Arguments:
    - als   - an ordered list of numbers
    - step  - step size, to define "consecutive runs" [ 1 ]
    """
    # first check if the input list is sorted in ascending order
    try:
        assert als == sorted(als)
    except AssertionError:
        print("input list needs to be in ascending order")
        return
        
    # initialize output list
    out = []
    # no. of elements
    n = len(als)
    
    # if there are no elements, terminate
    if n == 0:
        return []
    
    # last no. stored, initialize with one number less than the first
    lns = als[0]-1
    # for every number
    for x in range(n):
        # first element
        e1 = als[x]
        # if it is not greater than the last stored no., 
        # then nothing more to do for this no.
        if e1 <= lns:
            continue
        # initial range end-points for this no.
        r1 = e1
        r2 = e1
        # initialize expected no.
        expn = e1
        # compare this no. with all subsequent no.s
        for y in range(x+1,n):
            # second element
            e2 = als[y]
            # expected no. increment by step size in every loop
            expn += step
            # if the next no. matches the expectation
            if e2 == expn:
                # update right end-point of the range with the current no.
                r2 = e2
            # when you reach a no. that doesn't
            else:
                # update the same end-point with the previous no.
                r2 = als[y-1]
                break
        # store current end-points
        out.append((r1,r2))
        # update the last number stored
        lns = out[-1][1]
    return out
