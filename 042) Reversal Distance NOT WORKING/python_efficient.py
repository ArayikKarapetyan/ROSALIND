from collections import deque

N = 10
TARGET = tuple(range(1, N+1))
REVERSALS = [(i, j) for i in range(N) for j in range(i+1, N)]

def neighbors(perm):
    perm = list(perm)
    for i, j in REVERSALS:
        newp = perm[:]
        newp[i:j+1] = reversed(newp[i:j+1])
        yield tuple(newp)

def reversal_distance(perm):
    if perm == TARGET:
        return 0

    q1 = deque([(perm, 0)])
    q2 = deque([(TARGET, 0)])

    dist1 = {perm: 0}
    dist2 = {TARGET: 0}

    while q1 and q2:
        # expand smaller frontier
        if len(q1) <= len(q2):
            cur, d = q1.popleft()
            for nxt in neighbors(cur):
                if nxt in dist1:
                    continue
                if nxt in dist2:
                    return d + 1 + dist2[nxt]
                dist1[nxt] = d + 1
                q1.append((nxt, d + 1))
        else:
            cur, d = q2.popleft()
            for nxt in neighbors(cur):
                if nxt in dist2:
                    continue
                if nxt in dist1:
                    return d + 1 + dist1[nxt]
                dist2[nxt] = d + 1
                q2.append((nxt, d + 1))

    return -1

def normalize(p1, p2):
    pos = {v: i+1 for i, v in enumerate(p2)}
    return tuple(pos[x] for x in p1)

def main():
    with open("Dataset.txt") as f:
        lines = [line.strip() for line in f if line.strip()]

    results = []
    for i in range(0, len(lines), 2):
        p1 = list(map(int, lines[i].split()))
        p2 = list(map(int, lines[i+1].split()))
        perm = normalize(p1, p2)
        results.append(str(reversal_distance(perm)))

    print(" ".join(results))

if __name__ == "__main__":
    main()
: