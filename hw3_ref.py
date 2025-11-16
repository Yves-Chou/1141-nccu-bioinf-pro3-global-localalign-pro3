def alignment(input_path, score_path, output_path, aln, gap):
    # Read first two sequences from FASTA
    names, seqs, cur = [], [], []
    with open(input_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
                names.append(line[1:].strip())
            else:
                cur.append(line)
    if cur:
        seqs.append("".join(cur))
    if len(seqs) < 2:
        return

    s1, s2 = seqs[0], seqs[1]
    n1 = names[0]
    n2 = names[1] if len(names) > 1 else "seq2"

    # Read PAM score matrix
    lines = [l.strip() for l in open(score_path)
             if l.strip() and not l.lstrip().startswith("#")]
    header = lines[0].split()
    score = {a: {} for a in header}
    for row in lines[1:]:
        p = row.split()
        r = p[0]
        for j, v in enumerate(p[1:]):
            score[r][header[j]] = int(v)

    def sub(a, b):
        return score[a][b]

    mode = aln.lower()

    # GLOBAL: Needleman–Wunsch
    if mode == "global":
        m, n = len(s1), len(s2)
        S = [[0] * (n + 1) for _ in range(m + 1)]   # score matrix
        T = [[None] * (n + 1) for _ in range(m + 1)]  # traceback matrix

        # initialize first column (gaps in seq2)
        for i in range(1, m + 1):
            S[i][0] = S[i - 1][0] + gap
            T[i][0] = "up"

        # initialize first row (gaps in seq1)
        for j in range(1, n + 1):
            S[0][j] = S[0][j - 1] + gap
            T[0][j] = "left"

        # fill DP table, tie-breaking: diag >= up >= left
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                d = S[i - 1][j - 1] + sub(s1[i - 1], s2[j - 1])
                u = S[i - 1][j] + gap
                l = S[i][j - 1] + gap
                if d >= u and d >= l:
                    S[i][j] = d
                    T[i][j] = "diag"
                elif u >= l:
                    S[i][j] = u
                    T[i][j] = "up"
                else:
                    S[i][j] = l
                    T[i][j] = "left"

        # traceback from bottom-right
        i, j = len(s1), len(s2)
        a1, a2 = [], []
        while i > 0 or j > 0:
            t = T[i][j]
            if t == "diag":
                a1.append(s1[i - 1])
                a2.append(s2[j - 1])
                i -= 1
                j -= 1
            elif t == "up":
                a1.append(s1[i - 1])
                a2.append("-")
                i -= 1
            else:  # "left"
                a1.append("-")
                a2.append(s2[j - 1])
                j -= 1

        aln1 = "".join(reversed(a1))
        aln2 = "".join(reversed(a2))

        # write FASTA output
        with open(output_path, "w") as f:
            f.write(f">{n1}\n{aln1}\n>{n2}\n{aln2}\n")
        return

    # LOCAL: ungapped best substring
    if mode == "local":
        m, n = len(s1), len(s2)
        best = -10**18          # current best score
        wins = []               # (i, j, length) windows with best score

        # brute-force all ungapped alignments (continuous substrings)
        for i in range(m):
            for j in range(n):
                sc = 0
                k = 0
                while i + k < m and j + k < n:
                    sc += sub(s1[i + k], s2[j + k])
                    if sc > best:
                        best = sc
                        wins = [(i, j, k + 1)]
                    elif sc == best:
                        wins.append((i, j, k + 1))
                    k += 1

        if not wins:
            return

        # keep only windows with maximum length
        Lmax = max(w[2] for w in wins)
        cand = [w for w in wins if w[2] == Lmax]

        # sort by aligned strings: first protein1 segment, then protein2 segment
        cand.sort(key=lambda w: (s1[w[0]:w[0] + w[2]],
                                 s2[w[1]:w[1] + w[2]]))

        # write all local alignments in FASTA format
        with open(output_path, "w") as f:
            for i, j, L in cand:
                aln1 = s1[i:i + L]
                aln2 = s2[j:j + L]
                f.write(f">{n1}\n{aln1}\n>{n2}\n{aln2}\n")
        return
