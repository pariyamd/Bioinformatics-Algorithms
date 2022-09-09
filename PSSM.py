import numpy as np
import math
from itertools import chain, combinations

N = int(input())
seqs = []
for i in range(N):
    seqs.append(list(input()))
db = input()
n_seq = len(seqs[0])

chars = set(np.array(seqs).flatten())
chars.add("-")

freq = {ch: 0 for ch in chars}

seqs_n = {ch: [0] * n_seq for ch in chars}
for i in range(N):
    for j, ch in enumerate(seqs[i]):
        freq[ch] += 1
        seqs_n[ch][j] += 1

for i in range(n_seq):
    for ch in chars:
        seqs_n[ch][i] = (seqs_n[ch][i] + 2) / (N + (2 * len(chars)))

scores = {ch: [0] * n_seq for ch in chars}

score_mx = {i: 0 for i in range(n_seq)}
maxs = {i: "" for i in range(n_seq)}

for ch in chars:
    f = sum(seqs_n[ch])
    f /= n_seq
    for j in range(len(seqs_n[ch])):
        scores[ch][j] = round(math.log(seqs_n[ch][j] / f, 2), 3)
        if scores[ch][j] > score_mx[j]:
            score_mx[j] = scores[ch][j]
            maxs[j] = ch


def calc_score(db_slice):
    score = 0
    for i in range(len(db_slice)):
        score += scores[db_slice[i]][i]
    return score


def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


psets = powerset(range(n_seq))
dict_pset = {i: [] for i in range(n_seq + 1)}
for pset in psets:
    dict_pset[n_seq - len(pset)].append(pset)


def create_seqs(seq: str):
    seqs = []
    for pset in dict_pset[len(seq)]:
        seq_c = list(seq)
        for pos in pset:
            seq_c.insert(pos, "-")
        seqs.append("".join(seq_c))
    return seqs


max_score = -10000
max_slice = ""
# print("profile", profile)
for s in range(len(db) - n_seq):
    for i in range(1, n_seq):
        seq = db[s:s + i]
        seqs = create_seqs(seq)
        for seq in seqs:
            score = calc_score(seq)
            if score > max_score:
                max_score = score
                max_slice = seq

print(max_slice)
