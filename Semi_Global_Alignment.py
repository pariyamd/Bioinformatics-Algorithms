import numpy as np

PAM250 = {
    'A': {'A': 2, 'C': -2, 'D': 0, 'E': 0, 'F': -3, 'G': 1, 'H': -1, 'I': -1, 'K': -1, 'L': -2, 'M': -1, 'N': 0, 'P': 1,
          'Q': 0, 'R': -2, 'S': 1, 'T': 1, 'V': 0, 'W': -6, 'Y': -3},
    'C': {'A': -2, 'C': 12, 'D': -5, 'E': -5, 'F': -4, 'G': -3, 'H': -3, 'I': -2, 'K': -5, 'L': -6, 'M': -5, 'N': -4,
          'P': -3, 'Q': -5, 'R': -4, 'S': 0, 'T': -2, 'V': -2, 'W': -8, 'Y': 0},
    'D': {'A': 0, 'C': -5, 'D': 4, 'E': 3, 'F': -6, 'G': 1, 'H': 1, 'I': -2, 'K': 0, 'L': -4, 'M': -3, 'N': 2, 'P': -1,
          'Q': 2, 'R': -1, 'S': 0, 'T': 0, 'V': -2, 'W': -7, 'Y': -4},
    'E': {'A': 0, 'C': -5, 'D': 3, 'E': 4, 'F': -5, 'G': 0, 'H': 1, 'I': -2, 'K': 0, 'L': -3, 'M': -2, 'N': 1, 'P': -1,
          'Q': 2, 'R': -1, 'S': 0, 'T': 0, 'V': -2, 'W': -7, 'Y': -4},
    'F': {'A': -3, 'C': -4, 'D': -6, 'E': -5, 'F': 9, 'G': -5, 'H': -2, 'I': 1, 'K': -5, 'L': 2, 'M': 0, 'N': -3,
          'P': -5, 'Q': -5, 'R': -4, 'S': -3, 'T': -3, 'V': -1, 'W': 0, 'Y': 7},
    'G': {'A': 1, 'C': -3, 'D': 1, 'E': 0, 'F': -5, 'G': 5, 'H': -2, 'I': -3, 'K': -2, 'L': -4, 'M': -3, 'N': 0, 'P': 0,
          'Q': -1, 'R': -3, 'S': 1, 'T': 0, 'V': -1, 'W': -7, 'Y': -5},
    'H': {'A': -1, 'C': -3, 'D': 1, 'E': 1, 'F': -2, 'G': -2, 'H': 6, 'I': -2, 'K': 0, 'L': -2, 'M': -2, 'N': 2, 'P': 0,
          'Q': 3, 'R': 2, 'S': -1, 'T': -1, 'V': -2, 'W': -3, 'Y': 0},
    'I': {'A': -1, 'C': -2, 'D': -2, 'E': -2, 'F': 1, 'G': -3, 'H': -2, 'I': 5, 'K': -2, 'L': 2, 'M': 2, 'N': -2,
          'P': -2, 'Q': -2, 'R': -2, 'S': -1, 'T': 0, 'V': 4, 'W': -5, 'Y': -1},
    'K': {'A': -1, 'C': -5, 'D': 0, 'E': 0, 'F': -5, 'G': -2, 'H': 0, 'I': -2, 'K': 5, 'L': -3, 'M': 0, 'N': 1, 'P': -1,
          'Q': 1, 'R': 3, 'S': 0, 'T': 0, 'V': -2, 'W': -3, 'Y': -4},
    'L': {'A': -2, 'C': -6, 'D': -4, 'E': -3, 'F': 2, 'G': -4, 'H': -2, 'I': 2, 'K': -3, 'L': 6, 'M': 4, 'N': -3,
          'P': -3, 'Q': -2, 'R': -3, 'S': -3, 'T': -2, 'V': 2, 'W': -2, 'Y': -1},
    'M': {'A': -1, 'C': -5, 'D': -3, 'E': -2, 'F': 0, 'G': -3, 'H': -2, 'I': 2, 'K': 0, 'L': 4, 'M': 6, 'N': -2,
          'P': -2, 'Q': -1, 'R': 0, 'S': -2, 'T': -1, 'V': 2, 'W': -4, 'Y': -2},
    'N': {'A': 0, 'C': -4, 'D': 2, 'E': 1, 'F': -3, 'G': 0, 'H': 2, 'I': -2, 'K': 1, 'L': -3, 'M': -2, 'N': 2, 'P': 0,
          'Q': 1, 'R': 0, 'S': 1, 'T': 0, 'V': -2, 'W': -4, 'Y': -2},
    'P': {'A': 1, 'C': -3, 'D': -1, 'E': -1, 'F': -5, 'G': 0, 'H': 0, 'I': -2, 'K': -1, 'L': -3, 'M': -2, 'N': 0,
          'P': 6, 'Q': 0, 'R': 0, 'S': 1, 'T': 0, 'V': -1, 'W': -6, 'Y': -5},
    'Q': {'A': 0, 'C': -5, 'D': 2, 'E': 2, 'F': -5, 'G': -1, 'H': 3, 'I': -2, 'K': 1, 'L': -2, 'M': -1, 'N': 1, 'P': 0,
          'Q': 4, 'R': 1, 'S': -1, 'T': -1, 'V': -2, 'W': -5, 'Y': -4},
    'R': {'A': -2, 'C': -4, 'D': -1, 'E': -1, 'F': -4, 'G': -3, 'H': 2, 'I': -2, 'K': 3, 'L': -3, 'M': 0, 'N': 0,
          'P': 0, 'Q': 1, 'R': 6, 'S': 0, 'T': -1, 'V': -2, 'W': 2, 'Y': -4},
    'S': {'A': 1, 'C': 0, 'D': 0, 'E': 0, 'F': -3, 'G': 1, 'H': -1, 'I': -1, 'K': 0, 'L': -3, 'M': -2, 'N': 1, 'P': 1,
          'Q': -1, 'R': 0, 'S': 2, 'T': 1, 'V': -1, 'W': -2, 'Y': -3},
    'T': {'A': 1, 'C': -2, 'D': 0, 'E': 0, 'F': -3, 'G': 0, 'H': -1, 'I': 0, 'K': 0, 'L': -2, 'M': -1, 'N': 0, 'P': 0,
          'Q': -1, 'R': -1, 'S': 1, 'T': 3, 'V': 0, 'W': -5, 'Y': -3},
    'V': {'A': 0, 'C': -2, 'D': -2, 'E': -2, 'F': -1, 'G': -1, 'H': -2, 'I': 4, 'K': -2, 'L': 2, 'M': 2, 'N': -2,
          'P': -1, 'Q': -2, 'R': -2, 'S': -1, 'T': 0, 'V': 4, 'W': -6, 'Y': -2},
    'W': {'A': -6, 'C': -8, 'D': -7, 'E': -7, 'F': 0, 'G': -7, 'H': -3, 'I': -5, 'K': -3, 'L': -2, 'M': -4, 'N': -4,
          'P': -6, 'Q': -5, 'R': 2, 'S': -2, 'T': -5, 'V': -6, 'W': 17, 'Y': 0},
    'Y': {'A': -3, 'C': 0, 'D': -4, 'E': -4, 'F': 7, 'G': -5, 'H': 0, 'I': -1, 'K': -4, 'L': -1, 'M': -2, 'N': -2,
          'P': -5, 'Q': -4, 'R': -4, 'S': -3, 'T': -3, 'V': -2, 'W': 0, 'Y': 10}
}


class NeedlemanWansch():
    def __init__(self, seq_left, seq_top):

        self.seqs = []
        self.seq_left = list(seq_left)
        self.seq_top = list(seq_top)
        self.lleft = len(self.seq_left)
        self.ltop = len(self.seq_top)
        self.matrix = np.zeros((len(seq_left) + 1, len(seq_top) + 1))
        self.directions = [[[] for i in range(len(seq_top) + 1)] for j in range(len(seq_left) + 1)]
        self.maxs = []
        self.alignment_score = 0

    def calc_score(self, i, j):
        left_gap = self.matrix[i][j - 1] - 9
        right_gap = self.matrix[i - 1][j] - 9
        align = self.matrix[i - 1][j - 1] + PAM250[self.seq_left[i - 1]][self.seq_top[j - 1]]
        mx = max(left_gap, right_gap, align)
        if left_gap == mx:
            self.directions[i][j].append("h")
        if right_gap == mx:
            self.directions[i][j].append("v")
        if align == mx:
            self.directions[i][j].append("d")
        return mx

    def fill_matrix(self):
        mx = 0
        for i in range(1, self.lleft + 1):
            for j in range(1, self.ltop + 1):
                self.matrix[i][j] = self.calc_score(i, j)
                if i==self.lleft or j==self.ltop:
                    if mx < self.matrix[i][j]:
                        self.maxs = [(i, j)]
                        mx = self.matrix[i][j]
                    elif mx == self.matrix[i][j]:
                        self.maxs.append((i, j))
        self.alignment_score = mx
        # print(self.matrix)
        # for row in self.directions:
        #     print(row)
        # print(self.maxs)

    def recursion(self, seq_left, seq_top, i, j):
        if i == 0:
            # print("i==0",seq_left,seq_top)
            self.seqs.append((str("-" * j) + seq_left[::-1], "".join(self.seq_top[:j]) + seq_top[::-1]))
            return
        if j == 0:
            # print("j==0", seq_left, seq_top)
            self.seqs.append(("".join(self.seq_left[:i]) + seq_left[::-1], str("-" * i) + seq_top[::-1]))
            return

        for dir in self.directions[i][j]:
            if dir == "h":
                self.recursion(seq_left + "-", seq_top + self.seq_top[j - 1], i, j - 1)
            elif dir == "v":
                self.recursion(seq_left + self.seq_left[i - 1], seq_top + "-", i - 1, j)
            elif dir == "d":
                self.recursion(seq_left + self.seq_left[i - 1], seq_top + self.seq_top[j - 1], i - 1, j - 1)

    def printScores(self):
        sortedSeq = [i[0] + i[1] for i in self.seqs]
        sortedSeq.sort()
        print(int(self.alignment_score))
        for i in sortedSeq:
            print(i[0:int(len(i) / 2)])
            print(i[int(len(i) / 2):])

    def backtrack(self):
        for mx in self.maxs:
            left_chars = self.lleft - mx[0]  # 9-9=0
            top_chars = self.ltop - mx[1]  # 6-5=1
            max_chars = max(left_chars, top_chars)
            seq_left = str("-" * (max_chars - left_chars)) + str("".join(self.seq_left[mx[0]:]))[::-1]
            seq_top = str("-" * (max_chars - top_chars)) + str("".join(self.seq_top[mx[1]:]))[::-1]
            # print(seq_left, seq_top)
            self.recursion(seq_left, seq_top, mx[0], mx[1])
        self.printScores()


if __name__ == '__main__':
    dp = NeedlemanWansch(input(), input())
    dp.fill_matrix()
    dp.backtrack()

# HEAGAWGHE
# PAWHEA


# FKHMEDPLE
# FMDTPLNE

# FKHMEDPC
# FMDTPL

# ACGTACGTACGTCCCCCCCCC
# ACTGACGTCCCCCWWWWCCC

# ACTATATTATATATA
# ACTATATATATATA


# ACTATATTATA
# ACTATATATA


# GTCCCCCCCCC
# GTCCCCCWWWWCCC


# ACGTWWW
# ACGTCCC