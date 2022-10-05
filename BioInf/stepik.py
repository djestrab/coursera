import itertools
import numpy as np
import random as rand

# Course 1

def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = list(range(r))
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)

def pattern_count(pattern, text):
    """Returns a number of occurences of pattern in text."""

    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i+len(pattern)] == pattern:
            count = count + 1
    return count


def PatternStart(Pattern, Text):
    """ Returns starting indexes of Pattern in Text """

    StartIndex = Text.find(Pattern)
    StartingPositions = []
    while StartIndex < len(Text) - len(Pattern):
        FindResult = Text.find(Pattern, StartIndex)
        if FindResult != -1:
            StartingPositions.append(FindResult)
            StartIndex = FindResult + 1
        else:
            return StartingPositions
            #return (' '.join(str(e) for e in StartingPositions)

def PatternMatching(Pattern, Genome):
  l1 = len(Pattern)
  l2 = len(Genome)
  ans=[]
  for i in range(0, l2-l1+1):
    if( Pattern == Genome[i:i+l1]):
      ans.append(i)
  return ans


def FrequentWords(Text, k):
    """Returns all the most frequent k-long patterns in Text."""

    i = 0
    PatternFrequencies = {}
    Pattern = ""
    MaxCount = 0

    while i < len(Text) - k:
        Pattern = Text[i:i + k]
        PatternFrequencies[Pattern] = PatternCount(Pattern, Text)
        i += 1

    MaxCount = max(PatternFrequencies.items(), key=lambda x: x[1])
    MostFrequentKMers = list()
    for key, value in PatternFrequencies.items():
        if value == MaxCount[1]:
            MostFrequentKMers.append(key)

    return MostFrequentKMers


def ReverseComplement(Text):
    """ Make a reverse complement pattern for Text. """

    complementPattern = []
    Complements = {"C": "G", "G": "C", "A": "T", "T": "A"}
    for n in Text.upper():
        complementPattern.append(Complements[n])
    complementPattern.reverse()
    return complementPattern

def PatternToNumber(Pattern):
    """  """

    if len(Pattern) == 0:
        return 0

    NuklIndex = ["A", "C", "G", "T"]
    Symbol = Pattern[-1:]
    Prefix = Pattern[:-1]
    return (4 * PatternToNumber(Prefix)) + NuklIndex.index(Symbol)

def NumberToPattern(index, k):
    """  """

    NuklIndex = ["A", "C", "G", "T"]
    if k == 1:
        return NuklIndex[index]
    prefixIndex = int(index / 4)
    r = index % 4
    symbol = NuklIndex[r]
    PrefixPattern = NumberToPattern(prefixIndex, k-1)
    return PrefixPattern + symbol

def ComputingFrequencies(Text, k):
    """ Computes occurences of each k-long subtext of Text.
    In the first range creates array for each k-mer sorted lexicographically
    in which it in the second range stores number of occurences. """
    frequency_array = []
    i = 0
    while i < len(Text) - k:
        for i in range(0, 4**k):
            frequency_array.append(0)
        for i in range(0, len(Text) - (k-1)):
            pattern = Text[i:i+k]
            j = PatternToNumber(pattern)
            frequency_array[j] += 1
    return frequency_array


def ClumpFinding(Genome, k, L, t):
    """ Returns all k-long patterns forming (L,t)-clumps in Genome.
    L = clump lenght, t = number of occurences of Pattern in clump """

    FrequentPatterns = []
    Clump = []
    for i in range(0, 4**k):
        Clump.append(0)
    for i in range(0, len(Genome) - L + 1):
        Text = Genome[i:i+L]
        FrequencyArray = ComputingFrequencies(Text, k)
        for index in range(0, 4**k):
            if FrequencyArray[index] >= t:
                Clump[index] = 1
    for i in range(0, 4**k):
        if Clump[i] == 1:
            Pattern = NumberToPattern(i, k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns


def get_skew_minimum(genome):
    """ Returns indices of minimum values of G-C difference """

    skew = [0]
    skew_min = 0
    skew_min_indices = []
    for n in genome:
        if n == "G":
            skew.append(skew[-1] + 1)
        elif n == "C":
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])

        if skew[-1] == skew_min:
            skew_min_indices.append(len(skew) - 1)
        if skew[-1] < skew_min:
            skew_min_indices = [len(skew) - 1]
        skew_min = min(skew)
    return skew_min_indices


def hamming_distance(p, q):
    """ Returns number of mismatches between two strings """

    h = 0
    for i, char in enumerate(p):
        if char != q[i]:
            h += 1
    return h


def approx_pattern_matching(pattern, text, d):
    """ Returns all starting positions where pattern appears as a substring of text with at most d mismatches."""

    l1 = len(pattern)
    l2 = len(text)
    ans = []
    for i in range(0, l2-l1+1):
        if hamming_distance(pattern, text[i:i+l1]) <= d:
            ans.append(i)
    return ans


def approx_pattern_count(pattern, text, d):
    """ Returns number of times where pattern appears as a substring of text with at most d mismatches."""

    l1 = len(pattern)
    l2 = len(text)
    count = 0
    for i in range(0, l2-l1+1):
        if hamming_distance(pattern, text[i:i+l1]) <= d:
            count += 1
    return count


def immediate_neighbors(pattern):
    """  """

    nucleotides = "ACGT"
    neighborhood = [pattern]
    for i in range(0, len(pattern)):
        symbol = pattern[i]
        for n in nucleotides:
            if n != symbol:
                neighborhood.append(pattern[:i] + n + pattern[i + 1:])
    return neighborhood


def neighbors(pattern, d):
    """ Returns all strings similar to pattern with up to "d" mismatches. """

    nucleotides = "ACGT"
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return ["A", "C", "G", "T"]
    neighborhood = []
    suffix_neighbors = neighbors(pattern[1:], d)
    for Text in suffix_neighbors:
        if hamming_distance(pattern[1:], Text) < d:
            for n in nucleotides:
                neighborhood.append(n + Text)
        else:
            neighborhood.append(pattern[0] + Text)
    return neighborhood


def approx_frequent_words(text, k, d):
    """ Returns most frequent k-long patterns with its reverse complements in text with at most d mismatches."""

    frequent_patterns = []
    neighborhoods = []
    index = []
    count = []

    for i in range(0, len(text) - k + 1):  # add all neighbors to be checked
        neighborhoods += neighbors(text[i: i + k], d)
    for i in range(0, len(neighborhoods)):  # convert patterns to indices
        pattern = neighborhoods[i]
        index.append(PatternToNumber(pattern))
        count.append(1)
    index.sort()
    for i in range(0, len(neighborhoods) - 1):  # get counts for all patterns(indices)
        if index[i] == index[i + 1]:
            count[i + 1] = count[i] + 1
    max_count = max(count)
    for i in range(0, len(neighborhoods)):  # convert most frequent indices back to patterns
        if count[i] == max_count:
            pattern = NumberToPattern(index[i], k)
            frequent_patterns.append(pattern)
    return frequent_patterns

    # MotifEnumeration(Dna, k, d)
    #     Patterns ← an empty set
    #     for each k-mer Pattern in Dna
    #         for each k-mer Pattern’ differing from Pattern by at most d mismatches
    #             if Pattern' appears in each string from Dna with at most d mismatches
    #                 add Pattern' to Patterns
    #     remove duplicates from Patterns
    #     return Patterns


def motif_enumeration(dna, k, d):
    """ brute force motif finding """
    neighborhoods = []

    for string in dna:
        original_patterns = []
        neighborhood = []
        for i in range(0, len(string) - k + 1):
            pattern = (string[i:i+k])
            original_patterns.append(pattern)
        for p in original_patterns:
            neighborhood += neighbors(p, d)
        neighborhoods.append(neighborhood)

    patterns = set(neighborhoods[0])
    for s in neighborhoods[1:]:
        patterns.intersection_update(s)
    patterns = list(patterns)
    patterns.sort()
    return patterns


def pattern_string_distance(pattern, dna):
    """ Returns sum of hamming distances between pattern and strings of dna. Input dna as a list of strings """
    k = len(pattern)
    d = 0
    for string in dna:
        h_distance = 100
        kmers = []
        for i in range(0, len(string) - k + 1):
            kmers.append(string[i:i+k])
        for kmer in kmers:
            if h_distance > hamming_distance(pattern, kmer):
                h_distance = hamming_distance(pattern, kmer)
        d += h_distance
    return d


def median_string(dna, k):
    """ Returns a k-mer Pattern that minimizes d(Pattern, Dna) among all possible choices of k-mers.
    """

    d = 10000
    median = str
    for i in range(0, 4**k - 1):
        kmer = NumberToPattern(i, k)
        psd = pattern_string_distance(kmer, dna)
        if d > psd:
            d = psd
            median = kmer
    return median


def profile_probable_kmer(text, k, profile):
    """ finds most probable k-mer of length k in text matching the profile """

    winners = [text[0:k]]
    p_winner = 0
    nucleotides = ["A", "C", "G", "T"]

    for i in range(0, len(text)-k+1):
        kmer = text[i:i+k]
        p = 1
        for column, symbol in enumerate(kmer):
            row = nucleotides.index(symbol)
            p *= profile[row][column]
        if p > p_winner:
            p_winner = p
            winners = [kmer]
        elif p == p_winner:
            winners.append(kmer)
    return winners


def generate_profile(patterns, k):
    """Generates profile matrix showing probability of each nucleotide on every position for a list of k-long patterns."""
    nucleotides = ["A", "C", "G", "T"]
    profile = np.zeros([len(nucleotides), k], dtype=float)
    patterns = [[p for p in pattern] for pattern in patterns]
    patterns_transposed = np.array(patterns).transpose()
    for i, n in enumerate(nucleotides):
        for xi, p in enumerate(patterns_transposed):
            profile[i][xi] = (list(p).count(n) + 1) / (len(patterns) + 4)
    return profile


def score(motifs):
    """  """

    score_sum = 0
    motifs = [[p for p in motif] for motif in motifs]
    motifs_transposed = np.array(motifs).transpose()
    motifs_transposed = motifs_transposed.tolist()
    for motif in motifs_transposed:
        max_count = 0
        for letter in motif:
            if motif.count(letter) > max_count:
                max_count = motif.count(letter)
        score_sum += len(motif) - max_count
    return score_sum


def greedy_motif_search(dna, k, t):
    """
    :param dna:
    :param k:
    :param t:
    :return: Applicable for k < 12
    """

    best_motifs = []
    first_string = str

    for index, string in enumerate(dna):
        best_motifs.append(string[0:k])
        if index == 0:
            first_string = string
    for i in range(0, len(first_string) - k + 1):
        motif = first_string[i:i+k]
        motifs = [motif]
        for ix in range(1, t):
            profile = generate_profile(motifs, k)
            motifs.append(str(profile_probable_kmer(dna[ix], k, profile)[0]))
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


def randomized_motif_search_core(dna, k, t):
    """  """

    motifs = []
    for string in dna:
        r = rand.randint(0, len(string) - k)
        motifs.append(string[r:r+k])
    best_motifs = motifs
    while True:
        profile = generate_profile(motifs, k)
        temp_motifs = []
        for string in dna:
            temp_motifs.append(str(profile_probable_kmer(string, k, profile)[0]))
        motifs = temp_motifs
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs


def randomized_motif_search(dna, k, t, n):
    """
    Initializes randomized_motif_search function for n-number of times
    :param dna: list of dna strings
    :param k: k-mer
    :param t: number of dna strings
    :param n: number of iterations
    :return: list of motifs with lowest score found in n-iterations
    """

    all_best_motifs = []
    x_score = 1000

    for i in range(0, n):
        checked_motifs = randomized_motif_search_core(dna, k, t)
        r_score = score(checked_motifs)
        if r_score < x_score:
            print(str(r_score))
            x_score = r_score
            all_best_motifs = checked_motifs
    for m in all_best_motifs:
        print(m)
    return all_best_motifs


def random(p):
    """
    :param p: list of probabilities
    :return: index of the value chosen
    """

    return rand.choices([i for i, pop in enumerate(p)], p)


# GibbsSampler(Dna, k, t, N)
#         randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
#         BestMotifs ← Motifs
#         for j ← 1 to N
#             i ← Random(t)
#             Profile ← profile matrix constructed from all strings in Motifs except for Motifi
#             Motifi ← Profile-randomly generated k-mer in the i-th sequence
#             if Score(Motifs) < Score(BestMotifs)
#                 BestMotifs ← Motifs
#         return BestMotifs


def gibbs_sampler_core(dna, k, t, n):
    """
    :param dna: list of dna strings
    :param k: k-mer
    :param t: number of dna strings
    :param n: number of iterations
    :return:
    """

    motifs = []
    for string in dna:  # randomly select k-mer motifs from every DNA string
        r = rand.randint(0, len(string) - k)
        motifs.append(string[r:r + k])
    best_motifs = motifs
    for j in range(1, n):
        i = rand.randint(0, t-1)   # randomly select motif for removal
        motifs_for_pattern = list(motifs)
        motifs_for_pattern.pop(i)
        profile = generate_profile(motifs_for_pattern, k)   # generate pattern for remaining motifs
        motifs[i] = profile_probable_kmer(dna[i], k, profile)[0]    # replace the randomly selected motif with a "better" one
        if (score(motifs)) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


def gibbs_sampler(dna, k, t, n, *starts):
    all_best_motifs = []
    x_score = 1000

    if not starts:
        starts = 1

    while starts > 0:
        checked_motifs = gibbs_sampler_core(dna, k, t, n)
        r_score = score(checked_motifs)
        if r_score < x_score:
            x_score = r_score
            all_best_motifs = checked_motifs
        starts -= starts
    return all_best_motifs

# Course 2


def composition(text, k):
    """ Generates the k-mer composition of a string. """

    kmers = []
    for i in range(0, len(text) - k + 1):
        kmers.append(text[i:i+k])
#    kmers.sort()
    return kmers


def path_to_genome(path):
    """
    :param path: list of paths
    :return: initial dna string
    """

    for p in path:
        if p == path[0]:
            dna = path[0][0:]
        else:
            dna = dna + p[-1]
    return dna


def overlap(patterns):
    """
    :param patterns: list of kmers
    :return: adjacency list (dictionary of lists) of nodes and edges of an overlap-graph
    """
    patterns.sort()
    prefixes = list(map(lambda x: x[:-1], patterns))
    print(prefixes)
    adj_dict = {}

    for i, p in enumerate(patterns):
        if p[1:] in prefixes:
            idxs = [i for i, e in enumerate(prefixes) if e == p[1:]]
            for ix in idxs:
                if p in adj_dict.keys():
                    adj_dict[p].append(patterns[ix])
                else:
                    adj_dict[p] = [patterns[ix]]
    return adj_dict


def debruijn(text, *k):
    """
    :param text: String which we want to make De Bruijn graph of
    :param k: lenght of k-mer making the graph edges
    :return: adjacency list (dictionary of lists) of edges
    """

    adj_dict = {}

    if type(text) == str:
        edges = composition(text, k)    # make k-long edges
        nodes = [text[0:k - 1]]  # make first node (k-1 long)

        for edge in edges:  # make other nodes out of edges-suffixes (k-1 long)
            nodes.append(edge[1:])
        for i, node in enumerate(nodes[:-1]):  # take each node and find adjacent nodes (only the very next)
            if nodes[i + 1][:-1] == node[1:]:
                if node in adj_dict.keys():
                    adj_dict[node].append(nodes[i + 1])
                else:
                    adj_dict[node] = [nodes[i + 1]]

    elif type(text) == list:
        for t in text:
            if t[:-1] not in adj_dict.keys():
                adj_dict[t[:-1]] = [t[1:]]
            else:
                adj_dict[t[:-1]].append(t[1:])
    else:
        raise Exception("Argument Text must be either string or list of kmers")

    adj_dict = dict(sorted(adj_dict.items()))
    return adj_dict


def eulerian_cycle(graph):
    """
    :param graph: The adjacency list of an Eulerian directed graph (2d-list)
    :return: An Eulerian cycle in this graph = List of nodes(indices)
    """

    stack = []
    location = int
    circuit = []

    # Mark the nodes as "used" / subtract
    # While there are unused nodes left, move to next "multiple" node and iterate


with open("/Users/d.j/Downloads/dataset_199_6.txt", "r") as xinput:
    xlines = xinput.readlines()
    xlines = [x.strip() for x in xlines]
    xtext = ",".join(xlines[1:])

with open('/Users/d.j/Downloads/output.txt', 'w') as output:
    for key, value in overlap(xlines).items():
        output.writelines(key + " -> " + ", ".join(value) + "\n")


