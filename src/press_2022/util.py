from time import time
from numpy import *
import torch
import argparse
import re
import pickle

N_COS = 4  # number of cosine templates (usually 4)
alphabet = "ACGT"

# verify cuda and set device
if torch.cuda.is_available() :
    device = torch.device(f"cuda")
    cudaname = torch.cuda.get_device_name()
else:
    # raise RuntimeError("Required GPU not found! Exiting.")
    print("GPU not found! Falling back to CPU...")
    device = torch.device(f"cpu")


def decode(x):  # convert list of ints to string
    s = "".join([alphabet[xx] for xx in x])
    return s


def encode(st):  # convert a string into a list of ints
    x = [alphabet.index(ss) for ss in st]
    return x


def seqtomer(seq):  # return list of trimers in a seq
    ans = [int(16 * seq[k] + 4 * seq[k + 1] + seq[k + 2]) for k in range(len(seq) - 2)]
    return ans


def mertobin(mer):  # trimer list to occupancy uint64 bitmap
    ibin = 0
    for m in mer:
        ibin |= (1 << m)
    return ibin


def makeerrors(seq, srate, irate, drate):  # error mode applied to a sequence seq
    # note: modifies (also returns) seq
    n = len(seq)
    # substitutions
    ns = random.binomial(n, srate * 1.3333)  # 3/4 of substitutions are "effective"
    ndx = random.randint(low=0, high=n, size=ns)
    vals = random.randint(low=0, high=4, size=ns)
    seq[ndx] = vals
    # deletions
    nd = random.binomial(n, drate)
    ndx = random.randint(low=0, high=n, size=nd)
    seq = delete(seq, ndx)
    # insertions (at specified rate into smaller seq)
    ni = random.binomial(len(seq), irate)
    ndx = random.randint(low=0, high=len(seq) + 1, size=ni)
    vals = random.randint(low=0, high=4, size=ni)
    seq = insert(seq, ndx, vals)
    # pad or truncate to original length
    nn = len(seq)
    if nn > n:
        seq = seq[:n]
    elif nn < n:
        seq = concatenate((seq, random.randint(low=0, high=4, size=n - nn)))
    return seq


m0 = uint64(0x5555555555555555)  # binary: 0101...
m1 = uint64(0x3333333333333333)  # binary: 00110011..
m2 = uint64(0x0f0f0f0f0f0f0f0f)  # binary:  4 zeros,  4 ones ...
m3 = uint64(0x00ff00ff00ff00ff)  # binary:  8 zeros,  8 ones ...
m4 = uint64(0x0000ffff0000ffff)  # binary: 16 zeros, 16 ones ...
m5 = uint64(0x00000000ffffffff)  # binary: 32 zeros, 32 ones


def popcount(x):
    # https://github.com/google/jax/blob/6c8fc1b031275c85b02cb819c6caa5afa002fa1d/jax/lax_reference.py#L121-L150
    x = (x & m0) + ((x >> 1) & m0)  # put count of each  2 bits into those  2 bits
    x = (x & m1) + ((x >> 2) & m1)  # put count of each  4 bits into those  4 bits
    x = (x & m2) + ((x >> 4) & m2)  # put count of each  8 bits into those  8 bits
    x = (x & m3) + ((x >> 8) & m3)  # put count of each 16 bits into those 16 bits
    x = (x & m4) + ((x >> 16) & m4)  # put count of each 32 bits into those 32 bits
    x = (x & m5) + ((x >> 32) & m5)  # put count of each 64 bits into those 64 bits
    return x


try:
    assert torch.cuda.is_available()
    from popcll_torch import popcll

    # https://github.com/iamgroot42/popcll_torch , thanks to Anshuman Suri (as9rw@virginia.edu)!
    mypopcount = popcll
except:
    mypopcount = popcount
    print("popcll_torch not found, so using slower popcount")


def find_runs(x):
    # Find runs of consecutive items in an array.
    # credit: https://gist.github.com/alimanfoo/c5977e87111abe8127453b21204c1065
    n = x.shape[0]
    # find run starts
    loc_run_start = empty(n, dtype=bool)
    loc_run_start[0] = True
    not_equal(x[:-1], x[1:], out=loc_run_start[1:])
    run_starts = nonzero(loc_run_start)[0]
    # find run values
    run_values = x[loc_run_start]
    # find run lengths
    run_lengths = diff(append(run_starts, n))
    return run_values, run_starts, run_lengths


def find_max_run(seq):
    # Finds run with maximum number of repetitions of a substring
    # Added by Uphoff 2023
    # TODO optimize
    max_run = 0
    possible_runs = set()

    for i, val in enumerate(seq):
        new_runs = set()

        for run in possible_runs:
            # [start_id, substring_length, substring_position, repetitions]
            max_run = max(max_run, run[3])

            if run[0] + run[1] == i:
                new_runs.add((run[0], run[1] + 1, 0, 1))

            if seq[run[0] + run[2]] == val:
                new_runs.add((run[0], run[1], (run[2] + 1) % run[1], run[3] + int((run[2] + 1) / run[1])))

        del possible_runs
        possible_runs = new_runs
        possible_runs.add((i, 1, 0, 1))

    return max_run


def chemfilter(seq, homomax=3, atmax=22, cgmax=22):
    # returns whether seq satisfies chemistry constraints
    bc = bincount(seq, minlength=4)
    if (bc[0] + bc[3] > atmax) or (bc[1] + bc[2] > cgmax): return False
    # Original homopolymer restrictions
    # _, _, run_lengths = find_runs(seq)
    # if max(run_lengths) > homomax: return False
    # Altered, more strict restrictions
    if find_max_run(seq) > homomax: return False
    return True


def allcoses(mer, tcosvecs):  # correlate a mer against all the cosine templates
    mmer = torch.LongTensor(mer).to(device)
    ncos = tcosvecs.shape[0]
    cosvec = torch.zeros(ncos, 64, dtype=torch.float, device=device)
    for k in range(ncos):
        source = tcoses[k, torch.arange(len(mmer), dtype=torch.long, device=device)]
        cosvec[k, :].index_add_(0, mmer, source)
    return torch.sum(torch.unsqueeze(cosvec, dim=1) * tcosvecs, dim=2)


def prank(arr, descending=False):  # returns rank of each element in torch array
    argsrt = torch.argsort(arr, descending=descending)
    rank = torch.zeros(arr.shape, dtype=torch.float, device=device)
    rank[argsrt] = torch.arange(len(argsrt), dtype=torch.float, device=device)
    return rank


class ApproximateLevenshtein:
    def __init__(s, M, N, Q, zsub, zins, zdel, zskew):
        torch.set_grad_enabled(False)  # just in case not done elsewhere!
        s.M = M  # length of seq1
        s.N = N  # length of each seq2
        s.Q = Q  # number of seq2s
        (s.zsub, s.zins, s.zdel, s.zskew) = (zsub, zins, zdel, zskew)
        s.tab = torch.zeros(N + 1, Q, device=device)

    def __call__(s, seq1, seq2):
        assert (len(seq1) == s.M) and (seq2.shape[1] == s.N) and (seq2.shape[0] == s.Q)
        s.tab[:, :] = (s.zskew * torch.arange(s.N + 1., device=device)).unsqueeze(1)  # force broadcast
        for i in range(1, s.M + 1):
            diag = s.tab[:-1, :] + torch.where(seq1[i - 1] == seq2.t(), 0., s.zsub)  # diagonal move
            s.tab[0, :] += s.zskew
            s.tab[1:, :] += s.zdel  # down move
            s.tab[1:, :] = torch.minimum(s.tab[1:, :], diag)  # or diag if better
            s.tab[1:, :] = torch.minimum(s.tab[1:, :], s.tab[:-1, :] + s.zins)  # right move
            s.tab[1:, :] = torch.minimum(s.tab[1:, :], s.tab[:-1, :] + s.zins)  # repeat (>= 0 times) as you can afford
        # N.B.: M >= N gives better approx than N > M, so change arg order accordingly
        return s.tab[s.N, :]


class ParallelLevenshtein:
    def __init__(s, M, N, Q, zsub, zins, zdel, zskew):
        torch.set_grad_enabled(False)  # just in case not done elsewhere!
        s.M = M  # length of seq1
        s.N = N  # length of each seq2
        s.Q = Q  # number of seq2s
        (s.zsub, s.zins, s.zdel, s.zskew) = (zsub, zins, zdel, zskew)
        MN1 = M + N + 1
        s.blue = torch.zeros(Q, MN1, MN1, device=device)
        s.bluef = s.blue.view(Q, MN1 * MN1)
        s.ndxr = torch.zeros(M * N, dtype=int, device=device)  # index of mer matches array into flat blue
        for m in torch.arange(M, device=device):
            for n in torch.arange(N, device=device):
                s.ndxr[n + N * m] = (3 * M + 2 * N + 2) + (M + N) * m + (M + N + 2) * n
        s.lls = torch.zeros(MN1 + 1, dtype=torch.int, device=device)
        s.rrs = torch.zeros(MN1 + 1, dtype=torch.int, device=device)
        for i in range(2, MN1 + 1):
            s.lls[i] = abs(M - i + 1) + 1
            s.rrs[i] = (M + N - 1) - abs(- i + 1 + N)

    def __call__(s, seq1, sseq2):  # single seq1, tensor of sseq2s
        assert (len(seq1) == s.M) and (sseq2.shape[1] == s.N) and (sseq2.shape[0] == s.Q)
        (M1, N1, MN, MN1, MN2) = (s.M + 1, s.N + 1, s.M + s.N, s.M + s.N + 1, s.M + s.N + 2)
        abmatch = (seq1.view(1, s.M, 1) != sseq2.view(s.Q, 1, s.N)).type(torch.float) * s.zsub
        s.bluef[:, s.ndxr] = abmatch.view(s.Q, s.M * s.N)
        s.bluef[:, torch.arange(s.M, MN2 * N1, MN2)] = (s.zskew * torch.arange(N1, device=device)).unsqueeze(0)
        s.bluef[:, torch.arange(s.M, MN * M1, MN)] = (s.zskew * torch.arange(M1, device=device)).unsqueeze(0)
        for k in torch.arange(2, MN1, device=device):
            ll = s.lls[k]
            rr = s.rrs[k]
            slis = torch.arange(ll, rr + 1, 2, device=device)
            s.blue[:, k, slis] = torch.minimum(
                s.blue[:, k, slis] + s.blue[:, k - 2, slis],
                torch.minimum(
                    s.blue[:, k - 1, slis - 1] + s.zdel,
                    s.blue[:, k - 1, slis + 1] + s.zins
                )
            )
        return s.blue[:, -1, s.N]


def measure_execution_time(func):
    def wrapper_function(*args, **kwargs):
        start = time()
        res = func(*args, **kwargs)
        end = time()

        print(f"Function {func.__name__} took {(end - start):.4f} s")
        return res

    return wrapper_function
