from util import *


def generate_random_barcodes(filename, N=100000, M=34, DoChemFilter=True):
    """
    Generate and write the codes
    :param filename: output file path (produces file length ~ M*N)
    :param N: number of codewords (this is a small value for demo purposes, 1000000 is feasible)
    :param M: length of a codeword (nt)
    :param DoChemFilter: filter for homopolymer runs, AT and CG content?
    :return:
    """
    # Restrictions for chem. filter
    HOMO_MAX = 2    # Changed from 3 to 2
    AT_MAX = int(0.6 * M)  # Changed from 0.66 to 0.6
    GC_MAX = int(0.6 * M)  # Changed from 0.66 to 0.6

    torch.set_grad_enabled(False)
    if DoChemFilter :
        fac = 8. # should work in all expected cases !Changed from 2 to 8
        Ntrial = int(fac*N)
        allcands = random.randint(low=0,high=4,size=(Ntrial,M))
        goodcands = array([chemfilter(x,HOMO_MAX,AT_MAX,GC_MAX) for x in allcands], dtype=bool)
        allcodes = allcands[goodcands][:N] # if this or next throws an error, increase fac, but shouldn't happen
        if allcodes.shape != (N,M) : raise RuntimeError(f"Only {allcodes.shape[0]} barcodes found! see code, increase fac (this shouldn\'t happen')")
    else :
        allcodes = random.randint(low=0,high=4,size=(N,M))
    with open(filename, 'w') as OUT:
        for code in allcodes :
            OUT.write(decode(code) + '\n')
    print(f"wrote {N} codewords of length {M} in ascii to {filename}")


@measure_execution_time
def load_barcodes(filename):
    if filename[-4:] == ".pkl":
        # Pickle file with vectors
        with open(filename, 'rb') as IN:
            pickledict = pickle.load(IN)
        (N, M, allseqs, alltrimers, allbitmaps, coses, cosvecs) = \
            [pickledict[x] for x in ('N', 'M', 'allseqs', 'alltrimers', 'allbitmaps', 'coses', 'cosvecs')]
        print(f"Loaded pickled barcodes with {N} codewords of length {M} from {filename}.")
    else:
        with open(filename) as IN:
            codes = IN.readlines()  # have \n EOL
            for i, code in enumerate(codes):
                codes[i] = encode(code[:-1])

            pickledict = pickle_barcodes(codes)

    return pickledict


def pickle_barcodes(codes):
    N = len(codes)
    M = 0 if N == 0 else len(codes[0][:])

    print(f"Loaded {N} barcodes of length {M}...")

    allseqs = []
    alltrimers = []
    allbitmaps = zeros(N, dtype=uint64)
    cosvecs = torch.zeros((N_COS, N, 64), dtype=torch.float)
    coses = zeros((N_COS, M))
    for k in range(N_COS):
        coses[k, :] = cos(pi * arange(M) * (k + 1.) / (M - 1.))
    tcoses = torch.tensor(coses, dtype=torch.float)
    for i, code in enumerate(codes):
        seq = code
        allseqs.append(seq)
        mer = seqtomer(seq)
        mmer = torch.LongTensor(mer)
        alltrimers.append(mer)
        allbitmaps[i] = mertobin(mer)
        for k in range(N_COS):
            source = tcoses[k, arange(M - 2)]
            cosvecs[k, i, :].index_add_(0, mmer, source)

    pickledict = {"N": N, "M": M, "allseqs": allseqs,
                  "alltrimers": alltrimers, "allbitmaps": allbitmaps, "coses": coses, "cosvecs": cosvecs}

    print(f"Created vectors for barcodes...")

    return pickledict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='generate_barcodes',
        description='Generates random barcodes with given attributes',
        epilog='TODO')

    parser.add_argument("-o", "--outfile", required=True)
    parser.add_argument("-n", "--number", type=int, required=True)
    parser.add_argument("-l", "--length", type=int, required=True)
    parser.add_argument("-f", "--filter", required=False, action='store_true', default=True)

    args = parser.parse_args()

    generate_random_barcodes(args.outfile, N=args.number, M=args.length, DoChemFilter=args.filter)
