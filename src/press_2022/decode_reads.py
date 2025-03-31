from util import *
from generate_barcodes import load_barcodes
from timeit import default_timer as timer


def get_reads(readsfilename):
    with open(readsfilename, 'r') as IN:
        stringreads=IN.readlines()
    reads = []
    for read in stringreads :
        reads.append(encode(read[:-1])) # lose the \n
    reads = array(reads)

    return reads


def move_to_gpu(barcode_dict):
    """
    Move the code (assumed loaded above) to the GPU.
    :param readsfilename:
    :param barcode_dict:
    :return:
    """
    torch.cuda.empty_cache()  # TODO
    barcode_dict_gpu = dict()

    torch.set_grad_enabled(False)
    barcode_dict_gpu["allseqs"] = torch.tensor(array(barcode_dict["allseqs"]), device=device)
    barcode_dict_gpu["alltrimers"] = torch.tensor(array(barcode_dict["alltrimers"]), device=device)
    barcode_dict_gpu["allbitmaps"] = torch.tensor(barcode_dict["allbitmaps"].astype(int64), dtype=torch.int64, device=device)
    barcode_dict_gpu["coses"] = torch.tensor(barcode_dict["coses"], dtype=torch.float, device=device)
    barcode_dict_gpu["cosvecs"] = barcode_dict["cosvecs"].to(device)
    barcode_dict_gpu["N"] = barcode_dict["N"]
    barcode_dict_gpu["M"] = barcode_dict["M"]

    return barcode_dict_gpu


@measure_execution_time
def call_reads(reads, barcode_dict_gpu, pos=100, L=8):
    print(f"Calling reads with L = {L} and pos = {pos}...")
    N = barcode_dict_gpu["N"]
    M = barcode_dict_gpu["M"]
    r = reads.shape[0]

    # mydist = ApproximateLevenshtein(M, M, pos, 1., 1., 1., 1.)
    mydist = ParallelLevenshtein(M,M,pos, 1.,1.,1.,1.)

    # Performance tracking
    exec_time_triage = 0
    exec_time_levenshtein = 0

    Ncos = barcode_dict_gpu["cosvecs"].shape[0]
    dists = torch.zeros(Ncos + 1, N, dtype=torch.float, device=device)  # will hold distances for each read
    allrank = torch.zeros(Ncos + 1, N, dtype=torch.float, device=device)
    best = torch.zeros(r, dtype=torch.int32, device=device)
    best_dists = torch.zeros(r, dtype=torch.int32, device=device)

    for j, seq in enumerate(reads[:r]):
        start = timer()  # Performance tracking
        # primary and secondary triage
        mer = seqtomer(seq)
        foo = int64(uint64(mertobin(mer)))  # need to cast 64 bits to a type known to torch
        seqbin = torch.tensor(foo, dtype=torch.int64, device=device)
        xored = torch.bitwise_xor(seqbin, barcode_dict_gpu["allbitmaps"])
        dists[0, :] = 64. - mypopcount(xored)  # all Hamming distances
        cosvec = torch.zeros(Ncos, 64, dtype=torch.float, device=device)
        for k in range(Ncos):
            cosvec[k, mer] = barcode_dict_gpu["coses"][k, torch.arange(len(mer), dtype=torch.long, device=device)]
        dists[1:, :] = torch.sum(torch.unsqueeze(cosvec, dim=1) * barcode_dict_gpu["cosvecs"],
                                 dim=2)  # all cosine distances
        for k in range(Ncos + 1):
            allrank[k, :] = prank(dists[k, :], descending=True)  # rank them all
        offset = 1.
        fm = torch.prod(offset + allrank, dim=0)
        fmargsort = torch.argsort(fm)

        end = timer()  # Performance tracking
        exec_time_triage += end - start
        start = timer()  # Performance tracking

        # Levenshtein distance
        tseq1 = torch.tensor(seq, device=device)
        tseq2 = barcode_dict_gpu["allseqs"][fmargsort[:pos], :]
        ans = mydist(tseq1, tseq2)
        ia = torch.argmin(ans)  # index into fmargsort of best
        ibest = fmargsort[ia]  # index of best codeword in codes
        best[j] = (ibest if ans[ia] <= L else -1)  # erasures returned as -1
        best_dists[j] = ans[ia]

        end = timer()  # Performance tracking
        exec_time_levenshtein += end - start

    print("Execution times:")
    print(f"- Triages: {exec_time_triage:.4f} s ({exec_time_triage * 1000 / r:.4f} ms per read)")
    print(f"- Levenshtein: {exec_time_levenshtein:.4f} s ({exec_time_levenshtein * 1000 / r:.4f} ms per read)")

    return best.cpu().numpy(), best_dists.cpu().numpy(), exec_time_triage, exec_time_levenshtein


def decode_reads(barcode_filename, reads_filename, decodes_filename, dists_filename=None, pos=100, L=8):
    """
    DECODING Primary and secondary triage, followed by Levenshtein
    :param barcode_filename:
    :param reads_filename:
    :param decodes_filename:
    :param pos: number passed from triage to Levenshtein
    :param L: Levenshtein score greater than which is called an erasure
    :return:
    """
    barcode_dict = load_barcodes(barcode_filename)
    reads = get_reads(reads_filename)
    barcode_dict_gpu = move_to_gpu(barcode_dict)
    print("Moved barcodes to GPU")
    r = reads.shape[0]

    torch.set_grad_enabled(False)

    best, best_dists, _, _ = call_reads(reads, barcode_dict_gpu, pos=pos, L=L)

    print("done with decoding.")

    with open(decodes_filename, "w") as out_file:
        for i in range(r):
            out_file.write(str(best[i]) + "\n")

    if dists_filename is not None:
        with open(decodes_filename, "w") as out_file:
            for i in range(r):
                out_file.write(str(best_dists[i]) + "\n")

    print("wrote labels to file " + decodes_filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='decode_reads',
        description='Decodes reads and writes labels to file',
        epilog='TODO')

    parser.add_argument("-b", "--barcodes", required=True)
    parser.add_argument("-r", "--reads", required=True)
    parser.add_argument("-o", "--outfile", required=True)
    parser.add_argument("-d", "--distancefile", required=False, default=None)
    parser.add_argument("-p", "--pos", type=int, required=False, default=100)
    parser.add_argument("-L", "--max-dist", type=int, required=False, default=8)

    args = parser.parse_args()

    decode_reads(args.barcodes, args.reads, args.outfile,
                 dists_filename=args.distancefile, pos=args.pos, L=args.max_dist)
    torch.cuda.empty_cache()  # TODO
