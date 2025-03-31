from util import *
from generate_barcodes import load_barcodes
from decode_reads import move_to_gpu, get_reads


@measure_execution_time
def get_triage_positions(barcode_dict_gpu, reads, answers):
    torch.set_grad_enabled(False)

    N = barcode_dict_gpu["N"]
    M = barcode_dict_gpu["M"]
    r = reads.shape[0]

    Ncos = barcode_dict_gpu["cosvecs"].shape[0]
    dists = torch.zeros(Ncos + 1, N, dtype=torch.float, device=device)  # will hold distances for each read
    
    allrank = torch.zeros(Ncos + 1, N, dtype=torch.float, device=device)
    
    positions = torch.zeros(N, dtype=torch.int32, device=device)

    for j, seq in enumerate(reads[:r]):
        
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
        # Heuristic ordering
        ordering = fmargsort

        # positions[j] = torch.where(ordering == answers[j])[0]
        pos = (ordering == answers[j]).nonzero(as_tuple=True)[0]
        positions[pos] += 1

    return positions.cpu().numpy()


def get_index_recall(barcode_dict_gpu, reads, answers):
    N = barcode_dict_gpu["N"]
    r = reads.shape[0]

    frequencies_positions = get_triage_positions(barcode_dict_gpu, reads, answers)
    print("done with triages.")

    # Calc recall
    recall = zeros(N, dtype=double, )

    frequency_total = 0
    for i in range(N):
        frequency_total += frequencies_positions[i]
        recall[i] = frequency_total / r

    return recall


def write_triage_recall(barcode_filename, reads_filename, answers_filename, recall_filename):
    """
    DECODING Primary and secondary triage, followed by Levenshtein
    :param barcode_filename:
    :param reads_filename:
    :param recall_filename:
    :param answers_filename:
    :return:
    """
    barcode_dict = load_barcodes(barcode_filename)
    reads = get_reads(reads_filename)
    answers = genfromtxt(answers_filename, dtype=int)
    barcode_dict_gpu = move_to_gpu(barcode_dict)

    recall = get_index_recall(barcode_dict_gpu, reads, answers)

    with open(recall_filename, "w") as out_file:
        for i in range(N):
            out_file.write(str(recall[i]) + "\n")

    print("wrote recall per position to file " + recall_filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='calc_triage_recall',
        description='Calculates the recall per position achieved by the triage and write to file',
        epilog='TODO')

    parser.add_argument("-b", "--barcodes", required=True)
    parser.add_argument("-r", "--reads", required=True)
    parser.add_argument("-a", "--answers", required=True)
    parser.add_argument("-o", "--outfile", required=True)

    args = parser.parse_args()

    write_triage_recall(args.barcodes, args.reads, args.answers, args.outfile)
    torch.cuda.empty_cache()  # TODO
