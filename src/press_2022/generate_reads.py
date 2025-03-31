from util import *
from generate_barcodes import load_barcodes

# PRELIMINARIES Step 2: make simulated data, write "reads" to ascii file,


def generate_reads(barcodesfilename, readsfilename, answersfilename, Q=10000, srate=0.033, irate=0.033, drate=0.033, nave=4.):
    """
    "answers (indices of true codeword)" to another file. Assumes "barcodes" is loaded by Step 1a or 1b.
    :param barcodes:
    :param readsfilename: output reads path
    :param answersfilename: output labels path
    :param Q: desired number of simulated reads (10000 is a small value for demo purposes)
    :param srate: substitution rate
    :param irate: insertion rate
    :param drate: deletion rate
    :param nave: average Poisson number of reads from each randomly selected codeword
    :return:
    """
    barcodes_dict = load_barcodes(barcodesfilename)

    q = 0
    reads = []
    answers = []
    while q < Q :
        ansindex = random.randint(low=0,high=barcodes_dict["N"])
        ans = barcodes_dict["allseqs"][ansindex]
        n = min(Q-q, random.poisson(lam=nave))
        for i in range(n) :
            reads.append(makeerrors(copy(ans),srate,irate,drate))
            answers.append(ansindex)
        q += n
    with open(readsfilename, 'w') as OUT:
        for seq in reads :
            OUT.write(decode(seq) + '\n')
    with open(answersfilename, 'w') as OUT:
        for ans in answers :
            OUT.write(str(ans) + '\n')
    print(f"done creating {Q} simulated reads and writing to {readsfilename} \n(answers to {answersfilename})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='generate_reads',
        description='Generates reads from barcodes',
        epilog='TODO')

    parser.add_argument("-b", "--barcodes", required=True)
    parser.add_argument("-r", "--outfile-reads", required=True)
    parser.add_argument("-o", "--outfile-labels", required=True)
    parser.add_argument("-n", "--number", type=int, required=True)
    parser.add_argument("-s", "--p-sub", type=float, required=False, default=0.033)
    parser.add_argument("-i", "--p-ins", type=float, required=False, default=0.033)
    parser.add_argument("-d", "--p-del", type=float, required=False, default=0.033)

    args = parser.parse_args()

    generate_reads(args.barcodes, args.outfile_reads, args.outfile_labels, Q=args.number, srate=args.p_sub, irate=args.p_ins, drate=args.p_del)
